"""The built-in PIV, tested against displacements that are known beforehand.

There is no reference implementation to compare against here, which is the whole difficulty:
PIVlab's answer for a real test is itself a measurement, with its own noise, so it cannot
settle whether a correlation is right to a hundredth of a pixel. What can settle it is a
synthetic pair of images, where the displacement was *put there* and the answer is known
exactly. That is what most of this file does, and it is how the two real errors in the first
version were found:

* the correlation was circular rather than linear and nothing divided out the window
  overlap, which made every displacement come back 8 to 15 % small (0.377 px of error on a
  known shift, now 0.025);
* and the sub-pixel fit returned zero on the first test images, because those were white
  noise -- whose correlation peak is one pixel wide, with nothing for a three-point fit to
  work on. The images here are blurred for that reason, like real speckle.

The agreement with PIVlab on a real test is a separate matter and is not a unit test; it
needs the photographs, which are not in the repository. See ``experiments/end_to_end_piv.py``
and the figures recorded in ``docs/VALIDATION.md``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pivnp.moisture.images import gaussian_blur
from pivnp.piv.correlation import Field, analyse, one_pass, window_centres
from pivnp.piv.settings import DEFAULTS, PivSettingsError, parse
from pivnp.piv.source import ImageSource
from pivnp.piv.validation import find_outliers, replace, smooth
from pivnp.sources import DisplacementSource, available_sources

# A displacement well inside what two passes of a 32 px window can see, with a fractional
# part that is not a half, so a sub-pixel fit that quietly rounds would be caught.
SHIFT_X, SHIFT_Y = 2.3, -1.7


def speckle(shape=(256, 256), sigma=1.5, seed=0) -> np.ndarray:
    """An image that behaves like the surface of soil: random, but not white.

    Blurred noise has a correlation peak a few pixels wide, which is what gives the
    three-point sub-pixel fit something to fit. Real photographs of sand look like this.
    """
    rng = np.random.default_rng(seed)
    rough = rng.normal(128.0, 40.0, shape)
    return np.clip(gaussian_blur(rough, sigma), 0.0, 255.0)


def shifted(image: np.ndarray, dx: float, dy: float) -> np.ndarray:
    """The same image moved by a fractional number of pixels, by linear interpolation."""
    rows, cols = image.shape
    r = np.arange(rows)[:, None] - dy
    c = np.arange(cols)[None, :] - dx
    r0 = np.clip(np.floor(r).astype(int), 0, rows - 2)
    c0 = np.clip(np.floor(c).astype(int), 0, cols - 2)
    fr, fc = r - r0, c - c0
    return ((1 - fr) * ((1 - fc) * image[r0, c0] + fc * image[r0, c0 + 1])
            + fr * ((1 - fc) * image[r0 + 1, c0] + fc * image[r0 + 1, c0 + 1]))


# --- the grid of windows ---------------------------------------------------------------
def test_the_grid_steps_by_the_window_less_the_overlap():
    rows, cols, step = window_centres((128, 256), window=32, overlap=0.5)
    assert step == 16
    assert cols.size == len(range(0, 256 - 32 + 1, 16))
    assert rows.size == len(range(0, 128 - 32 + 1, 16))
    assert cols[0] == pytest.approx(15.5)          # the centre of a 32 px window at 0


def test_no_overlap_gives_windows_side_by_side():
    _, cols, step = window_centres((64, 64), window=32, overlap=0.0)
    assert step == 32
    assert cols.size == 2


def test_a_region_keeps_the_grid_off_the_background():
    _, cols, _ = window_centres((200, 200), window=32, overlap=0.5, region=(50, 50, 150, 150))
    assert cols.min() >= 50 + 15.5
    assert cols.max() <= 150


def test_a_window_too_big_for_the_region_says_so():
    with pytest.raises(ValueError, match="does not fit in the region"):
        window_centres((200, 200), window=64, overlap=0.5, region=(0, 0, 40, 40))


def test_an_impossible_overlap_is_refused():
    with pytest.raises(ValueError, match="overlap must be in"):
        window_centres((64, 64), window=32, overlap=1.0)


# --- the measurement itself ------------------------------------------------------------
def test_a_uniform_shift_is_recovered_to_a_fraction_of_a_pixel():
    """The numbers here were measured first and written down after, not chosen to pass.

    Over six different speckle patterns at this shift: worst single window 0.20 px, 99th
    percentile 0.15, median 0.05, and every window measured. The limits sit just above the
    worst of that, so a change that makes the correlation less accurate fails rather than
    quietly degrades.

    The gap between the median and the 99th percentile is the shape of PIV error and worth
    knowing: nearly every window is accurate to a hundredth of a pixel, and a handful locked
    onto the wrong correlation peak and are wrong by a tenth. Reading the windows between
    pixels improved the median fivefold and left that tail alone -- it cannot rescue a window
    whose first pass went to the wrong place.
    """
    first = speckle()
    second = shifted(first, SHIFT_X, SHIFT_Y)
    field = analyse(first, second, window=32, overlap=0.5, passes=2, smoothing=0.0)

    inside = field.measured
    assert inside.all(), "a clean synthetic pair should leave no window unmeasured"
    error = np.hypot(field.u[inside] - SHIFT_X, field.v[inside] - SHIFT_Y)
    assert np.max(error) < 0.4, "worst window"
    assert np.percentile(error, 99) < 0.2
    assert np.median(error) < 0.08, "the typical window, which is what the result rests on"


def test_the_error_does_not_depend_on_where_the_displacement_falls():
    """Peak locking, and that reading the windows between pixels is what removes it.

    Peak locking is the characteristic fault of correlation PIV: the three-point fit is
    exact when the peak sits on a sample and least certain when it sits between two, so the
    error swings with the *fractional part* of the displacement. An estimator that rounds
    its second-pass offset to a whole pixel has it in full. Measured on this pattern:

    ====== ========= ============
    shift  rounded   interpolated
    ====== ========= ============
    1.00   0.0000    0.0225
    1.25   0.0784    0.0254
    1.50   0.1439    0.0211
    1.75   0.0804    0.0182
    2.00   0.0000    0.0226
    ====== ========= ============

    Rounding is exact at a whole pixel and worst at a half; interpolating is flat. The trade
    is real and worth stating: interpolating is *worse* at an exactly whole displacement,
    0.023 against nothing, and five times better at the half that matters, 0.021 against
    0.144. What is tested is the flatness, because that is the property, and that the worst
    case across the fractions beats the rounded estimator's.
    """
    first = speckle()
    fractions = (1.0, 1.25, 1.5, 1.75, 2.0)
    worst = {}
    for rounded in (False, True):
        errors = []
        for shift in fractions:
            field = analyse(first, shifted(first, shift, 0.0), window=32, overlap=0.5,
                            passes=2, smoothing=0.0, between_pixels=not rounded)
            errors.append(np.max(np.abs(field.u[1:-1, 1:-1] - shift)))
        worst[rounded] = errors

    interpolated = worst[False]
    assert max(interpolated) < 0.04, "interpolating should be accurate at every fraction"
    assert max(interpolated) - min(interpolated) < 0.02, "and flat across them"
    # the rounded estimator swings instead, which is what peak locking looks like
    assert max(worst[True]) > 4 * max(interpolated), "the fault should still be visible"
    assert min(worst[True]) < 1e-9, "rounding is exact at a whole pixel, and only there"


def test_the_displacement_is_not_systematically_small():
    """The bug that mattered: a circular correlation reads every shift 8 to 15 % low.

    This is the test that would have caught it. It is stated as a *bias* -- the mean of the
    signed error -- because the fault was not scatter: every window was wrong in the same
    direction. And it is stated in pixels rather than as a percentage, because the bias of a
    correlation does not scale with the displacement; 0.02 px is 1 % of a 2 px shift and
    0.1 % of a 20 px one, and the pixel is the honest number. Measured: 0.0040 px at worst
    over six speckle patterns. Reading the windows between pixels brings it to 0.0040.
    """
    first = speckle()
    second = shifted(first, SHIFT_X, SHIFT_Y)
    field = analyse(first, second, window=32, overlap=0.5, passes=2, smoothing=0.0)
    inside = field.measured

    bias_x = np.mean(field.u[inside] - SHIFT_X)
    bias_y = np.mean(field.v[inside] - SHIFT_Y)
    assert abs(bias_x) < 0.03, f"x read {bias_x:+.4f} px off, in one direction"
    assert abs(bias_y) < 0.03, f"y read {bias_y:+.4f} px off, in one direction"


def test_a_window_pushed_past_the_edge_still_reports_the_right_displacement():
    """The bug this is named for, and how it hid.

    The second pass takes each window of the second image from where the first pass said the
    soil went. At the top of the image that offset can point outside the photograph, and the
    window is pulled back in -- so the shift the code asked for is not the shift it got. The
    first version added back the shift it *asked for*, which made the whole top row of the
    grid wrong by the full offset: a true -1.7 px came back as -3.7.

    It stayed hidden because the vectors were so wrong that the outlier test threw the row
    away, and a thrown-away row looks like the edge of the soil rather than like a mistake.
    So the check is on the measurement before any validation can bury it: the top row has to
    agree with the rest of the field.
    """
    first = speckle()
    second = shifted(first, 0.0, -1.7)        # upwards, so the top row's offset leaves the image
    field = analyse(first, second, window=32, overlap=0.5, passes=2, smoothing=0.0,
                    threshold=1e9)            # no rejection: nothing may hide the error
    np.testing.assert_allclose(field.v[0], -1.7, atol=0.2)
    np.testing.assert_allclose(field.v[-1], -1.7, atol=0.2)


def test_images_that_did_not_move_measure_no_movement():
    first = speckle()
    field = analyse(first, first.copy(), window=32, overlap=0.5, passes=2, smoothing=0.0)
    inside = field.measured
    assert np.nanmax(np.abs(field.u[inside])) < 1e-6
    assert np.nanmax(np.abs(field.v[inside])) < 1e-6


def test_a_shear_is_recovered_where_it_was_put():
    """A displacement that varies across the image, which is the case that matters.

    The first and last column of the grid are left out, and not to make the test pass:
    :func:`shifted` builds the second image by interpolation and has to clamp at the border,
    so the leftmost and rightmost few pixels of it are not a shifted copy of anything. Those
    columns are a flaw in the *image*, measured at 1.0 px of apparent error, and including
    them would be testing the test. Everywhere else the error is 0.14 px, against a
    displacement that changes by half a pixel across one window.
    """
    first = speckle(shape=(256, 256), seed=3)
    # u grows from 0 at the top to 4 px at the bottom; v stays zero
    second = np.empty_like(first)
    for row in range(256):
        second[row] = shifted(first, 4.0 * row / 255.0, 0.0)[row]

    field = analyse(first, second, window=32, overlap=0.5, passes=2, smoothing=0.0)
    expected = 4.0 * field.y / 255.0
    inside = field.measured
    inside[:, 0] = inside[:, -1] = False

    assert np.nanmax(np.abs(field.u[inside] - expected[inside])) < 0.2
    assert np.nanmax(np.abs(field.v[inside])) < 0.2
    # and the gradient itself, which is what the strain is made of, survives.
    # The slope is written out rather than fitted with polyfit, which goes through a matrix
    # solve: this environment loads several BLAS libraries at once and any of them crashes
    # the interpreter. Recorded in docs/STATUS.md.
    y, u = field.y[inside], field.u[inside]
    dy, du = y - y.mean(), u - u.mean()
    slope = (dy * du).sum() / (dy * dy).sum()
    assert slope == pytest.approx(4.0 / 255.0, rel=0.05)


def test_one_pass_can_be_told_where_to_look():
    """The offset is how the second pass works, and it has to come back in the total."""
    first = speckle()
    second = shifted(first, 6.0, 0.0)
    rows, cols, _ = window_centres(first.shape, 32, 0.5)
    corner_r = np.rint(rows - 15.5).astype(np.int64)
    corner_c = np.rint(cols - 15.5).astype(np.int64)

    offset = np.full((corner_r.size, corner_c.size), 5.0)
    u, v, _ = one_pass(first, second, corner_r, corner_c, 32,
                       offset_u=offset, offset_v=np.zeros_like(offset))
    # 5 px of offset plus the 1 px left to measure. The last column is left out: its offset
    # points past the right edge, so there the window is pulled back in -- the case the
    # window-clipping test covers.
    assert np.median(u[:, :-1]) == pytest.approx(6.0, abs=0.1)


def test_the_mask_keeps_the_measurement_off_the_background():
    """What the mask has to do is reject the background, and only the background.

    Stated in both directions on purpose. The first version of this test asked that
    everything inside the mask be measured, which is the wrong question -- the outlier test
    also rejects, for its own reasons -- and the vagueness is what let the window-clipping
    bug sit here unnoticed.
    """
    first = speckle()
    second = shifted(first, SHIFT_X, SHIFT_Y)
    mask = np.zeros(first.shape, dtype=bool)
    mask[:128] = True                                     # material only in the top half

    field = analyse(first, second, window=32, overlap=0.5, passes=2, mask=mask,
                    smoothing=0.0)
    assert not field.measured[field.y > 150].any(), "measured where there is no material"
    assert field.measured[field.y < 110].mean() > 0.9, "lost material that was there"


def test_images_of_different_sizes_are_refused():
    with pytest.raises(ValueError, match="differ in size"):
        analyse(speckle((64, 64)), speckle((64, 128)))


def test_at_least_one_pass():
    with pytest.raises(ValueError, match="passes must be at least 1"):
        analyse(speckle((64, 64)), speckle((64, 64)), passes=0)


def test_the_field_reports_where_it_measured():
    first = speckle((96, 96))
    field = analyse(first, shifted(first, 1.0, 0.0), window=32, overlap=0.5, passes=1,
                    smoothing=0.0)
    assert isinstance(field, Field)
    assert field.shape == field.u.shape == field.peak_ratio.shape
    assert field.measured.dtype == bool


# --- validation and smoothing ----------------------------------------------------------
def test_a_single_wrong_vector_is_found():
    u = np.full((9, 9), 1.0)
    v = np.full((9, 9), 0.5)
    u[4, 4] = 20.0                                        # a correlation on the wrong peak
    bad = find_outliers(u, v, threshold=2.0)
    assert bad[4, 4]
    assert bad.sum() == 1


def test_a_shear_band_is_not_mistaken_for_an_error():
    """The normalised median test is used precisely because a vector may differ for a reason."""
    u = np.zeros((11, 11))
    u[6:] = 3.0                                           # everything below the band moves
    bad = find_outliers(u, np.zeros((11, 11)), threshold=2.0)
    assert not bad.any()


def test_replacing_fills_from_the_neighbours():
    field = np.full((7, 7), 2.0)
    bad = np.zeros((7, 7), dtype=bool)
    bad[3, 3] = True
    out = replace(field, bad)
    assert out[3, 3] == pytest.approx(2.0)
    assert np.isfinite(out).all()


def test_smoothing_leaves_a_straight_field_alone():
    """A field with constant strain must keep that strain exactly, edges included.

    Smoothing is there to take out scatter. If it also flattened a real gradient it would
    take the strain with it, so the case to pin down is the one where the answer is known: a
    linear field has the same gradient everywhere, and smoothing a straight line must give
    back the same straight line.

    The first version of this test only looked at ``[2:-2, 2:-2]`` -- and passed, over an
    implementation that was flattening the gradient at every edge. The edge is the whole
    point, so the whole field is checked.
    """
    rows, cols = np.mgrid[0:15, 0:15]
    field = 0.3 + 0.7 * cols - 0.2 * rows
    np.testing.assert_allclose(smooth(field, 0.6), field, atol=1e-8)


def test_smoothing_keeps_the_gradient_at_a_hole_in_the_field():
    """The real case misses about 42 % of its points, so most edges are interior ones.

    A straight-line continuation past each hole is what keeps the gradient; averaging over
    whichever neighbours happen to exist does not. Measured over ten hole patterns on this
    field, the two ways of doing it come out:

    ====================  ==============  =============
    ..                    typical point   worst point
    ====================  ==============  =============
    continued outwards    0.012           0.147
    renormalised kernel   0.060           0.196
    ====================  ==============  =============

    which is five times better where it counts -- the typical point, since the worst one is
    a single place where a hole leaves almost nothing to extrapolate from. The limits below
    are set from the measurement, and the typical-point one is what would catch a quiet
    return to renormalising.
    """
    rng = np.random.default_rng(0)
    rows, cols = np.mgrid[0:20, 0:20]
    field = 0.7 * cols - 0.2 * rows
    holed = np.where(rng.random(field.shape) < 0.42, np.nan, field)

    kept = np.isfinite(holed)
    off = np.abs(smooth(holed, 0.6)[kept] - field[kept])
    assert np.median(off) < 0.02, f"the typical point was bent by {np.median(off):.4f}"
    assert np.max(off) < 0.2, f"the worst point was bent by {np.max(off):.4f}"


def test_smoothing_takes_out_scatter():
    rng = np.random.default_rng(1)
    truth = np.zeros((40, 40))
    noisy = truth + rng.normal(0.0, 0.1, truth.shape)
    out = smooth(noisy, 0.6)
    # what matters is the difference between neighbours, which is what the strain is made of
    before = np.median(np.abs(np.diff(noisy, axis=1)))
    after = np.median(np.abs(np.diff(out, axis=1)))
    assert after < 0.5 * before


def test_smoothing_does_not_invent_a_measurement():
    field = np.full((9, 9), 1.0)
    field[4, 4] = np.nan
    out = smooth(field, 0.6)
    assert np.isnan(out[4, 4])
    assert np.isfinite(out[np.isfinite(field)]).all()


def test_smoothing_near_a_hole_is_not_pulled_towards_zero():
    """The normalised convolution: a missing neighbour must be ignored, not counted as 0."""
    field = np.full((9, 9), 5.0)
    field[0, :] = np.nan                                   # the soil ends here
    out = smooth(field, 0.8)
    assert out[1, 4] == pytest.approx(5.0, abs=1e-9)


def test_no_smoothing_changes_nothing():
    field = np.array([[1.0, 5.0], [9.0, 2.0]])
    np.testing.assert_array_equal(smooth(field, 0.0), field)


# --- the settings file -----------------------------------------------------------------
MINIMAL = "IMAGES = shot_{n}.jpg\nSCALE = 0.0005\n"


def test_the_minimum_is_the_images_and_the_scale():
    settings = parse(MINIMAL)
    assert settings.image_pattern == "shot_{n}.jpg"
    assert settings.scale == 0.0005
    assert settings.window == 32                           # the rest comes from the defaults
    assert settings.passes == 2
    assert settings.smoothing == 0.6


def test_zero_padded_numbering_is_accepted():
    """``Masked_001.jpg`` is how the dam-break photographs are named, and it must work."""
    settings = parse("IMAGES = Masked_{n:03d}.jpg\nSCALE = 0.0005\n")
    assert settings.image_path(Path("d"), 1).name == "Masked_001.jpg"
    assert settings.image_path(Path("d"), 12).name == "Masked_012.jpg"


def test_a_pattern_without_a_number_is_refused():
    with pytest.raises(PivSettingsError, match="IMAGES must contain"):
        parse("IMAGES = only_one.jpg\nSCALE = 0.0005\n")


def test_the_first_image_can_be_numbered_anything():
    settings = parse(MINIMAL + "FIRST_IMAGE = 40\n")
    assert settings.image_path(Path("d"), 1).name == "shot_40.jpg"
    assert settings.image_path(Path("d"), 3).name == "shot_42.jpg"


def test_a_missing_scale_says_why_it_is_needed():
    with pytest.raises(PivSettingsError, match="metres per pixel"):
        parse("IMAGES = a_{n}.jpg\n")


def test_comments_and_order_do_not_matter():
    settings = parse("! a comment\nSCALE = 0.001  ! after the value\n"
                     "WINDOW = 16\nIMAGES = x_{n}.png\n")
    assert (settings.window, settings.scale) == (16, 0.001)


def test_a_line_without_a_value_is_refused():
    with pytest.raises(PivSettingsError, match="expected 'KEY = value'"):
        parse("IMAGES\n")


def test_a_window_that_is_not_a_power_of_two_is_refused():
    with pytest.raises(PivSettingsError, match="power of two"):
        parse(MINIMAL + "WINDOW = 24\n")


def test_a_window_too_small_to_correlate_is_refused():
    with pytest.raises(PivSettingsError, match="too small to correlate"):
        parse(MINIMAL + "WINDOW = 4\n")


def test_the_region_is_four_numbers():
    settings = parse(MINIMAL + "REGION = 10, 20, 300, 400\n")
    assert settings.region == (10, 20, 300, 400)
    with pytest.raises(PivSettingsError, match="REGION needs four numbers"):
        parse(MINIMAL + "REGION = 10, 20\n")


def test_a_region_smaller_than_a_window_is_refused():
    with pytest.raises(PivSettingsError, match="smaller than one window"):
        parse(MINIMAL + "WINDOW = 32\nREGION = 0, 0, 20, 20\n")


def test_smoothing_can_be_turned_off_but_not_set_absurdly():
    assert parse(MINIMAL + "SMOOTH = 0\n").smoothing == 0.0
    with pytest.raises(PivSettingsError, match="cannot be negative"):
        parse(MINIMAL + "SMOOTH = -1\n")
    with pytest.raises(PivSettingsError, match="flatten the real strain"):
        parse(MINIMAL + "SMOOTH = 8\n")


def test_two_ways_of_masking_at_once_is_refused():
    with pytest.raises(PivSettingsError, match="give only one"):
        parse(MINIMAL + "MASK_BELOW = 20\nMASK_IMAGE = mask.png\n")


def test_an_unrecognized_key_is_reported_not_silently_dropped():
    settings = parse(MINIMAL + "WIDNOW = 16\n")
    assert settings.unknown_keys == ("WIDNOW",)


def test_a_value_that_is_not_a_number_says_which_key():
    with pytest.raises(PivSettingsError, match="WINDOW must be a number"):
        parse(MINIMAL + "WINDOW = large\n")


def test_reading_between_pixels_is_off_unless_asked_for():
    """Off by default, because of what it does to the accumulated answer; see the module."""
    assert parse(MINIMAL).between_pixels is False
    assert parse(MINIMAL + "SUBPIXEL_OFFSET = 1\n").between_pixels is True
    assert parse(MINIMAL + "SUBPIXEL_OFFSET = 0\n").between_pixels is False


def test_every_key_has_a_default_or_is_required():
    assert set(DEFAULTS) >= {"IMAGES", "SCALE", "WINDOW", "OVERLAP", "PASSES", "SMOOTH"}


# --- the source ------------------------------------------------------------------------
def test_the_image_source_is_available_by_name():
    assert "images" in available_sources()


def test_the_image_source_satisfies_the_contract(workdir: Path):
    settings = parse(MINIMAL)
    source = ImageSource(workdir, settings, config=None)
    assert isinstance(source, DisplacementSource)
    assert source.frame_interval() is None                 # photographs do not carry a time
    assert source.moisture is False
