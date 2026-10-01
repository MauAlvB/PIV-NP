"""Particle image velocimetry: displacement between two photographs, by cross-correlation.

This is what PIVlab does, done inside PIV-NP so that a case can be analysed from the
photographs alone. It is deliberately a *basic* PIV -- two passes and a sub-pixel fit -- not
a replacement for a dedicated package. The point is that someone with a sequence of images
and no PIV experience can get a result, look at it, and decide later whether to learn PIVlab
and compare.

How it works, in the order it happens:

1. The image is cut into square interrogation windows on a regular grid, overlapping.
2. Each window of the first image is cross-correlated with the same window of the second.
   The correlation is done by FFT, for **every window at once** rather than one at a time:
   a 1920x1080 image gives about eight thousand windows per pass, and a Python loop over
   them would make the whole thing unusable.
3. The peak of each correlation says how far that patch of soil moved, to the nearest pixel.
   A Gaussian fit through the peak and its two neighbours refines it to a fraction of one.
4. A second pass repeats it with smaller windows, offsetting the second image by what the
   first pass found. This is where most of the accuracy comes from: the remaining
   displacement is near zero, which is where the correlation is sharpest and the sub-pixel
   estimator is least biased.

Between passes the field is validated and repaired, because one wrong vector used as the
offset for the next pass does not stay one wrong vector.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Field:
    """Displacement measured on a grid of interrogation windows, in pixels.

    ``x`` and ``y`` are the centres of the windows in the image, with ``y`` measured
    downwards as images are. ``u`` and ``v`` are the displacement in those same axes, and
    ``NaN`` where nothing could be measured -- outside the mask, or rejected as an outlier.
    """

    x: np.ndarray          # (rows, cols) column of the window centre, in pixels
    y: np.ndarray          # (rows, cols) row of the window centre
    u: np.ndarray          # (rows, cols) displacement along x, in pixels
    v: np.ndarray          # (rows, cols) displacement along y, in pixels
    peak_ratio: np.ndarray  # (rows, cols) how far the peak stands above the rest

    @property
    def shape(self) -> tuple[int, int]:
        return self.x.shape

    @property
    def measured(self) -> np.ndarray:
        return np.isfinite(self.u) & np.isfinite(self.v)


def window_centres(shape: tuple[int, int], window: int, overlap: float,
                   region: tuple[int, int, int, int] | None = None,
                   ) -> tuple[np.ndarray, np.ndarray, int]:
    """Where the interrogation windows sit, and the step between them.

    ``region`` is ``left, top, right, bottom`` in pixels and limits the grid to the part of
    the photograph worth measuring; without it the grid covers everything, including
    background that will never move.
    """
    if not 0.0 <= overlap < 1.0:
        raise ValueError(f"overlap must be in [0, 1) and it is {overlap}")
    step = max(1, int(round(window * (1.0 - overlap))))
    left, top, right, bottom = region or (0, 0, shape[1], shape[0])
    right, bottom = min(right, shape[1]), min(bottom, shape[0])
    rows = np.arange(max(0, top), bottom - window + 1, step)
    cols = np.arange(max(0, left), right - window + 1, step)
    if rows.size == 0 or cols.size == 0:
        raise ValueError(f"a window of {window} px does not fit in the region "
                         f"{left},{top} to {right},{bottom} of a "
                         f"{shape[1]}x{shape[0]} image")
    centre = (window - 1) / 2.0
    return rows + centre, cols + centre, step


def _corners(rows: np.ndarray, cols: np.ndarray, window: int, shape: tuple[int, int],
             offset_r: np.ndarray | None = None,
             offset_c: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Top-left corner of every window, kept inside the image.

    The offsets are what the second pass needs: each window of the second image is taken
    from where the first pass said the soil went, so that what is left to measure is close
    to zero. Near the border of the photograph an offset can point outside it, and then the
    window is pulled back in -- so the offset that is *applied* is not always the one asked
    for, and the caller has to work out the displacement from the corners that came back
    rather than from what it requested. Getting that wrong reads a displacement wrong by the
    whole offset; see the test named after it.
    """
    r0 = rows[:, None] + (offset_r if offset_r is not None else 0)
    c0 = cols[None, :] + (offset_c if offset_c is not None else 0)
    r0 = np.clip(np.broadcast_to(r0, (rows.size, cols.size)), 0, shape[0] - window)
    c0 = np.clip(np.broadcast_to(c0, (rows.size, cols.size)), 0, shape[1] - window)
    return r0, c0


def _windows(image: np.ndarray, corner_r: np.ndarray, corner_c: np.ndarray,
             window: int) -> np.ndarray:
    """Stack the interrogation window at every corner into one array."""
    inside = np.arange(window)
    r = corner_r[:, :, None, None] + inside[None, None, :, None]   # (rows, cols, window, 1)
    c = corner_c[:, :, None, None] + inside[None, None, None, :]   # (rows, cols, 1, window)
    return image[r, c]                                             # (rows, cols, w, w)


def _windows_between_pixels(image: np.ndarray, corner_r: np.ndarray, corner_c: np.ndarray,
                            window: int, fraction_r: np.ndarray,
                            fraction_c: np.ndarray) -> np.ndarray:
    """The window at each corner, taken a fraction of a pixel further along.

    This is what lets the second pass finish the job. The offset it gets from the first pass
    is a real number, and taking only its whole part leaves the fractional part still to be
    measured -- on a slow test that *is* the whole displacement, so the pass contributes
    nothing at all. Measured on a known shift of 0.16 px: one pass and two gave identical
    numbers until this existed.

    Each pixel of the window is read between the four pixels of the photograph around it,
    weighted by how close it falls to each. The four readings are accumulated one at a time
    rather than built and then summed, which keeps the memory to two of these stacks instead
    of five; on a 1920x1080 photograph one stack is already hundreds of megabytes.

    ``fraction_r`` and ``fraction_c`` are in ``[0, 1)``, so the four pixels are always the
    corner and its neighbours below and to the right.

    The interpolation is **cubic**, over four pixels along each axis, and that is not
    refinement for its own sake. Reading between pixels with a straight line average does to
    the photograph what a blur does: it takes out the fine detail the correlation lives on.
    The residual then comes back smaller than it is, so the offset that was already applied
    is confirmed rather than corrected, and whatever the first pass under-read stays
    under-read. Measured on a shear whose answer is known, the straight-line version turned
    an error of -1.1 % into -8.6 % while improving everything else -- it is the reason this
    is cubic and not two lines shorter.

    The indices are kept inside the photograph, which near an edge means a tap falls on the
    edge pixel twice. Reserving room by pulling the *window* in instead costs the last row
    and column of the grid a whole pixel of offset they never asked for, which then has to be
    measured and comes back with the 0.12 px that a one-pixel shift costs; the first version
    did that and the test on two identical photographs caught it.
    """
    inside = np.arange(window)
    r = corner_r[:, :, None, None] + inside[None, None, :, None]
    c = corner_c[:, :, None, None] + inside[None, None, None, :]
    rows, cols = image.shape

    weights_r = _cubic_weights(fraction_r)
    weights_c = _cubic_weights(fraction_c)
    # the sixteen terms are accumulated into one array, one gather at a time. A separable
    # version with a buffer per axis is tidier and needs a third array of this size, which
    # on a real photograph is another hundred megabytes for nothing.
    out = np.zeros(r.shape[:2] + (window, window))
    for i, wr in enumerate(weights_r):
        taken_r = np.clip(r + i - 1, 0, rows - 1)
        for j, wc in enumerate(weights_c):
            taken_c = np.clip(c + j - 1, 0, cols - 1)
            out += (wr * wc)[:, :, None, None] * image[taken_r, taken_c]
    return out


def _cubic_weights(fraction: np.ndarray) -> tuple[np.ndarray, ...]:
    """Catmull-Rom weights for the four pixels around a fractional position.

    The standard cubic convolution kernel, which passes through the samples and keeps the
    gradient continuous across them. The four weights sum to one at every fraction, so a
    uniform patch of photograph comes back unchanged.
    """
    t = fraction
    t2 = t * t
    t3 = t2 * t
    return (-0.5 * t3 + t2 - 0.5 * t,
            1.5 * t3 - 2.5 * t2 + 1.0,
            -1.5 * t3 + 2.0 * t2 + 0.5 * t,
            0.5 * t3 - 0.5 * t2)


def _overlap_weight(window: int) -> np.ndarray:
    """How many pixels two windows still share at each displacement.

    Two windows shifted by ``d`` only overlap on ``window - |d|`` rows or columns, so the
    correlation falls off towards the edges for a reason that has nothing to do with the
    soil. Dividing it out is what stops every displacement coming back smaller than it is.
    """
    shifts = np.arange(-window, window)
    along = np.maximum(window - np.abs(shifts), 0).astype(np.float64)
    return np.maximum(np.outer(along, along), 1.0)


def _correlate(first: np.ndarray, second: np.ndarray) -> np.ndarray:
    """Cross-correlation of every pair of windows at once, peak centred on zero displacement.

    Two things here are what separate a PIV that measures from one that nearly measures:

    * each window has its **mean removed**, or a window that is simply brighter than its
      neighbour correlates with everything;
    * the transform is taken over **twice the window**, so the correlation is linear rather
      than circular. Without that the window wraps around on itself and a displacement near
      the edge of the search range finds a peak that is not there. The price is four times
      the memory, which is why the caller works in chunks.

    What comes back is still weighted by how much the two windows overlap at each shift; the
    caller divides that out.
    """
    window = first.shape[-1]
    a = first - first.mean(axis=(-2, -1), keepdims=True)
    b = second - second.mean(axis=(-2, -1), keepdims=True)
    size = (2 * window, 2 * window)
    spectrum = (np.fft.rfft2(b, s=size, axes=(-2, -1))
                * np.conj(np.fft.rfft2(a, s=size, axes=(-2, -1))))
    power = np.fft.irfft2(spectrum, s=size, axes=(-2, -1))
    return np.fft.fftshift(power, axes=(-2, -1)) / _overlap_weight(window)


def _peak(power: np.ndarray) -> tuple[np.ndarray, ...]:
    """Integer peak of each correlation plane, and how far it stands above the rest."""
    rows, cols = power.shape[-2:]
    flat = power.reshape(*power.shape[:-2], rows * cols)
    index = np.argmax(flat, axis=-1)
    peak_row, peak_col = np.divmod(index, cols)
    highest = np.take_along_axis(flat, index[..., None], axis=-1)[..., 0]

    # the second highest, ignoring a 5x5 patch around the peak, is what the peak is judged
    # against: a peak that barely stands out is one the displacement can hide behind
    r = np.arange(rows).reshape(1, 1, rows, 1)
    c = np.arange(cols).reshape(1, 1, 1, cols)
    near = ((np.abs(r - peak_row[..., None, None]) <= 2)
            & (np.abs(c - peak_col[..., None, None]) <= 2))
    runner_up = np.where(near, -np.inf, power).reshape(*power.shape[:-2], rows * cols)
    second = runner_up.max(axis=-1)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(second > 0, highest / second, np.inf)
    return peak_row, peak_col, ratio


def _neighbourhood(power: np.ndarray, peak_row: np.ndarray,
                   peak_col: np.ndarray) -> np.ndarray:
    """The 3x3 of correlation values around each window's peak, as ``(rows, cols, 3, 3)``."""
    height, width = power.shape[-2:]
    around = np.array([-1, 0, 1])
    r = np.clip(peak_row[..., None] + around, 0, height - 1)     # (rows, cols, 3)
    c = np.clip(peak_col[..., None] + around, 0, width - 1)
    i = np.arange(power.shape[0])[:, None, None, None]
    j = np.arange(power.shape[1])[None, :, None, None]
    return power[i, j, r[:, :, :, None], c[:, :, None, :]]


def _gaussian_shift(low: np.ndarray, middle: np.ndarray, high: np.ndarray) -> np.ndarray:
    """Where the peak really is, between the sample at ``middle`` and its two neighbours.

    Fitting a parabola to the logarithm of three samples is exact when the correlation peak
    is Gaussian, which it nearly is -- this is the estimator every PIV package uses. Where a
    sample is not positive the logarithm has no meaning, and the integer peak is kept.
    """
    usable = (low > 0) & (middle > 0) & (high > 0)
    with np.errstate(divide="ignore", invalid="ignore"):
        a = np.log(np.where(usable, low, 1.0))
        b = np.log(np.where(usable, middle, 1.0))
        c = np.log(np.where(usable, high, 1.0))
        bottom = 2.0 * (a - 2.0 * b + c)
        shift = np.where(np.abs(bottom) > 1e-12, (a - c) / bottom, 0.0)
    return np.where(usable & np.isfinite(shift), np.clip(shift, -1.0, 1.0), 0.0)


def _sub_pixel(power: np.ndarray, peak_row: np.ndarray,
               peak_col: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Fractional part of the displacement, along each axis."""
    height, width = power.shape[-2:]
    patch = _neighbourhood(power, peak_row, peak_col)
    middle = patch[..., 1, 1]

    along_x = _gaussian_shift(patch[..., 1, 0], middle, patch[..., 1, 2])
    along_y = _gaussian_shift(patch[..., 0, 1], middle, patch[..., 2, 1])

    # a peak sitting on the edge of the plane has no neighbour on one side
    inside_x = (peak_col > 0) & (peak_col < width - 1)
    inside_y = (peak_row > 0) & (peak_row < height - 1)
    return np.where(inside_x, along_x, 0.0), np.where(inside_y, along_y, 0.0)


def _split(offset: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """An offset as a whole number of pixels and a fraction in ``[0, 1)``.

    A fraction that is a hair short of a whole pixel is counted as the whole pixel. Without
    that, an offset of -1e-16 -- which is what two identical photographs produce, from
    round-off in the transform -- comes out as one pixel back plus a fraction of 0.999...,
    and the window is then resampled where there was nothing to resample. The arithmetic
    still lands in the right place, but the interpolation no longer returns the pixel
    exactly, and two identical photographs stop measuring exactly no movement.
    """
    whole = np.floor(offset)
    fraction = offset - whole
    at_the_top = fraction > 1.0 - 1e-9
    whole = np.where(at_the_top, whole + 1.0, whole)
    fraction = np.where(at_the_top | (fraction < 1e-9), 0.0, fraction)
    return whole.astype(np.int64), fraction


#: How much memory one pass may use for its correlation planes, in bytes. The planes are
#: four times the window area per window, and a photograph of a few megapixels with a large
#: window asks for gigabytes of them at once -- enough to bring the interpreter down rather
#: than merely slow it. Splitting the grid into bands of rows costs nothing measurable: the
#: transform is already done for thousands of windows at a time within each band.
CORRELATION_BUDGET = 192 * 1024 * 1024


def _band_height(columns: int, window: int) -> int:
    """How many rows of the grid fit in one band, given the budget."""
    per_row = columns * (2 * window) ** 2 * 8      # the correlation plane, in bytes
    return max(1, int(CORRELATION_BUDGET // max(per_row, 1)))


def one_pass(first: np.ndarray, second: np.ndarray, rows: np.ndarray, cols: np.ndarray,
             window: int, offset_u: np.ndarray | None = None,
             offset_v: np.ndarray | None = None,
             between_pixels: bool = False) -> tuple[np.ndarray, ...]:
    """One correlation pass, returning the displacement in pixels and the peak ratio.

    ``offset_u``/``offset_v`` are what a previous pass found. The result includes whatever
    part of them was applied, so the caller always gets a total displacement.

    With ``between_pixels`` the offset is applied in full, the fractional part included, by
    reading each window between the pixels of the photograph. What is then left to measure is
    near zero, which is where the correlation peak is sharpest and the three-point fit least
    biased. Without it only the whole pixels are applied, which is what the first version
    did and which leaves a slow test no better off for the second pass at all.

    The grid is worked through in bands of rows, so that the correlation planes of a large
    window on a large photograph do not all have to exist at once.
    """
    height = _band_height(cols.size, window)
    if height < rows.size:
        pieces = []
        for start in range(0, rows.size, height):
            stop = min(start + height, rows.size)
            pieces.append(one_pass(
                first, second, rows[start:stop], cols, window,
                offset_u=None if offset_u is None else offset_u[start:stop],
                offset_v=None if offset_v is None else offset_v[start:stop],
                between_pixels=between_pixels))
        return tuple(np.concatenate([piece[which] for piece in pieces], axis=0)
                     for which in range(3))

    whole_c = whole_r = None
    extra_c = extra_r = 0.0
    if offset_u is not None:
        if between_pixels:
            # floor rather than rint, so the fraction is never negative and the four pixels
            # read are always the corner and its neighbours below and to the right
            whole_c, extra_c = _split(offset_u)
            whole_r, extra_r = _split(offset_v)
        else:
            whole_c = np.rint(offset_u).astype(np.int64)
            whole_r = np.rint(offset_v).astype(np.int64)

    base_r, base_c = _corners(rows, cols, window, first.shape)
    a = _windows(first, base_r, base_c, window)
    if whole_c is None:
        moved_r, moved_c = _corners(rows, cols, window, second.shape)
        b = _windows(second, moved_r, moved_c, window)
    elif between_pixels:
        moved_r, moved_c = _corners(rows, cols, window, second.shape, whole_r, whole_c)
        b = _windows_between_pixels(second, moved_r, moved_c, window, extra_r, extra_c)
    else:
        moved_r, moved_c = _corners(rows, cols, window, second.shape, whole_r, whole_c)
        b = _windows(second, moved_r, moved_c, window)
    power = _correlate(a, b)

    # Only the middle of the plane is searched. Past half a window the two windows barely
    # overlap, so what is there is noise divided by a small number: the classic rule that a
    # displacement should stay under a quarter of the window is the same statement.
    middle_r, middle_c = power.shape[-2] // 2, power.shape[-1] // 2
    reach = max(2, window // 2)
    searched = np.full(power.shape[-2:], -np.inf)
    searched[middle_r - reach:middle_r + reach + 1,
             middle_c - reach:middle_c + reach + 1] = 0.0
    power = power + searched

    peak_row, peak_col, ratio = _peak(power)
    fine_c, fine_r = _sub_pixel(power, peak_row, peak_col)

    # What is added back is the distance between the two windows that were actually
    # correlated, which near the border is not the offset that was asked for: the window
    # gets pulled back inside the photograph and the whole-pixel part of the offset is lost
    # with it. The fraction survives, because it is applied by interpolation and not by
    # where the window was taken from.
    u = (peak_col - middle_c) + fine_c + (moved_c - base_c) + extra_c
    v = (peak_row - middle_r) + fine_r + (moved_r - base_r) + extra_r
    return u.astype(np.float64), v.astype(np.float64), ratio


def analyse(first: np.ndarray, second: np.ndarray, window: int = 32, overlap: float = 0.5,
            passes: int = 2, mask: np.ndarray | None = None, threshold: float = 2.0,
            region: tuple[int, int, int, int] | None = None,
            smoothing: float = 0.6, between_pixels: bool = False) -> Field:
    """Displacement between two photographs, in pixels, on a grid of windows.

    The first pass uses a window twice as large for every pass that follows, which is what
    lets it find a displacement it knows nothing about; each pass after that offsets the
    second image by what the last one found, so the remaining displacement is near zero --
    where the correlation is sharpest and the sub-pixel estimator least biased.

    ``mask`` is True where there is material to measure. A window is kept when most of it is
    inside; the rest come back as ``NaN``, which is what the analysis expects to see at the
    edge of the soil.

    ``smoothing`` is the width, in grid points, of the Gaussian applied to the finished
    field; see :func:`~pivnp.piv.validation.smooth` for why it is on by default and set to
    this value. Zero returns the raw correlation result.

    ``between_pixels`` reads each window between the pixels of the photograph, so the offset
    a pass inherits is applied in full instead of rounded to a whole pixel. On the field it
    produces, it is better by a wide margin: four to five times less bias and scatter on
    known shifts, and it removes peak locking, the error that otherwise swings with where the
    displacement happens to fall between two pixels.

    **It is off by default anyway**, and the reason is the one measurement that matters more
    than those. On the example whose accumulated answer is known, it makes the accumulated
    shear *worse* -- 8.4 % low against 1.1 % -- while making the artifacts that should be
    zero about five times smaller. Both estimators read the gradient of each single step to
    within 1.3 %, so the loss is somewhere in ten steps of accumulation and is not
    understood. Until it is, the better per-step field is not worth a worse answer.
    """
    if passes < 1:
        raise ValueError(f"passes must be at least 1 and it is {passes}")
    if first.shape != second.shape:
        raise ValueError(f"the two images differ in size: {first.shape} and {second.shape}")

    centres_y, centres_x, _ = window_centres(first.shape, window, overlap, region)
    corner_r = np.rint(centres_y - (window - 1) / 2.0).astype(np.int64)
    corner_c = np.rint(centres_x - (window - 1) / 2.0).astype(np.int64)

    u = v = ratio = None
    for step in range(passes):
        size = window * (2 ** (passes - 1 - step))
        # every pass sits on the same grid, so a window of a different size is centred on it
        grow = (size - window) // 2
        rows = np.clip(corner_r - grow, 0, first.shape[0] - size)
        cols = np.clip(corner_c - grow, 0, first.shape[1] - size)
        u, v, ratio = one_pass(first, second, rows, cols, size, offset_u=u, offset_v=v,
                               between_pixels=between_pixels)
        if step < passes - 1:
            from .validation import find_outliers, replace
            bad = find_outliers(u, v, threshold)
            u, v = replace(u, bad), replace(v, bad)

    from .validation import find_outliers, smooth
    outside = _mask_rejects(mask, corner_r, corner_c, window)
    bad = find_outliers(u, v, threshold) | outside
    u = smooth(np.where(bad, np.nan, u), smoothing)
    v = smooth(np.where(bad, np.nan, v), smoothing)
    x, y = np.meshgrid(centres_x, centres_y)
    return Field(x=x, y=y, u=u, v=v, peak_ratio=np.where(bad, np.nan, ratio))


def _mask_rejects(mask: np.ndarray | None, rows: np.ndarray, cols: np.ndarray,
                  window: int) -> np.ndarray:
    """Windows with too little material in them to be worth correlating."""
    if mask is None:
        return np.zeros((rows.size, cols.size), dtype=bool)
    corner_r, corner_c = _corners(rows, cols, window, mask.shape)
    inside = _windows(mask.astype(np.float64), corner_r, corner_c, window)
    return inside.mean(axis=(-2, -1)) < 0.5
