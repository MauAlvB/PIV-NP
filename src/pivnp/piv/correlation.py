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


def one_pass(first: np.ndarray, second: np.ndarray, rows: np.ndarray, cols: np.ndarray,
             window: int, offset_u: np.ndarray | None = None,
             offset_v: np.ndarray | None = None) -> tuple[np.ndarray, ...]:
    """One correlation pass, returning the displacement in pixels and the peak ratio.

    ``offset_u``/``offset_v`` are what a previous pass found, as whole pixels. The result
    includes them, so the caller always gets a total displacement.
    """
    shift_c = None if offset_u is None else np.rint(offset_u).astype(np.int64)
    shift_r = None if offset_v is None else np.rint(offset_v).astype(np.int64)

    base_r, base_c = _corners(rows, cols, window, first.shape)
    moved_r, moved_c = _corners(rows, cols, window, second.shape, shift_r, shift_c)
    a = _windows(first, base_r, base_c, window)
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

    # what is added back is the distance between the two windows that were actually
    # correlated, which near the border is not the offset that was asked for
    u = (peak_col - middle_c) + fine_c + (moved_c - base_c)
    v = (peak_row - middle_r) + fine_r + (moved_r - base_r)
    return u.astype(np.float64), v.astype(np.float64), ratio


def analyse(first: np.ndarray, second: np.ndarray, window: int = 32, overlap: float = 0.5,
            passes: int = 2, mask: np.ndarray | None = None, threshold: float = 2.0,
            region: tuple[int, int, int, int] | None = None,
            smoothing: float = 0.6) -> Field:
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
        u, v, ratio = one_pass(first, second, rows, cols, size, offset_u=u, offset_v=v)
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
