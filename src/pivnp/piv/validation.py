"""Finding and replacing the vectors that are simply wrong, and smoothing the small scatter.

Not all PIV error is small scatter. A correlation sometimes locks onto the wrong peak and
returns a displacement that has nothing to do with the soil, and one of those is worse than
it looks: between passes it becomes the offset the next pass searches around, so a single bad
vector can take its neighbourhood with it.

The test used is the **normalised median test** (Westerweel & Scarano, 2005), which is what
PIV packages use and what the field expects to see reported. A vector is judged against the
median of its eight neighbours, with the spread of those neighbours as the yardstick, so it
survives in a shear band -- where a vector genuinely differs from its surroundings -- and
fails where it differs from them for no reason.

What is left after that is scatter of a fraction of a pixel, which barely shows in the
displacement and dominates the strain, because the strain is a *difference between
neighbouring vectors* and the scatter is largest exactly there. So the field is smoothed
before it is handed on. See :func:`smooth`.
"""

from __future__ import annotations

import numpy as np


def _neighbour_stack(field: np.ndarray) -> np.ndarray:
    """The eight neighbours of every point, as ``(8, rows, cols)``, ``NaN`` past the edge."""
    out = []
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            if di == 0 and dj == 0:
                continue
            shifted = np.full(field.shape, np.nan)
            rows, cols = field.shape
            src_r = slice(max(0, -di), rows - max(0, di))
            dst_r = slice(max(0, di), rows - max(0, -di))
            src_c = slice(max(0, -dj), cols - max(0, dj))
            dst_c = slice(max(0, dj), cols - max(0, -dj))
            shifted[dst_r, dst_c] = field[src_r, src_c]
            out.append(shifted)
    return np.stack(out)


def normalised_median_residual(field: np.ndarray, floor: float = 0.1) -> np.ndarray:
    """How far each vector is from its neighbours, measured in their own spread.

    ``floor`` is the noise level of the measurement in pixels, and it is what keeps the test
    from rejecting everything in a region that is barely moving, where the spread of the
    neighbours is near zero and any difference would look enormous.
    """
    neighbours = _neighbour_stack(field)
    with np.errstate(invalid="ignore"):
        middle = np.nanmedian(neighbours, axis=0)
        spread = np.nanmedian(np.abs(neighbours - middle), axis=0)
        return np.abs(field - middle) / (spread + floor)


def find_outliers(u: np.ndarray, v: np.ndarray, threshold: float = 2.0,
                  floor: float = 0.1) -> np.ndarray:
    """Vectors that fail the normalised median test in either component."""
    with np.errstate(invalid="ignore"):
        bad = ((normalised_median_residual(u, floor) > threshold)
               | (normalised_median_residual(v, floor) > threshold))
    return np.where(np.isfinite(u) & np.isfinite(v), bad, False)


def _blur_along(field: np.ndarray, weights: np.ndarray, axis: int) -> np.ndarray:
    """One axis of a separable blur, repeating the edge value past the border."""
    radius = len(weights) // 2
    length = field.shape[axis]
    edge = np.clip(np.arange(-radius, length + radius), 0, length - 1)
    padded = np.take(field, edge, axis=axis)
    out = np.zeros_like(field)
    for offset, weight in enumerate(weights):
        out = out + weight * np.take(padded, np.arange(offset, offset + length), axis=axis)
    return out


def _shift(field: np.ndarray, di: int, dj: int) -> np.ndarray:
    """The field moved by ``di, dj``, with zeros coming in from outside."""
    out = np.zeros_like(field)
    rows, cols = field.shape
    src_r = slice(max(0, -di), rows - max(0, di))
    dst_r = slice(max(0, di), rows - max(0, -di))
    src_c = slice(max(0, -dj), cols - max(0, dj))
    dst_c = slice(max(0, dj), cols - max(0, -dj))
    out[dst_r, dst_c] = field[src_r, src_c]
    return out


def _continue_outwards(field: np.ndarray, rings: int) -> np.ndarray:
    """Extend the field into the empty space around it, one ring of points at a time.

    Each new point is the straight-line continuation of the two behind it, ``2a - b``, which
    leaves a constant gradient exactly as it was. That property is the whole purpose: see
    :func:`smooth`.
    """
    out = np.array(field, dtype=np.float64)
    for _ in range(rings):
        known = np.isfinite(out)
        filled = np.where(known, out, 0.0)
        weight = known.astype(np.float64)
        estimate = np.zeros_like(filled)
        votes = np.zeros_like(weight)
        for di, dj in ((0, 1), (0, -1), (1, 0), (-1, 0)):
            near, near_ok = _shift(filled, di, dj), _shift(weight, di, dj)
            far, far_ok = _shift(filled, 2 * di, 2 * dj), _shift(weight, 2 * di, 2 * dj)
            line = np.where(far_ok > 0, 2.0 * near - far, near)
            estimate += np.where(near_ok > 0, line, 0.0)
            votes += (near_ok > 0).astype(np.float64)
        grown = np.where((~known) & (votes > 0), estimate / np.maximum(votes, 1.0), np.nan)
        out = np.where(known, out, grown)
    return out


def _fill_the_rest(field: np.ndarray, rounds: int) -> np.ndarray:
    """Give any point still empty the average of its neighbours, so none of them is zero.

    :func:`_continue_outwards` grows a diamond, four neighbours at a time, so after ``n``
    rings the diagonal corners are still empty while the sides are full. Those corners are
    within reach of the blur at a corner of the material, and leaving them at zero pulls
    that corner towards zero -- which is the fault the continuation exists to avoid, put
    back in by the back door. So they are filled with something plausible instead.
    """
    out = field
    for _ in range(rounds):
        empty = ~np.isfinite(out)
        if not empty.any():
            break
        # the mean of the neighbours that exist, written out rather than with nanmean,
        # which warns on the points that have none
        around = _neighbour_stack(out)
        there = np.isfinite(around)
        votes = there.sum(axis=0)
        total = np.where(there, around, 0.0).sum(axis=0)
        with np.errstate(invalid="ignore", divide="ignore"):
            near = np.where(votes > 0, total / np.maximum(votes, 1), np.nan)
        out = np.where(empty & np.isfinite(near), near, out)
    return np.where(np.isfinite(out), out, 0.0)


def smooth(field: np.ndarray, sigma: float, rings: int | None = None) -> np.ndarray:
    """Blur the field with a Gaussian, without flattening the gradient at the edges.

    Why smooth at all: the displacement between two photographs of a slow test is a fraction
    of a pixel, and PIV resolves about a tenth of one. The leftover scatter is invisible in
    the displacement -- it averages away over the steps -- but the strain is the *difference
    between one vector and the next*, so there it does not average away, it dominates.
    Measured on the dam-break test, the raw field disagrees with PIVlab's by 50 % in strain
    and 8 % in displacement; smoothing brings the strain to within 4 % *and* moves each
    vector closer to PIVlab's, which is what says the scatter was noise and not soil.

    Why it is done this way matters as much. The obvious way -- averaging over whichever
    neighbours happen to exist -- is wrong at an edge, and this project has already measured
    how wrong: a kernel with no neighbours on one side leans inwards and flattens the very
    gradient being measured. So the field is first *continued outwards* past every edge, the
    blur then runs with a full set of neighbours everywhere, and the continuation is thrown
    away. A straight-line continuation leaves a constant gradient untouched, which is the
    property the strain depends on. On a field with a known gradient and 42 % of its points
    missing, as the real case has, the worst point comes back 9 % off this way against 23 %
    the other, and the median gradient is exact.

    Edges are of two kinds and both are handled: the holes inside the grid where nothing
    could be measured, and the border of the grid itself, where the material simply ends.
    Points that had no measurement stay ``NaN`` -- smoothing does not invent soil.

    ``sigma`` is in grid points, not pixels, so it means the same thing whatever the window
    and the overlap are. Zero or less returns the field untouched.
    """
    if sigma <= 0:
        return field
    radius = int(np.ceil(3 * sigma))
    offsets = np.arange(-radius, radius + 1)
    weights = np.exp(-0.5 * (offsets / sigma) ** 2)
    weights = weights / weights.sum()
    # the continuation grows a diamond, so reaching every square corner the kernel can see
    # takes twice the radius in rings
    if rings is None:
        rings = 2 * radius

    measured = np.isfinite(field)
    rows, cols = field.shape
    roomy = np.full((rows + 2 * rings, cols + 2 * rings), np.nan)
    roomy[rings:rings + rows, rings:rings + cols] = field
    grown = _fill_the_rest(_continue_outwards(roomy, rings), rings)
    for axis in (0, 1):
        grown = _blur_along(grown, weights, axis)
    blurred = grown[rings:rings + rows, rings:rings + cols]
    return np.where(measured, blurred, np.nan)


def replace(field: np.ndarray, bad: np.ndarray, rounds: int = 3) -> np.ndarray:
    """Fill the rejected vectors with the median of the neighbours that were kept.

    Done because the next pass needs a displacement to search around everywhere, not because
    an invented vector is a measurement. Which ones were replaced is reported separately so
    they can be left out of what is published.
    """
    out = np.where(bad, np.nan, field)
    for _ in range(rounds):
        missing = ~np.isfinite(out)
        if not missing.any():
            break
        with np.errstate(invalid="ignore"):
            around = np.nanmedian(_neighbour_stack(out), axis=0)
        out = np.where(missing & np.isfinite(around), around, out)
    return out
