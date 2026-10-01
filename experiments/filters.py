"""Candidate filters, and the source that applies one to another source's velocities.

A filter is a function of the velocity field laid out on the PIVlab grid, with ``NaN``
where nothing was measured. Every one of them has to leave the gaps alone: averaging a
boundary point with the void outside the material would drag it towards zero, which is the
very mistake the contour correction exists to undo.

Nothing here is part of PIV-NP. It plugs in through ``DisplacementSource``, so the solver
never knows a filter was applied.
"""
from __future__ import annotations

from collections.abc import Callable, Iterator
from dataclasses import replace

import numpy as np

from pivnp.pivlab_io import pivlab_to_node
from pivnp.sources import Frame

#: A filter takes the field as a (rows, cols) array with NaN in the gaps, and returns it
#: with the same shape and the same gaps.
Filter = Callable[[np.ndarray], np.ndarray]


def no_filter(field: np.ndarray) -> np.ndarray:
    """The baseline everything is measured against."""
    return field


def _weighted_pass(field: np.ndarray, kernel: np.ndarray) -> np.ndarray:
    """Convolve ignoring the gaps: every point is the weighted mean of the neighbours it has.

    Where a point has no measured neighbour at all it is left as it was, and the gaps stay
    gaps, so the edge of the material does not get pulled towards zero.
    """
    measured = np.isfinite(field)
    values = np.where(measured, field, 0.0)
    weights = measured.astype(np.float64)
    half = kernel.shape[0] // 2

    total = np.zeros_like(values)
    total_weight = np.zeros_like(weights)
    for i in range(-half, half + 1):
        for j in range(-half, half + 1):
            w = kernel[i + half, j + half]
            if w == 0.0:
                continue
            total += w * _shift(values, i, j)
            total_weight += w * _shift(weights, i, j)

    out = np.where(total_weight > 0, total / np.maximum(total_weight, 1e-300), field)
    return np.where(measured, out, np.nan)


def _shift(field: np.ndarray, rows: int, cols: int) -> np.ndarray:
    """Shift without wrapping around: what leaves the grid is zero, not the far edge."""
    out = np.zeros_like(field)
    src_rows = slice(max(0, -rows), field.shape[0] - max(0, rows))
    dst_rows = slice(max(0, rows), field.shape[0] - max(0, -rows))
    src_cols = slice(max(0, -cols), field.shape[1] - max(0, cols))
    dst_cols = slice(max(0, cols), field.shape[1] - max(0, -cols))
    out[dst_rows, dst_cols] = field[src_rows, src_cols]
    return out


def moving_average(passes: int = 1) -> Filter:
    """The floor of the exercise: a plain 3x3 mean, applied ``passes`` times."""
    kernel = np.ones((3, 3))

    def apply(field: np.ndarray) -> np.ndarray:
        out = field
        for _ in range(passes):
            out = _weighted_pass(out, kernel)
        return out

    return apply


def gaussian(sigma: float) -> Filter:
    """Gaussian of width ``sigma`` in grid cells, truncated at three sigma."""
    half = max(1, int(np.ceil(3 * sigma)))
    offsets = np.arange(-half, half + 1)
    line = np.exp(-0.5 * (offsets / sigma) ** 2)
    kernel = np.outer(line, line)

    def apply(field: np.ndarray) -> np.ndarray:
        return _weighted_pass(field, kernel)

    return apply


def median(half_width: int = 1) -> Filter:
    """Median of the neighbourhood: the one filter that removes a spurious vector instead
    of spreading it.

    PIV noise is not only small scatter; it is also the occasional vector that is simply
    wrong, from a bad correlation peak. An average drags that error into its neighbours,
    a median throws it away.
    """
    def apply(field: np.ndarray) -> np.ndarray:
        rows, cols = field.shape
        stack = []
        for i in range(-half_width, half_width + 1):
            for j in range(-half_width, half_width + 1):
                shifted = np.full(field.shape, np.nan)
                src_r = slice(max(0, -i), rows - max(0, i))
                dst_r = slice(max(0, i), rows - max(0, -i))
                src_c = slice(max(0, -j), cols - max(0, j))
                dst_c = slice(max(0, j), cols - max(0, -j))
                shifted[dst_r, dst_c] = field[src_r, src_c]
                stack.append(shifted)
        neighbourhood = np.stack(stack)
        with np.errstate(invalid="ignore"):
            out = np.nanmedian(neighbourhood, axis=0)
        return np.where(np.isfinite(field), out, np.nan)

    return apply


def _product(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """``a.T @ b`` without matrix multiplication.

    This machine's numpy crashes (0xc06d007f) on any matmul, even 3x3: the environment has
    three BLAS libraries and four OpenMP runtimes fighting each other. The matrices here are
    6x6, so element-wise sums cost nothing and sidestep the problem entirely.
    """
    out = np.empty((a.shape[1], b.shape[1]))
    for i in range(a.shape[1]):
        for j in range(b.shape[1]):
            out[i, j] = np.sum(a[:, i] * b[:, j])
    return out


def _solve(matrix: np.ndarray, rhs: np.ndarray) -> np.ndarray:
    """Gaussian elimination with partial pivoting, for a small dense system."""
    a = np.array(matrix, dtype=np.float64)
    b = np.array(rhs, dtype=np.float64)
    n = a.shape[0]
    for column in range(n):
        pivot = column + int(np.argmax(np.abs(a[column:, column])))
        if pivot != column:
            a[[column, pivot]] = a[[pivot, column]]
            b[[column, pivot]] = b[[pivot, column]]
        if abs(a[column, column]) < 1e-300:
            raise ValueError("singular system")
        for row in range(column + 1, n):
            factor = a[row, column] / a[column, column]
            a[row, column:] -= factor * a[column, column:]
            b[row] -= factor * b[column]
    out = np.zeros_like(b)
    for row in range(n - 1, -1, -1):
        behind = float(np.sum(a[row, row + 1:] * out[row + 1:]))
        out[row] = (b[row] - behind) / a[row, row]
    return out


def savitzky_golay(half_width: int = 2) -> Filter:
    """Fit a quadratic surface to the neighbourhood and keep its value at the centre.

    Unlike an average, a polynomial fit does not flatten a genuine gradient, which is why
    it is the usual choice when the field is going to be differentiated afterwards. The
    weights do not depend on the data, so they are worked out once and applied as a
    convolution; where the window runs into a gap or the edge there are not enough points
    to fit a quadratic, and a plain weighted average is used instead.
    """
    offsets = np.arange(-half_width, half_width + 1)
    dy, dx = np.meshgrid(offsets, offsets, indexing="ij")
    dy, dx = dy.ravel().astype(float), dx.ravel().astype(float)
    design = np.stack([np.ones_like(dx), dx, dy, dx * dx, dx * dy, dy * dy], axis=1)
    # the value at the centre is the first coefficient: row 0 of (A'A)^-1 A'
    normal = _product(design, design)
    weights = np.array([_solve(normal, design[k, :])[0] for k in range(design.shape[0])])
    kernel = weights.reshape(dy.size // offsets.size, offsets.size)
    fallback = gaussian(half_width / 2.0)

    def apply(field: np.ndarray) -> np.ndarray:
        measured = np.isfinite(field)
        values = np.where(measured, field, 0.0)
        half = half_width
        total = np.zeros_like(values)
        complete = np.ones(field.shape, dtype=bool)
        for i in range(-half, half + 1):
            for j in range(-half, half + 1):
                total += kernel[i + half, j + half] * _shift(values, i, j)
                complete &= _shift(measured.astype(float), i, j) > 0
        smooth = fallback(field)
        return np.where(measured, np.where(complete, total, smooth), np.nan)

    return apply


class FilteredSource:
    """Wraps a displacement source and smooths its velocities before handing them over."""

    def __init__(self, inner, n_cols: int, n_rows: int, velocity_filter: Filter) -> None:
        self.inner = inner
        self.filter = velocity_filter
        self.shape = (n_rows + 1, n_cols + 1)
        self.to_node = pivlab_to_node(n_cols, n_rows)
        self.from_node = np.argsort(self.to_node)

    # the attributes the analysis sets on a source
    @property
    def moisture(self) -> bool:
        return self.inner.moisture

    @moisture.setter
    def moisture(self, value: bool) -> None:
        self.inner.moisture = value

    @property
    def images(self):
        return self.inner.images

    @images.setter
    def images(self, value) -> None:
        self.inner.images = value

    def mesh_in_metres(self, step: int = 1):
        return self.inner.mesh_in_metres(step)

    def frame_interval(self) -> float | None:
        return self.inner.frame_interval()

    def _smooth(self, values: np.ndarray) -> np.ndarray:
        """From the PIVlab point order to the grid, filter, and back."""
        grid = np.full(self.shape, np.nan)
        grid.ravel()[self.to_node] = values
        filtered = self.filter(grid)
        return filtered.ravel()[self.to_node]

    def frames(self, steps: range) -> Iterator[Frame]:
        for frame in self.inner.frames(steps):
            yield replace(frame, u=self._smooth(frame.u), v=self._smooth(frame.v))


def _extrapolate_outwards(field: np.ndarray, rings: int) -> np.ndarray:
    """Continue the field into the void, one ring of cells at a time.

    A symmetric kernel at the edge of the material has no neighbours on one side. Averaging
    only over the side it does have leans the result inwards and flattens the gradient,
    which is exactly what was measured. Giving it a plausible continuation to chew on
    instead removes the reason for the bias: each new point is the linear continuation of
    the two behind it, which leaves a constant gradient untouched.
    """
    out = np.array(field, dtype=np.float64)
    for _ in range(rings):
        known = np.isfinite(out)
        filled = np.where(known, out, 0.0)
        weight = known.astype(np.float64)
        estimate = np.zeros_like(filled)
        votes = np.zeros_like(weight)
        for i, j in ((0, 1), (0, -1), (1, 0), (-1, 0)):
            near = _shift(filled, i, j)
            near_ok = _shift(weight, i, j)
            far = _shift(filled, 2 * i, 2 * j)
            far_ok = _shift(weight, 2 * i, 2 * j)
            linear = np.where(far_ok > 0, 2.0 * near - far, near)   # 2a - b, or just a
            estimate += np.where(near_ok > 0, linear, 0.0)
            votes += (near_ok > 0).astype(np.float64)
        new = np.where((~known) & (votes > 0), estimate / np.maximum(votes, 1.0), np.nan)
        out = np.where(known, out, new)
    return out


def across_the_edge(inner: Filter, rings: int = 3) -> Filter:
    """Apply a filter after continuing the field past the edge of the material.

    The edge is of two kinds and both have to be handled: the gaps inside the grid, where
    PIVlab measured nothing, and the border of the grid itself, which in a case like
    ``shear-block`` is where the material ends. So the field is first laid into a larger
    array with room to grow, the continuation is built, the filter runs where it now has
    neighbours on every side, and the result is cropped back.
    """
    def apply(field: np.ndarray) -> np.ndarray:
        measured = np.isfinite(field)
        rows, cols = field.shape
        roomy = np.full((rows + 2 * rings, cols + 2 * rings), np.nan)
        roomy[rings:rings + rows, rings:rings + cols] = field
        widened = _extrapolate_outwards(roomy, rings)
        smoothed = inner(widened)[rings:rings + rows, rings:rings + cols]
        return np.where(measured, smoothed, np.nan)

    return apply


#: The candidates of phase 1, by the name they are reported under.
CANDIDATES: dict[str, Filter] = {
    "none": no_filter,
    "median 3x3": median(1),
    "average 3x3": moving_average(1),
    "average 3x3 x2": moving_average(2),
    "gaussian s=0.7": gaussian(0.7),
    "gaussian s=1.0": gaussian(1.0),
    "gaussian s=1.5": gaussian(1.5),
    "savitzky-golay 5x5": savitzky_golay(2),
    # the same filters, but told what lies past the edge of the material
    "edge: median 3x3": across_the_edge(median(1)),
    "edge: median 5x5": across_the_edge(median(2), rings=4),
    "edge: median then avg": across_the_edge(
        lambda f: moving_average(1)(median(1)(f)), rings=4),
    "edge: average 3x3": across_the_edge(moving_average(1)),
    "edge: average x2": across_the_edge(moving_average(2)),
    "edge: gaussian s=1.0": across_the_edge(gaussian(1.0)),
    "edge: gaussian s=1.5": across_the_edge(gaussian(1.5), rings=5),
    "edge: sav-golay 5x5": across_the_edge(savitzky_golay(2)),
}
