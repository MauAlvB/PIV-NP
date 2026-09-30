"""Structured rectangular grids and point location (formerly ``UCELDA``).

Conventions (0-based indices):

* Cells are numbered by rows, bottom to top and left to right:
  ``cell = row * n_cols + col``.
* Nodes too: ``node = row * (n_cols + 1) + col``.
* The 4 nodes of a cell follow the order of the bilinear element
  (bottom-left, bottom-right, top-left, top-right).

There are two possible grids (``IVERSION``):

1. *PIV-NP*: the same as the PIVlab grid; its nodes are the PIVlab points.
2. *Staggered*: it has one more row and one more column and is shifted half a cell, so that
   every PIVlab point sits at the centre of a cell.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property

import numpy as np
from numba import njit, prange

from .config import CaseConfig


@dataclass(frozen=True)
class Grid:
    n_cols: int
    n_rows: int
    dx: float
    dy: float
    x0: float
    y0: float

    @property
    def n_cells(self) -> int:
        return self.n_cols * self.n_rows

    @property
    def n_nodes(self) -> int:
        return (self.n_cols + 1) * (self.n_rows + 1)

    @cached_property
    def row_y(self) -> np.ndarray:
        """Ordinate of the base of each row (``YF`` in the original), ``n_rows + 1`` values."""
        return self.y0 + np.arange(self.n_rows + 1, dtype=np.float64) * self.dy

    def cell_nodes(self, cells: np.ndarray) -> np.ndarray:
        """Nodes (n, 4) of each cell."""
        cells = np.asarray(cells)
        row, col = np.divmod(cells, self.n_cols)
        first = row * (self.n_cols + 1) + col
        above = first + self.n_cols + 1
        return np.stack([first, first + 1, above, above + 1], axis=-1)

    def locate(self, points: np.ndarray) -> np.ndarray:
        """Cell containing each point (n, 2), or -1 if it falls outside the grid."""
        cells, _, _ = locate_points(np.ascontiguousarray(points, dtype=np.float64),
                                    self.x0, self.dx, self.n_cols, self.row_y)
        return cells


def pivnp_grid(config: CaseConfig) -> Grid:
    """PIV-NP grid, the same as the PIVlab one (IVERSION = 1)."""
    return Grid(config.n_cols, config.n_rows, config.cell_width, config.cell_height, 0.0, 0.0)


def staggered_grid(config: CaseConfig) -> Grid:
    """Grid shifted half a cell, with the PIVlab points at the centres (IVERSION = 2)."""
    return Grid(
        config.n_cols + 1,
        config.n_rows + 1,
        config.cell_width,
        config.cell_height,
        -config.cell_width / 2,
        -config.cell_height / 2,
    )


def particle_grid(config: CaseConfig) -> Grid:
    """Grid the particles move on, according to IVERSION."""
    return pivnp_grid(config) if config.mesh_version == 1 else staggered_grid(config)


@njit(cache=True, nogil=True)
def locate_point(x, y, x0, dx, n_cols, row_y):
    """Return ``(cell, x_left, y_base)`` or ``(-1, 0, 0)`` if the point is outside.

    It reproduces the linear search of the original, including its exact floating-point
    comparisons, but in O(1): the row/column is estimated and the neighbours are checked.
    Row ``r`` contains ``y`` if ``row_y[r] <= y < row_y[r+1]``; column ``c`` if
    ``x_c <= x < x_c + dx`` with ``x_c = x0 + dx*c`` (the first match wins).
    """
    n_rows = row_y.size - 1
    if not (y >= row_y[0] and y < row_y[n_rows]):  # also rejects NaN
        return -1, 0.0, 0.0

    row = int((y - row_y[0]) / (row_y[1] - row_y[0]))
    row = min(max(row, 0), n_rows - 1)
    while y < row_y[row]:
        row -= 1
    while y >= row_y[row + 1]:
        row += 1

    guess = (x - x0) / dx
    if not (guess > -2.0 and guess < n_cols + 2.0):  # outside or NaN
        return -1, 0.0, 0.0
    first = max(int(np.floor(guess)) - 1, 0)
    last = min(int(np.floor(guess)) + 1, n_cols - 1)
    for col in range(first, last + 1):
        xc = x0 + dx * col
        if x >= xc and x < xc + dx:
            return row * n_cols + col, xc, row_y[row]
    return -1, 0.0, 0.0


@njit(parallel=True, cache=True)
def locate_points(points, x0, dx, n_cols, row_y):
    """Vectorized, parallel version of :func:`locate_point`."""
    n = points.shape[0]
    cells = np.empty(n, dtype=np.int64)
    xc = np.empty(n, dtype=np.float64)
    yc = np.empty(n, dtype=np.float64)
    for i in prange(n):
        cells[i], xc[i], yc[i] = locate_point(points[i, 0], points[i, 1], x0, dx, n_cols, row_y)
    return cells, xc, yc


@njit(cache=True, nogil=True, inline="always")
def cell_node(cell, n_cols, local):
    """Global node ``local`` (0..3) of ``cell``."""
    row = cell // n_cols
    first = row * (n_cols + 1) + cell - row * n_cols
    if local == 0:
        return first
    if local == 1:
        return first + 1
    if local == 2:
        return first + n_cols + 1
    return first + n_cols + 2
