"""Mallas rectangulares estructuradas y localización de puntos (antes ``UCELDA``).

Convenciones (índices base 0):

* Las celdas se numeran por filas, de abajo a arriba y de izquierda a derecha:
  ``cell = row * n_cols + col``.
* Los nodos también: ``node = row * (n_cols + 1) + col``.
* Los 4 nodos de una celda van en el orden del elemento bilineal
  (abajo-izq, abajo-der, arriba-izq, arriba-der).

Hay dos mallas posibles (``IVERSION``):

1. *PIV-NP*: coincide con la malla PIVlab; sus nodos son los puntos PIVlab.
2. *Desplazada*: tiene una fila y una columna más y está desplazada media celda, de modo
   que cada punto PIVlab queda en el centro de una celda.
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
        """Ordenada de la base de cada fila (``YF`` en el original), ``n_rows + 1`` valores."""
        return self.y0 + np.arange(self.n_rows + 1, dtype=np.float64) * self.dy

    def cell_nodes(self, cells: np.ndarray) -> np.ndarray:
        """Nodos (n, 4) de cada celda."""
        cells = np.asarray(cells)
        row, col = np.divmod(cells, self.n_cols)
        first = row * (self.n_cols + 1) + col
        above = first + self.n_cols + 1
        return np.stack([first, first + 1, above, above + 1], axis=-1)

    def locate(self, points: np.ndarray) -> np.ndarray:
        """Celda que contiene cada punto (n, 2), o -1 si queda fuera de la malla."""
        cells, _, _ = locate_points(np.ascontiguousarray(points, dtype=np.float64),
                                    self.x0, self.dx, self.n_cols, self.row_y)
        return cells


def pivnp_grid(config: CaseConfig) -> Grid:
    """Malla PIV-NP, coincidente con la de PIVlab (IVERSION = 1)."""
    return Grid(config.n_cols, config.n_rows, config.cell_width, config.cell_height, 0.0, 0.0)


def staggered_grid(config: CaseConfig) -> Grid:
    """Malla desplazada media celda, con los puntos PIVlab en los centros (IVERSION = 2)."""
    return Grid(
        config.n_cols + 1,
        config.n_rows + 1,
        config.cell_width,
        config.cell_height,
        -config.cell_width / 2,
        -config.cell_height / 2,
    )


def particle_grid(config: CaseConfig) -> Grid:
    """Malla sobre la que se mueven las partículas según IVERSION."""
    return pivnp_grid(config) if config.mesh_version == 1 else staggered_grid(config)


@njit(cache=True, nogil=True)
def locate_point(x, y, x0, dx, n_cols, row_y):
    """Devuelve ``(celda, x_izq, y_base)`` o ``(-1, 0, 0)`` si el punto está fuera.

    Reproduce la búsqueda lineal del original, incluidas sus comparaciones exactas en
    coma flotante, pero en O(1): se estima la fila/columna y se comprueban los vecinos.
    La fila ``r`` contiene ``y`` si ``row_y[r] <= y < row_y[r+1]``; la columna ``c`` si
    ``x_c <= x < x_c + dx`` con ``x_c = x0 + dx*c`` (se elige la primera que cumpla).
    """
    n_rows = row_y.size - 1
    if not (y >= row_y[0] and y < row_y[n_rows]):  # también descarta NaN
        return -1, 0.0, 0.0

    row = int((y - row_y[0]) / (row_y[1] - row_y[0]))
    row = min(max(row, 0), n_rows - 1)
    while y < row_y[row]:
        row -= 1
    while y >= row_y[row + 1]:
        row += 1

    guess = (x - x0) / dx
    if not (guess > -2.0 and guess < n_cols + 2.0):  # fuera o NaN
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
    """Versión vectorizada y paralela de :func:`locate_point`."""
    n = points.shape[0]
    cells = np.empty(n, dtype=np.int64)
    xc = np.empty(n, dtype=np.float64)
    yc = np.empty(n, dtype=np.float64)
    for i in prange(n):
        cells[i], xc[i], yc[i] = locate_point(points[i, 0], points[i, 1], x0, dx, n_cols, row_y)
    return cells, xc, yc


@njit(cache=True, nogil=True, inline="always")
def cell_node(cell, n_cols, local):
    """Nodo global ``local`` (0..3) de ``cell``."""
    row = cell // n_cols
    first = row * (n_cols + 1) + cell - row * n_cols
    if local == 0:
        return first
    if local == 1:
        return first + 1
    if local == 2:
        return first + n_cols + 1
    return first + n_cols + 2
