"""Corrección de velocidades en el contorno del material (subrutina ``CONTOUR``).

En el original ``CONTOUR`` se llama pero no existe todavía. El problema que pretende
resolver: PIVlab no da velocidad (NaN) en los puntos cuya ventana de interrogación cae
fuera del material. Esos nodos entran en la interpolación con velocidad 0 (o con la del
paso anterior), así que las partículas de las celdas del borde se mueven menos de lo que
deberían.

Propuesta (:class:`NeighborAverageCorrection`): antes de calcular la cantidad de movimiento
nodal, a cada nodo sin datos que tenga al menos ``min_neighbors`` vecinos válidos (de sus 8
vecinos) se le asigna la media de sus velocidades. Se puede repetir ``layers`` veces para
avanzar varias filas de nodos hacia el exterior.

Decisiones de diseño:

* Se promedian **nodos** y no partículas: la velocidad de las partículas se interpola desde
  los nodos, así que promediar partículas sería circular y más caro.
* Los nodos rellenados se marcan en ``Nodes.filled`` pero **no** cambian ``Nodes.is_nan``:
  la decisión de qué partículas están fuera del material (``nan_initial``) sigue usando los
  datos medidos, así que la corrección no activa partículas en el aire.
* Con ``ICONTOUR = 0`` se usa :class:`NoContourCorrection`, que no hace nada y reproduce
  exactamente el resultado del original.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

import numpy as np

from .state import Nodes

_NEIGHBOR_OFFSETS = [(dr, dc) for dr in (-1, 0, 1) for dc in (-1, 0, 1) if (dr, dc) != (0, 0)]


class ContourCorrection(Protocol):
    def apply(self, nodes: Nodes, n_cols: int, n_rows: int) -> None:
        """Modifica ``nodes.velocity`` y ``nodes.filled`` in situ."""


class NoContourCorrection:
    """Sin corrección (comportamiento del original)."""

    def apply(self, nodes: Nodes, n_cols: int, n_rows: int) -> None:
        return None


@dataclass(frozen=True)
class NeighborAverageCorrection:
    """Rellena nodos sin datos del borde con la media de sus vecinos con datos."""

    min_neighbors: int = 3
    layers: int = 1

    def __post_init__(self) -> None:
        if not 1 <= self.min_neighbors <= 8:
            raise ValueError("min_neighbors debe estar entre 1 y 8")
        if self.layers < 1:
            raise ValueError("layers debe ser >= 1")

    def apply(self, nodes: Nodes, n_cols: int, n_rows: int) -> None:
        shape = (n_rows + 1, n_cols + 1)
        velocity = nodes.velocity.reshape(*shape, 2)  # vista: escribe en nodes.velocity
        valid = (nodes.is_nan == 0).reshape(shape).copy()
        filled = np.zeros(shape, dtype=bool)

        for _ in range(self.layers):
            total, count = _neighbor_sums(velocity, valid)
            candidates = ~valid & (count >= self.min_neighbors)
            if not candidates.any():
                break
            velocity[candidates] = total[candidates] / count[candidates, None]
            valid |= candidates
            filled |= candidates

        nodes.filled[:] = filled.ravel()


def _neighbor_sums(velocity: np.ndarray, valid: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Suma de velocidades y número de vecinos válidos (8-conectividad) de cada nodo."""
    rows, cols = valid.shape
    padded_v = np.zeros((rows + 2, cols + 2, 2))
    padded_ok = np.zeros((rows + 2, cols + 2), dtype=bool)
    padded_v[1:-1, 1:-1] = np.where(valid[..., None], velocity, 0.0)
    padded_ok[1:-1, 1:-1] = valid

    total = np.zeros((rows, cols, 2))
    count = np.zeros((rows, cols), dtype=np.int64)
    for dr, dc in _NEIGHBOR_OFFSETS:
        window = (slice(1 + dr, rows + 1 + dr), slice(1 + dc, cols + 1 + dc))
        total += padded_v[window]
        count += padded_ok[window]
    return total, count


def make_contour_correction(icontour: int, min_neighbors: int = 3,
                            layers: int = 1) -> ContourCorrection:
    """Corrección según el parámetro ICONTOUR del ``.PAR``."""
    if icontour == 0:
        return NoContourCorrection()
    return NeighborAverageCorrection(min_neighbors=min_neighbors, layers=layers)
