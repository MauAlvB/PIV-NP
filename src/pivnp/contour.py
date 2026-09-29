"""Corrección de velocidades en el contorno del material (subrutina ``CONTOUR``).

PIVlab no da velocidad (NaN) en los puntos cuya ventana de interrogación cae parcialmente
fuera del material. Eso provoca dos efectos en el borde:

* las partículas del borde interpolan con nodos a velocidad 0, así que se mueven menos de
  lo que deberían (H-02: además, el original conserva ahí la velocidad del último paso con
  dato, lo que hace que algunas partículas "salgan volando");
* con IVERSION=2, los nodos que reciben menos de cuatro aportaciones se quedan con una
  fracción de la velocidad que les corresponde (H-23).

Este módulo implementa tres formas de reconstruir la velocidad de los puntos sin dato, para
poder compararlas. Se eligen con ``ICONTOUR`` en el ``.PAR``:

===========  ==========================================================================
ICONTOUR     Método
===========  ==========================================================================
0            Ninguno (comportamiento del original)
1            Media de los nodos vecinos con dato (8 vecinos, ``layers`` capas)
2            Media de las velocidades de las partículas de los elementos de alrededor
3            Extrapolación lineal desde el interior hacia el exterior
===========  ==========================================================================

Además, con cualquier método distinto de 0 la malla desplazada (IVERSION=2) normaliza las
aportaciones de cada nodo, que es la corrección de H-23.

Decisiones comunes a los tres métodos:

* Los nodos reconstruidos se marcan en ``Nodes.filled`` pero **no** cambian ``Nodes.is_nan``:
  la decisión de qué partículas están fuera del material (``nan_initial``) sigue usando los
  datos medidos, así que la corrección no activa partículas en el aire.
* Con ``ICONTOUR = 0`` el resultado es exactamente el del original.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

import numpy as np

from .config import CaseConfig
from .mesh import Grid
from .state import Nodes, Particles

_NEIGHBOR_OFFSETS = [(dr, dc) for dr in (-1, 0, 1) for dc in (-1, 0, 1) if (dr, dc) != (0, 0)]


@dataclass(frozen=True)
class ContourContext:
    """Todo lo que puede necesitar una corrección de contorno."""

    nodes: Nodes
    grid: Grid  # malla sobre la que se mueven las partículas
    particles: Particles
    config: CaseConfig

    @property
    def measurement_shape(self) -> tuple[int, int]:
        """Filas y columnas de la rejilla de puntos PIVlab."""
        return self.config.n_rows + 1, self.config.n_cols + 1


class ContourCorrection(Protocol):
    def apply(self, ctx: ContourContext) -> None:
        """Modifica ``ctx.nodes.velocity`` y ``ctx.nodes.filled`` in situ."""

    @property
    def normalizes_staggered(self) -> bool:
        """Si corrige también el reparto en la malla desplazada (H-23)."""


class NoContourCorrection:
    """Sin corrección (comportamiento del original)."""

    normalizes_staggered = False

    def apply(self, ctx: ContourContext) -> None:
        return None


@dataclass(frozen=True)
class NeighborAverageCorrection:
    """ICONTOUR=1: media de las velocidades de los nodos vecinos con dato."""

    min_neighbors: int = 3
    layers: int = 1
    normalizes_staggered = True

    def __post_init__(self) -> None:
        if not 1 <= self.min_neighbors <= 8:
            raise ValueError("min_neighbors debe estar entre 1 y 8")
        if self.layers < 1:
            raise ValueError("layers debe ser >= 1")

    def apply(self, ctx: ContourContext) -> None:
        shape = ctx.measurement_shape
        velocity = ctx.nodes.velocity.reshape(*shape, 2)  # vista: escribe en nodes.velocity
        valid = (ctx.nodes.is_nan == 0).reshape(shape).copy()
        filled = np.zeros(shape, dtype=bool)

        for _ in range(self.layers):
            total, count = _neighbor_sums(velocity, valid)
            candidates = ~valid & (count >= self.min_neighbors)
            if not candidates.any():
                break
            velocity[candidates] = total[candidates] / count[candidates, None]
            valid |= candidates
            filled |= candidates

        ctx.nodes.filled[:] = filled.ravel()


@dataclass(frozen=True)
class ParticleAverageCorrection:
    """ICONTOUR=2: media de las velocidades de las partículas de los elementos de alrededor.

    A un punto PIVlab sin dato se le asigna la velocidad media de las partículas que hay en
    los elementos que lo rodean (con IVERSION=2, las del elemento del que es centro). Usa la
    velocidad que las partículas traen del paso anterior, así que en el primer paso, cuando
    todavía valen cero, no reconstruye nada.
    """

    min_particles: int = 1
    normalizes_staggered = True

    def __post_init__(self) -> None:
        if self.min_particles < 1:
            raise ValueError("min_particles debe ser >= 1")

    def apply(self, ctx: ContourContext) -> None:
        nodes, grid, p = ctx.nodes, ctx.grid, ctx.particles
        missing = np.flatnonzero(nodes.is_nan != 0)
        if missing.size == 0:
            nodes.filled[:] = False
            return

        alive = ~p.lost[: p.position.shape[0]]
        cells = grid.locate(p.position)
        usable = alive & (cells >= 0)
        total = np.zeros((nodes.velocity.shape[0], 2))
        count = np.zeros(nodes.velocity.shape[0], dtype=np.int64)
        # Cada partícula aporta a los puntos de medida de su elemento: los 4 nodos con
        # IVERSION=1, el punto central (= número de celda) con IVERSION=2.
        if ctx.config.mesh_version == 1:
            targets = grid.cell_nodes(cells[usable])  # (n, 4)
            for k in range(4):
                np.add.at(total, targets[:, k], p.velocity[usable])
                np.add.at(count, targets[:, k], 1)
        else:
            np.add.at(total, cells[usable], p.velocity[usable])
            np.add.at(count, cells[usable], 1)

        enough = missing[count[missing] >= self.min_particles]
        nodes.velocity[enough] = total[enough] / count[enough, None]
        nodes.filled[:] = False
        nodes.filled[enough] = True


@dataclass(frozen=True)
class ExtrapolationCorrection:
    """ICONTOUR=3: extrapolación lineal desde el interior hacia el exterior.

    A cada punto sin dato con un vecino con dato que a su vez tiene otro vecino con dato en
    la misma dirección se le asigna ``2·v₁ − v₂`` (continuación de la pendiente), promediando
    todas las direcciones disponibles. Donde no hay dos puntos alineados se usa la media de
    los vecinos, como en ICONTOUR=1.
    """

    layers: int = 1
    normalizes_staggered = True

    def __post_init__(self) -> None:
        if self.layers < 1:
            raise ValueError("layers debe ser >= 1")

    def apply(self, ctx: ContourContext) -> None:
        shape = ctx.measurement_shape
        velocity = ctx.nodes.velocity.reshape(*shape, 2)
        valid = (ctx.nodes.is_nan == 0).reshape(shape).copy()
        filled = np.zeros(shape, dtype=bool)

        for _ in range(self.layers):
            total = np.zeros((*shape, 2))
            count = np.zeros(shape, dtype=np.int64)
            for dr, dc in _NEIGHBOR_OFFSETS:
                first, ok1 = _shift(velocity, valid, dr, dc)
                second, ok2 = _shift(velocity, valid, 2 * dr, 2 * dc)
                usable = ok1 & ok2
                total += np.where(usable[..., None], 2.0 * first - second, 0.0)
                count += usable
            fallback_total, fallback_count = _neighbor_sums(velocity, valid)
            candidates = ~valid & ((count > 0) | (fallback_count > 0))
            if not candidates.any():
                break
            extrapolated = candidates & (count > 0)
            averaged = candidates & (count == 0)
            velocity[extrapolated] = total[extrapolated] / count[extrapolated, None]
            velocity[averaged] = (fallback_total[averaged]
                                  / fallback_count[averaged, None])
            valid |= candidates
            filled |= candidates

        ctx.nodes.filled[:] = filled.ravel()


def _shift(values: np.ndarray, valid: np.ndarray, dr: int, dc: int):
    """Valores y validez del vecino en la dirección (dr, dc), con ceros fuera de la malla."""
    rows, cols = valid.shape
    out = np.zeros_like(values)
    ok = np.zeros_like(valid)
    src = (slice(max(dr, 0), rows + min(dr, 0)), slice(max(dc, 0), cols + min(dc, 0)))
    dst = (slice(max(-dr, 0), rows + min(-dr, 0)), slice(max(-dc, 0), cols + min(-dc, 0)))
    out[dst] = np.where(valid[src][..., None], values[src], 0.0)
    ok[dst] = valid[src]
    return out, ok


def _neighbor_sums(velocity: np.ndarray, valid: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Suma de velocidades y número de vecinos con dato (8-conectividad) de cada punto."""
    total = np.zeros_like(velocity)
    count = np.zeros(valid.shape, dtype=np.int64)
    for dr, dc in _NEIGHBOR_OFFSETS:
        values, ok = _shift(velocity, valid, dr, dc)
        total += values
        count += ok
    return total, count


def make_contour_correction(icontour: int, min_neighbors: int = 3, layers: int = 1,
                            min_particles: int = 1) -> ContourCorrection:
    """Corrección según el parámetro ICONTOUR del ``.PAR``."""
    if icontour == 0:
        return NoContourCorrection()
    if icontour == 1:
        return NeighborAverageCorrection(min_neighbors=min_neighbors, layers=layers)
    if icontour == 2:
        return ParticleAverageCorrection(min_particles=min_particles)
    if icontour == 3:
        return ExtrapolationCorrection(layers=layers)
    raise ValueError(f"ICONTOUR={icontour} debe estar entre 0 y 3")
