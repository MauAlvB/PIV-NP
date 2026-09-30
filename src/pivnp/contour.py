"""Velocity correction at the material boundary (the ``CONTOUR`` subroutine).

PIVlab reports no velocity (NaN) at points whose interrogation window falls partly outside
the material. That has two effects along the boundary:

* boundary particles interpolate against nodes holding zero velocity, so they move less than
  they should (and the original kept the last valid velocity there, which is what makes some
  particles "fly off");
* with IVERSION=2, nodes receiving fewer than four contributions keep only a fraction of the
  velocity they should have.

This module implements three ways of rebuilding the velocity of the points without data, so
that they can be compared. They are selected with ``ICONTOUR`` in the ``.PAR``:

===========  ==========================================================================
ICONTOUR     Method
===========  ==========================================================================
0            None (the original behaviour)
1            Average of the neighbouring nodes that do have data (8 neighbours, ``layers``)
2            Average of the particle velocities in the surrounding elements
3            Linear extrapolation from the interior towards the exterior
===========  ==========================================================================

With any method other than 0, the staggered grid (IVERSION=2) also averages the
contributions each node receives instead of leaving the sum of quarters.

Which one to use, measured on real cases (see ``docs/validacion-casos.md``): hide nodes that
*do* have a measurement and sit next to the boundary, rebuild them with each method and
compare against what PIVlab measured. Across five case/step combinations, **method 1 came
out best every time**, getting closer to the truth than leaving the zero in 76 % to 98 % of
the nodes. Method 3 follows closely on smooth fields and falls behind on abrupt ones. Method
2 is the weakest and is often worse than no correction at all: it averages particle
velocities that were themselves interpolated from the zeroed boundary nodes, so it feeds the
error back in.

Decisions shared by the three methods:

* Rebuilt nodes are flagged in ``Nodes.filled`` but do **not** change ``Nodes.is_nan``: which
  particles lie outside the material (``nan_initial``) is still decided from the measured
  data, so the correction never brings particles in mid-air to life.
* With ``ICONTOUR = 0`` the result is exactly the original one.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Protocol

import numpy as np

from .config import CaseConfig
from .mesh import Grid
from .state import Nodes, Particles

log = logging.getLogger("pivnp")

_NEIGHBOR_OFFSETS = [(dr, dc) for dr in (-1, 0, 1) for dc in (-1, 0, 1) if (dr, dc) != (0, 0)]


@dataclass(frozen=True)
class ContourContext:
    """Everything a boundary correction may need."""

    nodes: Nodes
    grid: Grid  # the grid the particles move on
    particles: Particles
    config: CaseConfig

    @property
    def measurement_shape(self) -> tuple[int, int]:
        """Rows and columns of the PIVlab point grid."""
        return self.config.n_rows + 1, self.config.n_cols + 1


class ContourCorrection(Protocol):
    def apply(self, ctx: ContourContext) -> None:
        """Modifies ``ctx.nodes.velocity`` and ``ctx.nodes.filled`` in place."""

    @property
    def normalizes_staggered(self) -> bool:
        """Whether it also averages the staggered-grid contributions."""


class NoContourCorrection:
    """No correction (the original behaviour)."""

    normalizes_staggered = False

    def apply(self, ctx: ContourContext) -> None:
        return None


@dataclass(frozen=True)
class NeighborAverageCorrection:
    """ICONTOUR=1: average of the velocities of the neighbouring nodes that have data.

    The one that measured best on real cases: it got closer to the truth than leaving the
    zero in 76 % to 98 % of the boundary nodes, depending on how smooth the field is there.
    """

    min_neighbors: int = 3
    layers: int = 1
    normalizes_staggered = True

    def __post_init__(self) -> None:
        if not 1 <= self.min_neighbors <= 8:
            raise ValueError("min_neighbors must be between 1 and 8")
        if self.layers < 1:
            raise ValueError("layers must be >= 1")

    def apply(self, ctx: ContourContext) -> None:
        shape = ctx.measurement_shape
        velocity = ctx.nodes.velocity.reshape(*shape, 2)  # a view: writes reach nodes.velocity
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
    """ICONTOUR=2: average of the particle velocities in the surrounding elements.

    A PIVlab point without data takes the mean velocity of the particles sitting in the
    elements around it (with IVERSION=2, those of the element it is the centre of). It uses
    the velocity the particles carry from the previous step, so on the first step, when those
    are still zero, it rebuilds nothing.

    Measured on real cases it is the weakest of the three, and often worse than applying no
    correction at all: those particle velocities were themselves interpolated from the zeroed
    boundary nodes, so the method feeds the boundary error back into the boundary.
    """

    min_particles: int = 1
    normalizes_staggered = True

    def __post_init__(self) -> None:
        if self.min_particles < 1:
            raise ValueError("min_particles must be >= 1")

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
        # Each particle contributes to the measurement points of its element: the 4 nodes
        # with IVERSION=1, the central point (= the cell number) with IVERSION=2.
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
    """ICONTOUR=3: linear extrapolation from the interior towards the exterior.

    Every point without data that has a neighbour with data which in turn has another one
    with data in the same direction takes ``2·v₁ − v₂`` (the slope continued), averaged over
    all available directions. Where no two aligned points exist it falls back to the
    neighbour average of ICONTOUR=1.

    On real cases it comes second: as good as the neighbour average where the field is
    smooth, clearly behind where the motion is abrupt, since continuing a slope amplifies
    whatever noise the boundary has.
    """

    layers: int = 1
    normalizes_staggered = True

    def __post_init__(self) -> None:
        if self.layers < 1:
            raise ValueError("layers must be >= 1")

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
    """Values and validity of the neighbour at (dr, dc), zero-filled outside the grid."""
    rows, cols = valid.shape
    out = np.zeros_like(values)
    ok = np.zeros_like(valid)
    src = (slice(max(dr, 0), rows + min(dr, 0)), slice(max(dc, 0), cols + min(dc, 0)))
    dst = (slice(max(-dr, 0), rows + min(-dr, 0)), slice(max(-dc, 0), cols + min(-dc, 0)))
    out[dst] = np.where(valid[src][..., None], values[src], 0.0)
    ok[dst] = valid[src]
    return out, ok


def _neighbor_sums(velocity: np.ndarray, valid: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Velocity sum and count of neighbours with data (8-connectivity) for every point."""
    total = np.zeros_like(velocity)
    count = np.zeros(valid.shape, dtype=np.int64)
    for dr, dc in _NEIGHBOR_OFFSETS:
        values, ok = _shift(velocity, valid, dr, dc)
        total += values
        count += ok
    return total, count


def make_contour_correction(icontour: int, min_neighbors: int = 3, layers: int = 1,
                            min_particles: int = 1) -> ContourCorrection:
    """Build the correction named by the ``ICONTOUR`` parameter of the ``.PAR``."""
    if icontour == 0:
        return NoContourCorrection()
    if icontour == 1:
        return NeighborAverageCorrection(min_neighbors=min_neighbors, layers=layers)
    if icontour == 2:
        log.warning(
            "ICONTOUR=2 (particle average) measured worse than applying no correction at all "
            "on the dam case, and rebuilds nothing on some steps: it averages particle "
            "velocities that were themselves interpolated from the zeroed boundary nodes. "
            "ICONTOUR=1 was the best of the three on every case tried.")
        return ParticleAverageCorrection(min_particles=min_particles)
    if icontour == 3:
        return ExtrapolationCorrection(layers=layers)
    raise ValueError(f"ICONTOUR={icontour} must be between 0 and 3")
