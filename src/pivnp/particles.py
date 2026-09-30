"""Generation of the numerical particles, and tracking whether they are still in the grid."""

from __future__ import annotations

import numpy as np
from numba import njit, prange

from .config import MAX_PARTICLES_PER_SIDE, CaseConfig
from .constants import GRAVITY, LEGACY_PARTICLE_LOCAL_COORDS, PARTICLE_LOCAL_COORDS
from .mesh import Grid
from .state import Particles


def local_coordinates(npc: int, legacy_compat: bool = False) -> np.ndarray:
    """Local coordinates in [-1, 1] of the particles along each side of the cell."""
    table = LEGACY_PARTICLE_LOCAL_COORDS if legacy_compat else PARTICLE_LOCAL_COORDS
    if npc not in table:
        raise ValueError(f"NPC={npc} must be between 1 and {MAX_PARTICLES_PER_SIDE}")
    return np.array(table[npc])


def seed_positions(grid: Grid, npc: int, legacy_compat: bool = False) -> np.ndarray:
    """Positions (n_cells * npc², 2) of the initial particles.

    Order: cell row, cell within the row, particle row, particle column (the same as the
    original). The operations keep the order of the original so as to reproduce its rounding.
    """
    cols = np.arange(grid.n_cols, dtype=np.float64)
    cell_left = grid.x0 + cols * grid.dx  # XF + (J - NCF) * AXC
    cell_bottom = grid.row_y[:-1]  # YF

    g = local_coordinates(npc, legacy_compat)
    x = (cell_left + grid.dx / 2.0)[:, None] + (g * grid.dx) / 2.0
    y = (cell_bottom + grid.dy / 2.0)[:, None] + (g * grid.dy) / 2.0

    shape = (grid.n_rows, grid.n_cols, npc, npc)  # row, column, iy, ix
    px = np.broadcast_to(x[None, :, None, :], shape)
    py = np.broadcast_to(y[:, None, :, None], shape)
    return np.stack([px.ravel(), py.ravel()], axis=1)


def cell_centers(grid: Grid) -> np.ndarray:
    """Cell centres (``XP2`` in the original), in cell order."""
    x = (grid.x0 + np.arange(grid.n_cols, dtype=np.float64) * grid.dx) + grid.dx / 2.0
    y = grid.row_y[:-1] + grid.dy / 2.0
    xx, yy = np.meshgrid(x, y)
    return np.stack([xx.ravel(), yy.ravel()], axis=1)


def create_particles(config: CaseConfig, grid: Grid,
                     legacy_compat: bool = False) -> Particles:
    """Create and initialize the particles (part of ``PIVLAB_DATA``).

    On a restart (IREC = 1) the positions are loaded afterwards from the ``.REC`` file.
    """
    particles = Particles.zeros(config.n_particles, n_lost=config.n_nodes)
    if not config.restart:
        particles.position[:] = seed_positions(grid, config.particles_per_side, legacy_compat)
        particles.initial_position[:] = particles.position

    particles.potential_energy[:] = particles.mass * GRAVITY * particles.position[:, 1]
    particles.total_energy[:] = particles.potential_energy + 0.0 + 0.0
    return particles


@njit(parallel=True, cache=True)
def update_lost_flags(lost, cells, step):
    """Update ``lost`` after locating the points, with the rule of ``UCELDA``.

    * Point outside the grid: it is marked as lost.
    * Point inside and step other than 1: it stops being lost.
    * Point inside on step 1: it keeps its previous state (the rule of the original).
    """
    for i in prange(cells.size):
        if cells[i] < 0:
            lost[i] = True
        elif step != 1:
            lost[i] = False
