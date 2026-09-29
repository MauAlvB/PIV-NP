"""Generación de partículas numéricas y seguimiento de si siguen dentro de la malla."""

from __future__ import annotations

import numpy as np
from numba import njit, prange

from .config import MAX_PARTICLES_PER_SIDE, CaseConfig
from .constants import GRAVITY, LEGACY_PARTICLE_LOCAL_COORDS, PARTICLE_LOCAL_COORDS
from .mesh import Grid
from .state import Particles


def local_coordinates(npc: int, legacy_compat: bool = False) -> np.ndarray:
    """Coordenadas locales en [-1, 1] de las partículas de cada lado de la celda."""
    table = LEGACY_PARTICLE_LOCAL_COORDS if legacy_compat else PARTICLE_LOCAL_COORDS
    if npc not in table:
        raise ValueError(f"NPC={npc} debe estar entre 1 y {MAX_PARTICLES_PER_SIDE}")
    return np.array(table[npc])


def seed_positions(grid: Grid, npc: int, legacy_compat: bool = False) -> np.ndarray:
    """Posiciones (n_cells * npc², 2) de las partículas iniciales.

    Orden: fila de celdas, celda dentro de la fila, fila de partículas, columna de partículas
    (el mismo que el original). Las operaciones mantienen el orden del original para
    reproducir su redondeo.
    """
    cols = np.arange(grid.n_cols, dtype=np.float64)
    cell_left = grid.x0 + cols * grid.dx  # XF + (J - NCF) * AXC
    cell_bottom = grid.row_y[:-1]  # YF

    g = local_coordinates(npc, legacy_compat)
    x = (cell_left + grid.dx / 2.0)[:, None] + (g * grid.dx) / 2.0
    y = (cell_bottom + grid.dy / 2.0)[:, None] + (g * grid.dy) / 2.0

    shape = (grid.n_rows, grid.n_cols, npc, npc)  # fila, columna, iy, ix
    px = np.broadcast_to(x[None, :, None, :], shape)
    py = np.broadcast_to(y[:, None, :, None], shape)
    return np.stack([px.ravel(), py.ravel()], axis=1)


def cell_centers(grid: Grid) -> np.ndarray:
    """Centros de las celdas (``XP2`` en el original), en orden de celda."""
    x = (grid.x0 + np.arange(grid.n_cols, dtype=np.float64) * grid.dx) + grid.dx / 2.0
    y = grid.row_y[:-1] + grid.dy / 2.0
    xx, yy = np.meshgrid(x, y)
    return np.stack([xx.ravel(), yy.ravel()], axis=1)


def create_particles(config: CaseConfig, grid: Grid,
                     legacy_compat: bool = False) -> Particles:
    """Crea e inicializa las partículas (parte de ``PIVLAB_DATA``).

    En un reinicio (IREC = 1) las posiciones se cargan después desde el archivo ``.REC``.
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
    """Actualiza ``lost`` tras localizar los puntos, con la regla de ``UCELDA``.

    * Punto fuera de la malla: se marca como perdido.
    * Punto dentro y paso distinto de 1: deja de estar perdido.
    * Punto dentro en el paso 1: conserva su estado anterior (regla del original).
    """
    for i in prange(cells.size):
        if cells[i] < 0:
            lost[i] = True
        elif step != 1:
            lost[i] = False
