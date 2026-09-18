"""Generación de partículas numéricas y seguimiento de si siguen dentro de la malla."""

from __future__ import annotations

import numpy as np
from numba import njit, prange

from .config import CaseConfig
from .constants import GRAVITY_INITIAL, MAX_GAUSS_NPC, PARTICLE_LOCAL_COORDS
from .mesh import Grid
from .state import Particles


def local_coordinates(npc: int) -> np.ndarray:
    """Coordenadas locales en [-1, 1] de las partículas de cada lado de la celda.

    Para NPC entre 7 y 10 el original usa un array sin inicializar (queda a 0 y todas las
    partículas caen en el centro de la celda); se reproduce tal cual (ver H-05).
    """
    if npc in PARTICLE_LOCAL_COORDS:
        return np.array(PARTICLE_LOCAL_COORDS[npc])
    return np.zeros(npc)


def seed_positions(grid: Grid, npc: int) -> np.ndarray:
    """Posiciones (n_cells * npc², 2) de las partículas iniciales.

    Orden: fila de celdas, celda dentro de la fila, fila de partículas, columna de partículas
    (el mismo que el original). Las operaciones mantienen el orden del original para
    reproducir su redondeo.
    """
    cols = np.arange(grid.n_cols, dtype=np.float64)
    cell_left = grid.x0 + cols * grid.dx  # XF + (J - NCF) * AXC
    cell_bottom = grid.row_y[:-1]  # YF

    if npc <= MAX_GAUSS_NPC:
        g = local_coordinates(npc)
        x = (cell_left + grid.dx / 2.0)[:, None] + (g * grid.dx) / 2.0
        y = (cell_bottom + grid.dy / 2.0)[:, None] + (g * grid.dy) / 2.0
    else:
        step_x, step_y = grid.dx / npc, grid.dy / npc
        k = np.arange(npc, dtype=np.float64)
        x = (cell_left[:, None] + k * step_x) + step_x / 2.0
        y = (cell_bottom[:, None] + k * step_y) + step_y / 2.0

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


def create_particles(config: CaseConfig, grid: Grid) -> Particles:
    """Crea e inicializa las partículas (parte de ``PIVLAB_DATA``).

    En un reinicio (IREC = 1) las posiciones de la malla se cargan después desde el
    archivo ``.REC``; aquí solo se colocan las partículas de seguimiento PTV.
    """
    n0 = config.n_base_particles
    particles = Particles.zeros(config.n_particles, n_lost=config.n_nodes)
    if not config.restart:
        particles.position[:n0] = seed_positions(grid, config.particles_per_side)
    if config.tracking:
        particles.position[n0:] = config.tracking.positions()

    particles.potential_energy[:] = particles.mass * GRAVITY_INITIAL * particles.position[:, 1]
    particles.total_energy[:] = particles.potential_energy + 0.0 + 0.0
    return particles


@njit(parallel=True, cache=True)
def update_lost_flags(lost, cells, step):
    """Actualiza ``lost`` tras localizar los puntos, con la regla de ``UCELDA``.

    * Punto fuera de la malla: se marca como perdido.
    * Punto dentro y paso distinto de 1: deja de estar perdido.
    * Punto dentro en el paso 1: conserva su estado anterior (ver H-09).
    """
    for i in prange(cells.size):
        if cells[i] < 0:
            lost[i] = True
        elif step != 1:
            lost[i] = False
