"""Actualización de las partículas en cada paso (antes subrutinas ``SOLMOV`` e ``INVAR2``).

Cada partícula se calcula de forma independiente, así que los bucles se ejecutan en
paralelo con ``numba.prange``. Las expresiones mantienen el orden de operaciones del
original para reproducir su redondeo.
"""

from __future__ import annotations

import math

import numpy as np
from numba import njit, prange

from .config import CaseConfig
from .constants import (
    GRAVITY,
    J2_THRESHOLD,
    LEGACY_GRAVITY_STEP,
    LEGACY_J2_THRESHOLD,
    MACHINE_EPSILON,
    NODE_SIGN_X,
    NODE_SIGN_Y,
)
from .mesh import Grid, cell_node
from .nodal import grid_locate
from .particles import update_lost_flags
from .state import Nodes, Particles

#: Valores de ``Particles.nan_initial``.
ACTIVE, NAN_AT_START = 0, 1


@njit(cache=True, nogil=True, inline="always")
def deviatoric_q(sx, sy, sz, sxy, threshold=J2_THRESHOLD):
    """Tensión (o deformación) desviadora q = sqrt(3·J2); 0 si J2 <= ``threshold`` (``INVAR2``)."""
    mean = (sx + sy + sz) / 3.0
    dx, dy, dz = sx - mean, sy - mean, sz - mean
    j2 = (dx * dx + dy * dy + dz * dz) / 2.0 + sxy * sxy
    return math.sqrt(3.0 * j2) if j2 > threshold else 0.0


@njit(parallel=True, cache=True)
def _advance_kernel(lost, cells, cell_x, cell_y, n_cols, dx, dy, dt,
                    momentum, momentum_increment, nodal_mass, node_is_nan,
                    nodal_moisture, nodal_saturation,
                    position, velocity, acceleration, increment, displacement, step_disp,
                    moisture, saturation, nan_initial, nan_step,
                    step, restart, with_moisture, mesh_version, legacy_restart_displacement):
    for i in prange(position.shape[0]):
        if not lost[i]:
            cell = cells[i]
            xx = 2 * (position[i, 0] - (cell_x[i] + dx / 2.0)) / dx
            yy = 2 * (position[i, 1] - (cell_y[i] + dy / 2.0)) / dy
            for k in range(2):
                increment[i, k] = 0.0
                velocity[i, k] = 0.0
                acceleration[i, k] = 0.0
                step_disp[i, k] = 0.0
            nan_step[i] = 0
            if step == 1 and not restart:
                nan_initial[i] = 0
            if step == 1 and legacy_restart_displacement:
                # El original sumaba el desplazamiento acumulado al "instantáneo" en el
                # primer paso tras un reinicio.
                step_disp[i, 0] = displacement[i, 0]
                step_disp[i, 1] = displacement[i, 1]
            nan_corners = 0
            moisture[i] = 0.0
            saturation[i] = 0.0

            for j in range(4):
                fn = (1.0 + xx * NODE_SIGN_X[j]) * (1.0 + yy * NODE_SIGN_Y[j]) / 4.0
                node = cell_node(cell, n_cols, j)
                if nodal_mass[node] >= MACHINE_EPSILON:
                    f1 = fn
                    f2 = fn * dt
                else:
                    f1 = 0.0
                    f2 = 0.0
                for k in range(2):
                    velocity[i, k] = velocity[i, k] + momentum[node, k] * f1
                    acceleration[i, k] = acceleration[i, k] + momentum_increment[node, k] * f1 / dt
                    increment[i, k] = increment[i, k] + momentum[node, k] * f2
                    step_disp[i, k] = step_disp[i, k] + momentum[node, k] * f2
                    displacement[i, k] = displacement[i, k] + momentum[node, k] * f2
                if with_moisture:
                    moisture[i] = moisture[i] + nodal_moisture[node] * f1
                    saturation[i] = saturation[i] + nodal_saturation[node] * f1
                if mesh_version == 1:
                    nan_corners += node_is_nan[node]
                elif node_is_nan[cell] == 1:  # celda desplazada centrada en un nodo sin datos
                    if step == 1 and not restart:
                        nan_initial[i] = NAN_AT_START
                    nan_step[i] = 1

            if mesh_version == 1 and nan_corners == 4:  # los 4 nodos sin datos
                if step == 1 and not restart:
                    nan_initial[i] = NAN_AT_START
                else:
                    nan_step[i] = 1

        if nan_initial[i] == NAN_AT_START:  # sin datos desde el inicio: no se mueve
            increment[i, 0] = 0.0
            increment[i, 1] = 0.0
            displacement[i, 0] = 0.0
            displacement[i, 1] = 0.0
            lost[i] = True


def advance_particles(particles: Particles, nodes: Nodes, grid: Grid, config: CaseConfig,
                      step: int, legacy_compat: bool = False) -> None:
    """Interpola velocidad, aceleración y desplazamiento de las partículas desde los nodos.

    Primera mitad de ``SOLMOV``. Marca como ``NAN_AT_START`` las partículas que en el primer
    paso están en celdas sin ningún dato (fuera del material).
    """
    p = particles
    n = p.position.shape[0]
    cells, cell_x, cell_y = grid_locate(p.position, grid)
    update_lost_flags(p.lost[:n], cells, step)
    _advance_kernel(
        p.lost, cells, cell_x, cell_y, grid.n_cols, grid.dx, grid.dy, config.dt,
        nodes.momentum, nodes.momentum_increment, nodes.mass, nodes.is_nan,
        nodes.moisture, nodes.saturation,
        p.position, p.velocity, p.acceleration, p.position_increment, p.displacement,
        p.step_displacement, p.moisture, p.saturation, p.nan_initial, p.nan_step,
        step, config.restart, config.moisture, config.mesh_version,
        config.restart and legacy_compat,
    )


@njit(parallel=True, cache=True)
def _strain_kernel(lost, cells, n_cols, dx, dy, dt, momentum, nodal_mass,
                   position, increment, velocity, mass, strain, strain_inc,
                   vol_strain, vol_strain_inc, potential, kinetic, total, moisture,
                   eq_strain, eq_strain_inc, legacy_divide_by_mass, gravity, j2_threshold):
    for i in prange(position.shape[0]):
        if lost[i]:
            continue
        cell = cells[i]
        d0 = 0.0
        d1 = 0.0
        d2 = 0.0
        for j in range(4):
            node = cell_node(cell, n_cols, j)
            # Derivadas de las funciones de forma en el centro del elemento: la deformación
            # es uniforme en cada celda.
            dndx = NODE_SIGN_X[j] * 0.5 / dx
            dndy = NODE_SIGN_Y[j] * 0.5 / dy
            if nodal_mass[node] >= MACHINE_EPSILON:
                # El original dividía por la masa nodal, que con IVERSION=2 vale
                # aproximadamente NPC² y deja las deformaciones a escala de la velocidad
                # dividida por NPC². Solo se reproduce en modo compatibilidad.
                if legacy_divide_by_mass:
                    fx = dndx * dt / nodal_mass[node]
                    fy = dndy * dt / nodal_mass[node]
                else:
                    fx = dndx * dt
                    fy = dndy * dt
            else:
                fx = 0.0
                fy = 0.0
            d0 = d0 + momentum[node, 0] * fx
            d1 = d1 + momentum[node, 1] * fy
            d2 = d2 + momentum[node, 0] * fy + momentum[node, 1] * fx
        if abs(d0) < MACHINE_EPSILON and d0 != 0.0:
            d0 = 0.0
        if abs(d1) < MACHINE_EPSILON and d1 != 0.0:
            d1 = 0.0
        if abs(d2) < MACHINE_EPSILON and d2 != 0.0:
            d2 = 0.0
        strain_inc[i, 0] = d0
        strain_inc[i, 1] = d1
        strain_inc[i, 2] = d2
        strain[i, 0] = strain[i, 0] + d0
        strain[i, 1] = strain[i, 1] + d1
        strain[i, 2] = strain[i, 2] + d2

        vol_strain[i] = strain[i, 0] + strain[i, 1]
        vol_strain_inc[i] = d0 + d1
        potential[i] = mass[i] * gravity * position[i, 1]
        kinetic[i, 0] = 0.5 * mass[i] * velocity[i, 0] * velocity[i, 0]
        kinetic[i, 1] = 0.5 * mass[i] * velocity[i, 1] * velocity[i, 1]
        total[i] = potential[i] + kinetic[i, 0] + kinetic[i, 1]
        if moisture[i] < 0.0:
            moisture[i] = 0.0

        position[i, 0] = position[i, 0] + increment[i, 0]
        position[i, 1] = position[i, 1] + increment[i, 1]

        # Deformación de corte equivalente: 2/3·q (con la deformación de corte ingenieril / 2)
        eq_strain[i] = 2.0 * deviatoric_q(strain[i, 0], strain[i, 1], strain[i, 3],
                                          strain[i, 2] / 2.0, j2_threshold) / 3.0
        eq_strain_inc[i] = 2.0 * deviatoric_q(d0, d1, 0.0, d2 / 2.0, j2_threshold) / 3.0


def update_strains(particles: Particles, nodes: Nodes, grid: Grid, config: CaseConfig,
                   step: int, legacy_compat: bool = False) -> None:
    """Deformaciones, energías y nueva posición de cada partícula (segunda mitad de SOLMOV)."""
    p = particles
    n = p.position.shape[0]
    cells, _, _ = grid_locate(p.position, grid)
    update_lost_flags(p.lost[:n], cells, step)
    _strain_kernel(
        p.lost, cells, grid.n_cols, grid.dx, grid.dy, config.dt, nodes.momentum, nodes.mass,
        p.position, p.position_increment, p.velocity, p.mass, p.strain, p.strain_increment,
        p.vol_strain, p.vol_strain_increment, p.potential_energy, p.kinetic_energy,
        p.total_energy, p.moisture, p.eq_strain, p.eq_strain_increment, legacy_compat,
        LEGACY_GRAVITY_STEP if legacy_compat else GRAVITY,
        LEGACY_J2_THRESHOLD if legacy_compat else J2_THRESHOLD,
    )


def output_mask(particles: Particles, grid: Grid, step: int) -> np.ndarray:
    """Partículas localizadas en la malla al imprimir (``IDONDE /= -1`` en IMPRES_GiD)."""
    n = particles.position.shape[0]
    cells, _, _ = grid_locate(particles.position, grid)
    update_lost_flags(particles.lost[:n], cells, step)
    return ~particles.lost[:n]


@njit(parallel=True, cache=True)
def _count_nan_kernel(lost, cells, n_cols, node_is_nan, mesh_version, out):
    for i in prange(lost.size):
        if lost[i]:
            out[i] = 0
            continue
        if mesh_version == 1:  # nodos del elemento sin dato (0 a 4)
            total = 0
            for j in range(4):
                total += node_is_nan[cell_node(cells[i], n_cols, j)]
            out[i] = total
        else:  # el punto PIVlab del centro del elemento no tiene dato (0 o 1)
            out[i] = node_is_nan[cells[i]]


def count_nan_nodes(particles: Particles, nodes: Nodes, grid: Grid, mesh_version: int,
                    step: int) -> np.ndarray:
    """Cuántos datos faltan alrededor de cada partícula (resultado ``NaNs``)."""
    n = particles.position.shape[0]
    cells, _, _ = grid_locate(particles.position, grid)
    update_lost_flags(particles.lost[:n], cells, step)
    _count_nan_kernel(particles.lost[:n], cells, grid.n_cols, nodes.is_nan, mesh_version,
                      particles.missing_data)
    return particles.missing_data
