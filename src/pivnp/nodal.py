"""Campos nodales a partir de los datos PIVlab (antes subrutina ``VELOCIDADES``).

El paso se divide en tres funciones para poder insertar entre ellas la corrección de
contorno (``CONTOUR``):

1. :func:`load_measurements`: vuelca velocidades, humedad y marcas NaN del instante en
   los nodos PIV-NP.
2. (opcional) corrección de contorno, ver :mod:`pivnp.contour`.
3. :func:`compute_nodal_momentum_v1` / :func:`compute_nodal_momentum_v2`: cantidad de
   movimiento y su incremento en los nodos de la malla de partículas.
"""

from __future__ import annotations

import numpy as np
from numba import njit

from .constants import NODE_SIGN_X, NODE_SIGN_Y
from .mesh import Grid, cell_node, locate_points
from .particles import update_lost_flags
from .pivlab_io import Frame
from .state import Nodes, Particles


@njit(cache=True, nogil=True)
def _load_measurements_kernel(u, v, moisture_in, saturation_in, point_to_node,
                              velocity, previous_velocity, is_nan,
                              moisture_measured, saturation_measured, legacy_previous):
    for p in range(u.size):
        if legacy_previous:
            # H-01: el original copia la "velocidad anterior" con el índice del punto
            # PIVlab en vez del índice del nodo, dentro del mismo bucle que va
            # sobrescribiendo las velocidades, así que en la mitad de los nodos guarda la
            # velocidad nueva. Solo se reproduce en modo compatibilidad.
            previous_velocity[p, 0] = velocity[p, 0]
            previous_velocity[p, 1] = velocity[p, 1]
        node = point_to_node[p]
        if np.isnan(u[p]) or np.isnan(v[p]):
            velocity[node, 0] = 0.0
            velocity[node, 1] = 0.0
            is_nan[node] = 1
        else:
            velocity[node, 0] = u[p]
            velocity[node, 1] = -v[p]
            is_nan[node] = 0

        if np.isnan(moisture_in[p]) or np.isnan(saturation_in[p]):
            moisture_measured[node] = 0.0
            saturation_measured[node] = 0.0
        else:
            moisture_measured[node] = 0.0 if moisture_in[p] < 0.0 else moisture_in[p]
            saturation_measured[node] = saturation_in[p]


def load_measurements(frame: Frame, point_to_node: np.ndarray, nodes: Nodes,
                      legacy_compat: bool = False) -> None:
    """Pasa las medidas del instante (orden PIVlab) a los nodos PIV-NP.

    Guarda antes las velocidades del paso anterior, necesarias para la aceleración.
    """
    nodes.filled[:] = False
    if not legacy_compat:
        nodes.previous_velocity[:] = nodes.velocity
    _load_measurements_kernel(
        frame.u, frame.v, frame.moisture, frame.saturation, point_to_node,
        nodes.velocity, nodes.previous_velocity, nodes.is_nan,
        nodes.moisture_measured, nodes.saturation_measured, legacy_compat,
    )


@njit(cache=True, nogil=True)
def _accumulate_particle_mass(position, mass, lost, cells, cell_x, cell_y,
                              n_cols, dx, dy, nodal_mass):
    # Secuencial: la suma en el mismo orden que el original da el mismo redondeo.
    nodal_mass[:] = 0.0
    for i in range(position.shape[0]):
        if lost[i]:
            continue
        xx = 2.0 * (position[i, 0] - (cell_x[i] + dx / 2.0)) / dx
        yy = 2.0 * (position[i, 1] - (cell_y[i] + dy / 2.0)) / dy
        for j in range(4):
            fn = (1.0 + xx * NODE_SIGN_X[j]) * (1.0 + yy * NODE_SIGN_Y[j]) / 4.0
            node = cell_node(cells[i], n_cols, j)
            nodal_mass[node] = nodal_mass[node] + mass[i] * fn


def accumulate_particle_mass(particles: Particles, grid: Grid, nodes: Nodes, step: int,
                             accumulate: bool = True) -> None:
    """Relocaliza las partículas y acumula su masa en los nodos de la malla.

    Con IVERSION = 1 la masa nodal se sobrescribe con 1 más adelante, así que basta con
    ``accumulate=False``: solo se actualiza la marca ``lost``.
    """
    n = particles.position.shape[0]
    cells, cell_x, cell_y = grid_locate(particles.position, grid)
    update_lost_flags(particles.lost[:n], cells, step)
    if accumulate:
        _accumulate_particle_mass(
            particles.position, particles.mass, particles.lost, cells, cell_x, cell_y,
            grid.n_cols, grid.dx, grid.dy, nodes.mass,
        )


def grid_locate(points: np.ndarray, grid: Grid):
    """Localiza puntos en ``grid``: devuelve (celdas, x_izq, y_base)."""
    return locate_points(points, grid.x0, grid.dx, grid.n_cols, grid.row_y)


def _has_velocity(nodes: Nodes) -> np.ndarray:
    return (nodes.is_nan == 0) | nodes.filled


def compute_nodal_momentum_v1(nodes: Nodes) -> None:
    """IVERSION = 1: los nodos de cálculo son los propios puntos PIVlab (masa nodal = 1).

    Los nodos sin velocidad conservan la cantidad de movimiento del paso anterior (H-02).
    """
    n = nodes.is_nan.size
    nodes.mass[:n] = 1.0
    ok = _has_velocity(nodes)
    mass = nodes.mass[:n, None]
    nodes.momentum[:n][ok] = (nodes.velocity * mass)[ok]
    nodes.momentum_increment[:n][ok] = ((nodes.velocity - nodes.previous_velocity) * mass)[ok]
    nodes.moisture[:n][ok] = (nodes.moisture_measured * nodes.mass[:n])[ok]
    nodes.saturation[:n][ok] = (nodes.saturation_measured * nodes.mass[:n])[ok]


@njit(cache=True, nogil=True)
def _distribute_to_staggered(velocity, previous_velocity, has_velocity, lost, cells, n_cols,
                             moisture_measured, saturation_measured,
                             momentum, momentum_increment, active_count, moisture, saturation,
                             normalize):
    momentum[:] = 0.0
    momentum_increment[:] = 0.0
    active_count[:] = 0
    moisture[:] = 0.0
    saturation[:] = 0.0
    weight = 0.25
    for i in range(velocity.shape[0]):
        if lost[i] or not has_velocity[i]:
            continue
        for j in range(4):
            node = cell_node(cells[i], n_cols, j)
            for k in range(2):
                momentum[node, k] = momentum[node, k] + velocity[i, k] * weight
                momentum_increment[node, k] = (
                    momentum_increment[node, k]
                    + (velocity[i, k] - previous_velocity[i, k]) * weight
                )
            active_count[node] += 1
            moisture[node] = moisture[node] + moisture_measured[i] * weight
            saturation[node] = saturation[node] + saturation_measured[i] * weight

    if normalize:
        # H-23: un nodo del borde recibe menos de 4 aportaciones y se quedaría con una
        # fracción de la velocidad. Dividir por el peso acumulado (1 en el interior, donde
        # por tanto no cambia nada) lo convierte en la media de los puntos que sí aportan.
        for node in range(momentum.shape[0]):
            total = weight * active_count[node]
            if total > 0.0 and total != 1.0:
                momentum[node, 0] /= total
                momentum[node, 1] /= total
                momentum_increment[node, 0] /= total
                momentum_increment[node, 1] /= total
                moisture[node] /= total
                saturation[node] /= total


def compute_nodal_momentum_v2(nodes: Nodes, grid: Grid, centers: np.ndarray,
                              particles: Particles, step: int,
                              normalize: bool = False) -> None:
    """IVERSION = 2: cada punto PIVlab reparte 1/4 de su velocidad a los nodos de la celda
    desplazada que lo contiene.

    Reproduce que el original localiza los centros de celda (``XP2``) usando la marca de
    partícula ``IDONDE`` del mismo índice (H-10).
    """
    n_points = nodes.is_nan.size
    cells, _, _ = grid_locate(centers[:n_points], grid)
    update_lost_flags(particles.lost[:n_points], cells, step)
    _distribute_to_staggered(
        nodes.velocity, nodes.previous_velocity, _has_velocity(nodes), particles.lost, cells,
        grid.n_cols, nodes.moisture_measured, nodes.saturation_measured,
        nodes.momentum, nodes.momentum_increment, nodes.active_count,
        nodes.moisture, nodes.saturation, normalize,
    )
