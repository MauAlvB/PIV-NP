import numpy as np

from pivnp.mesh import Grid
from pivnp.nodal import (
    compute_nodal_momentum_v1,
    compute_nodal_momentum_v2,
    load_measurements,
)
from pivnp.particles import cell_centers
from pivnp.pivlab_io import Frame, pivlab_to_node
from pivnp.state import Nodes, Particles

NAN = np.nan


def _frame(u, v, moisture=None, saturation=None):
    n = len(u)
    zeros = np.zeros(n)
    return Frame(1, None, np.array(u, float), np.array(v, float),
                 zeros if moisture is None else np.array(moisture, float),
                 zeros if saturation is None else np.array(saturation, float))


def _legacy_velocidades_mapping(u, v, conn, vel, vel_prev):
    """Bucle de VELOCIDADES tal cual (índices base 0)."""
    vel, vel_prev = vel.copy(), vel_prev.copy()
    for i in range(len(u)):
        vel_prev[i] = vel[i]
        if np.isnan(u[i]) or np.isnan(v[i]):
            vel[conn[i]] = 0.0
        else:
            vel[conn[i]] = (u[i], -v[i])
    return vel, vel_prev


def test_load_measurements_reproduces_legacy_previous_velocity_quirk():
    conn = pivlab_to_node(2, 2)  # 9 nodos
    rng = np.random.default_rng(3)
    nodes = Nodes.zeros(9, 9)
    nodes.velocity[:] = rng.normal(size=(9, 2))
    before = nodes.velocity.copy()
    u, v = rng.normal(size=9), rng.normal(size=9)
    u[4] = NAN
    expected_vel, expected_prev = _legacy_velocidades_mapping(u, v, conn, before,
                                                              nodes.previous_velocity)
    load_measurements(_frame(u, v), conn, nodes, legacy_compat=True)
    np.testing.assert_array_equal(nodes.velocity, expected_vel)
    np.testing.assert_array_equal(nodes.previous_velocity, expected_prev)
    # H-01: la "velocidad anterior" de algún nodo ya es la del paso actual
    assert (nodes.previous_velocity != before).any()
    assert nodes.is_nan[conn[4]] == 1 and nodes.is_nan.sum() == 1


def test_load_measurements_moisture_clamp_and_nan():
    conn = np.arange(3)
    nodes = Nodes.zeros(3, 3)
    frame = _frame([0, 0, 0], [0, 0, 0], moisture=[-0.2, 0.3, NAN], saturation=[0.1, 0.2, 0.3])
    load_measurements(frame, conn, nodes)
    assert nodes.moisture_measured.tolist() == [0.0, 0.3, 0.0]
    assert nodes.saturation_measured.tolist() == [0.1, 0.2, 0.0]


def test_momentum_v1_keeps_previous_value_for_nan_nodes():
    nodes = Nodes.zeros(2, 2)
    nodes.velocity[:] = [[1.0, 2.0], [3.0, 4.0]]
    compute_nodal_momentum_v1(nodes)
    nodes.velocity[:] = [[5.0, 6.0], [0.0, 0.0]]
    nodes.previous_velocity[:] = [[1.0, 2.0], [3.0, 4.0]]
    nodes.is_nan[1] = 1
    compute_nodal_momentum_v1(nodes)
    assert nodes.momentum.tolist() == [[5.0, 6.0], [3.0, 4.0]]  # H-02: nodo 1 conserva
    assert nodes.momentum_increment.tolist() == [[4.0, 4.0], [3.0, 4.0]]  # también conserva
    assert nodes.mass.tolist() == [1.0, 1.0]


def test_momentum_v1_uses_filled_nodes():
    nodes = Nodes.zeros(1, 1)
    nodes.velocity[:] = [[2.0, 0.0]]
    nodes.is_nan[0] = 1
    nodes.filled[0] = True
    compute_nodal_momentum_v1(nodes)
    assert nodes.momentum.tolist() == [[2.0, 0.0]]


def test_momentum_v2_averages_four_surrounding_points():
    # 2x2 celdas PIVlab -> 3x3 puntos -> malla desplazada de 3x3 celdas y 4x4 nodos.
    grid = Grid(3, 3, 1.0, 1.0, -0.5, -0.5)
    nodes = Nodes.zeros(9, grid.n_nodes)
    nodes.velocity[:, 0] = np.arange(9.0)
    particles = Particles.zeros(1, n_lost=9)
    compute_nodal_momentum_v2(nodes, grid, cell_centers(grid), particles, step=2)
    # El nodo 5 de la malla desplazada (fila 1, col 1) toca los puntos 0, 1, 3 y 4.
    assert nodes.momentum[5, 0] == 0.25 * (0 + 1 + 3 + 4)
    assert nodes.active_count[5] == 4 and nodes.active_count[0] == 1
    # Los puntos sin datos no aportan
    nodes.is_nan[4] = 1
    compute_nodal_momentum_v2(nodes, grid, cell_centers(grid), particles, step=2)
    assert nodes.momentum[5, 0] == 0.25 * (0 + 1 + 3)
