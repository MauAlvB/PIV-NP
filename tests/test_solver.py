import math

import numpy as np
import pytest

from pivnp.config import parse_par
from pivnp.constants import J2_THRESHOLD
from pivnp.mesh import particle_grid
from pivnp.nodal import compute_nodal_momentum_v1
from pivnp.particles import create_particles
from pivnp.solver import (
    ACTIVE,
    NAN_AT_START,
    advance_particles,
    deviatoric_q,
    output_mask,
    update_strains,
)
from pivnp.state import Nodes

from .legacy_reference import invar2

# Malla de 3x2 celdas de 1 m, 2x2 partículas por celda, dt = 0.5
PAR = """titulo
bloque 2
6 12 2 2 1.0 1.0
bloque 3
0.5 10 1 0 1 1 0 0 0
bloque 4
2000 0.4
"""


@pytest.fixture
def case():
    cfg = parse_par(PAR)
    grid = particle_grid(cfg)
    particles = create_particles(cfg, grid)
    nodes = Nodes.zeros(cfg.n_nodes, grid.n_nodes)
    return cfg, grid, particles, nodes


def _node_xy(grid):
    j, i = np.divmod(np.arange(grid.n_nodes), grid.n_cols + 1)
    return i * grid.dx, j * grid.dy


def _step(cfg, grid, particles, nodes, step):
    compute_nodal_momentum_v1(nodes)
    advance_particles(particles, nodes, grid, cfg, step)
    update_strains(particles, nodes, grid, cfg, step)


@pytest.mark.parametrize("values", [
    (1e-3, -2e-3, 0.0, 5e-4), (0.0, 0.0, 0.0, 0.0), (1e-6, 1e-6, 0.0, 0.0), (3.0, 1.0, 2.0, -1.0),
])
def test_deviatoric_q_matches_invar2(values):
    assert deviatoric_q(*values) == invar2(*values, J2_THRESHOLD)


def test_deviatoric_q_uniaxial():
    # Deformación uniaxial e: q = sqrt(3 J2) = e
    assert deviatoric_q(0.3, 0.0, 0.0, 0.0) == pytest.approx(0.3)


def test_uniform_velocity_translates_without_strain(case):
    cfg, grid, p, nodes = case
    start = p.position.copy()
    nodes.velocity[:] = [0.2, -0.1]
    _step(cfg, grid, p, nodes, 1)
    np.testing.assert_allclose(p.displacement, np.tile([0.1, -0.05], (p.position.shape[0], 1)))
    np.testing.assert_allclose(p.position, start + [0.1, -0.05])
    np.testing.assert_allclose(p.velocity, np.tile([0.2, -0.1], (p.position.shape[0], 1)))
    assert not p.strain.any() and not p.eq_strain.any()


def test_linear_velocity_field_gives_uniform_strain(case):
    cfg, grid, p, nodes = case
    x, _ = _node_xy(grid)
    rate = 0.01  # u = rate * x  ->  eps_xx = rate * dt por paso
    nodes.velocity[:, 0] = rate * x
    for step in (1, 2):
        _step(cfg, grid, p, nodes, step)
    located = output_mask(p, grid, 2)
    np.testing.assert_allclose(p.strain[located, 0], 2 * rate * cfg.dt)
    np.testing.assert_allclose(p.vol_strain[located], 2 * rate * cfg.dt)
    q = math.sqrt(3 * ((2 / 3) ** 2 + 2 * (1 / 3) ** 2) / 2) * 2 * rate * cfg.dt
    np.testing.assert_allclose(p.eq_strain[located], 2 * q / 3)


def test_particles_in_all_nan_cells_are_frozen(case):
    cfg, grid, p, nodes = case
    nodes.velocity[:] = [1.0, 0.0]
    nodes.is_nan[[0, 1, 4, 5]] = 1  # los 4 nodos de la celda 0
    cell0 = grid.locate(p.position) == 0
    start = p.position.copy()
    _step(cfg, grid, p, nodes, 1)
    assert (p.nan_initial[cell0] == NAN_AT_START).all()
    assert (p.nan_initial[~cell0] == ACTIVE).all()
    np.testing.assert_array_equal(p.position[cell0], start[cell0])
    assert not p.displacement[cell0].any()


def test_particle_leaving_the_mesh_is_lost(case):
    cfg, grid, p, nodes = case
    nodes.velocity[:] = [10.0, 0.0]  # 5 m en un paso: todas salen por la derecha
    _step(cfg, grid, p, nodes, 1)
    assert not output_mask(p, grid, 2).any()


def test_strain_is_uniform_inside_each_cell(case):
    cfg, grid, p, nodes = case
    rng = np.random.default_rng(5)
    nodes.velocity[:] = rng.normal(size=(cfg.n_nodes, 2)) * 0.01
    _step(cfg, grid, p, nodes, 1)
    cells = grid.locate(p.position - p.position_increment)
    for cell in np.unique(cells[cells >= 0]):
        same = cells == cell
        assert np.ptp(p.strain_increment[same], axis=0).max() == 0.0
