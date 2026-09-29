import numpy as np
import pytest

from pivnp.config import parse_par
from pivnp.constants import as_fortran_real4
from pivnp.mesh import Grid, particle_grid
from pivnp.particles import (
    cell_centers,
    create_particles,
    local_coordinates,
    seed_positions,
    update_lost_flags,
)

from .legacy_reference import generate_particles
from .test_config import PAR

NC, NFIL, AXC, AYC = 12, 3, 0.212115, 0.18


@pytest.mark.parametrize("staggered", [False, True])
@pytest.mark.parametrize("npc", [1, 2, 3, 4, 5, 6])
def test_seed_positions_match_legacy_bitwise(npc, staggered):
    nch = NC // NFIL
    if staggered:
        grid = Grid(nch + 1, NFIL + 1, AXC, AYC, -AXC / 2, -AYC / 2)
    else:
        grid = Grid(nch, NFIL, AXC, AYC, 0.0, 0.0)
    expected = generate_particles(NC, NFIL, npc, AXC, AYC, local_coordinates(npc), staggered)
    np.testing.assert_array_equal(seed_positions(grid, npc), expected)


def test_local_coordinates_use_double_precision():
    assert local_coordinates(3)[0] == -2.0 / 3.0  # H-07 corregido
    assert local_coordinates(2).tolist() == [-0.5, 0.5]
    # el modo compatibilidad conserva los literales REAL*4 del original
    assert local_coordinates(3, legacy_compat=True)[0] == as_fortran_real4(-0.66666666666667)
    assert local_coordinates(4, legacy_compat=True)[0] != local_coordinates(4)[0]


def test_cell_centers_are_pivlab_points_in_staggered_grid():
    grid = Grid(3, 2, 1.0, 2.0, -0.5, -1.0)
    centers = cell_centers(grid)
    assert centers.tolist() == [[0, 0], [1, 0], [2, 0], [0, 2], [1, 2], [2, 2]]
    assert grid.locate(centers).tolist() == list(range(6))


def test_create_particles_initial_state():
    cfg = parse_par(PAR)
    p = create_particles(cfg, particle_grid(cfg))
    assert p.position.shape == (cfg.n_particles, 2)
    np.testing.assert_array_equal(p.potential_energy, 9.81 * p.position[:, 1])
    np.testing.assert_array_equal(p.mass, 1.0)
    assert not p.displacement.any() and not p.strain.any()
    assert p.lost.size >= cfg.n_nodes


def test_restart_does_not_seed_positions():
    cfg = parse_par(PAR.replace("0\t0\t0\n", "0\t1\t0\n"))
    p = create_particles(cfg, particle_grid(cfg))
    assert not p.position.any()


@pytest.mark.parametrize(("before", "found", "step", "after"), [
    (False, False, 1, True),   # fuera de la malla: perdida
    (False, False, 5, True),
    (True, True, 5, False),    # vuelve a entrar
    (True, True, 1, True),     # en el paso 1 se conserva el estado (H-09)
    (False, True, 1, False),
])
def test_lost_flag_rule(before, found, step, after):
    lost = np.array([before])
    update_lost_flags(lost, np.array([3 if found else -1]), step)
    assert lost[0] == after
