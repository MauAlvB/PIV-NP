"""Strain from the deformation gradient, published beside the incremental one.

The point of the finite measure is that it is exact where the incremental one is only
approximate, so the tests are the cases whose answer is known by hand: a rigid rotation
deforms nothing, simple shear has a closed-form Green-Lagrange strain, and a uniform
stretch changes the area by a known amount.
"""

from __future__ import annotations

import math
import shutil
from pathlib import Path

import numpy as np
import pytest

from pivnp.config import parse_par
from pivnp.contour import make_contour_correction
from pivnp.finite_strain import (
    area_change,
    deformation_gradient,
    equivalent_shear,
    green_lagrange,
    initial_lattice,
    rotation_degrees,
)
from pivnp.gid_writer import result_specs
from pivnp.mesh import particle_grid
from pivnp.nodal import compute_nodal_momentum_v1
from pivnp.particles import create_particles
from pivnp.simulation import Simulation, run_case
from pivnp.solver import ACTIVE, advance_particles, output_mask, update_strains
from pivnp.state import Nodes

SYNTHETIC = Path(__file__).parent / "data" / "synthetic"

# 6x4 cells of 1 m, 2x2 particles per cell, dt = 0.5
PAR = """title
block 2
24 35 2 4 1.0 1.0
block 3
0.5 {steps} 1 0 1 1 0 0 0
block 4
2000 0.4
"""


def _case(steps: int = 1):
    cfg = parse_par(PAR.format(steps=steps))
    grid = particle_grid(cfg)
    return cfg, grid, create_particles(cfg, grid), Nodes.zeros(cfg.n_nodes, grid.n_nodes)


def _node_xy(grid):
    j, i = np.divmod(np.arange(grid.n_nodes), grid.n_cols + 1)
    return i * grid.dx, j * grid.dy


def _run(cfg, grid, particles, nodes, steps: int):
    for step in range(1, steps + 1):
        compute_nodal_momentum_v1(nodes)
        advance_particles(particles, nodes, grid, cfg, step)
        update_strains(particles, nodes, grid, cfg, step)
    return output_mask(particles, grid, steps) & (particles.nan_initial == ACTIVE)


def _interior(particles, usable):
    """Away from the lattice edge, where one-sided differences are less accurate."""
    row, col, _, _ = initial_lattice(particles.initial_position)
    return usable & (row > 0) & (row < row.max()) & (col > 0) & (col < col.max())


# --- the lattice the whole thing rests on ------------------------------------------------
def test_the_initial_lattice_is_found_from_the_positions():
    cfg, grid, p, _ = _case()
    row, col, dx, dy = initial_lattice(p.initial_position)
    assert row.max() + 1 == cfg.n_rows * cfg.particles_per_side
    assert col.max() + 1 == cfg.n_cols * cfg.particles_per_side
    assert dx == pytest.approx(cfg.cell_width / cfg.particles_per_side)
    assert dy == pytest.approx(cfg.cell_height / cfg.particles_per_side)
    # every lattice place is taken exactly once
    assert len({(int(r), int(c)) for r, c in zip(row, col, strict=True)}) == row.size


# --- rigid motions, which must not deform anything ---------------------------------------
def test_a_translation_deforms_nothing():
    cfg, grid, p, nodes = _case()
    nodes.velocity[:] = [0.2, -0.1]
    usable = _run(cfg, grid, p, nodes, 1)
    inside = _interior(p, usable)

    strain = green_lagrange(deformation_gradient(p, usable))
    for component in strain:
        np.testing.assert_allclose(component[inside], 0.0, atol=1e-12)


def test_a_rigid_rotation_shows_only_what_the_time_stepping_put_there():
    """A rotation deforms nothing, and what little is left is the stepping, not the measure.

    Moving particles through a fixed rotational velocity field makes them trace a polygon
    rather than a circle, so each step pushes them outwards: after ``n`` steps of ``dtheta``
    the radius has grown by ``(1 + dtheta^2)^(n/2)``. That is forward Euler, not the strain
    measure, and the finite strain reports exactly that much and no more -- which is the
    sharpest way to show it is adding nothing of its own. On ``rotation_1P``, where the
    velocities are read from data instead of integrated here, it comes out at zero.
    """
    cfg, grid, p, nodes = _case(steps=4)
    omega, steps = 0.05, 4
    x, y = _node_xy(grid)
    nodes.velocity[:, 0] = -omega * (y - y.mean())
    nodes.velocity[:, 1] = omega * (x - x.mean())
    usable = _run(cfg, grid, p, nodes, steps)
    inside = _interior(p, usable)

    gradient = deformation_gradient(p, usable)
    xx, yy, xy = green_lagrange(gradient)

    turn = omega * cfg.dt                                   # radians per step
    stretch = (1.0 + turn * turn) ** (steps / 2.0)          # what Euler adds to the radius
    expected = (stretch * stretch - 1.0) / 2.0

    # an equal stretch in both directions and no shear: a rotation, plus the stepping
    np.testing.assert_allclose(xx[inside], expected, rtol=1e-3)
    np.testing.assert_allclose(yy[inside], expected, rtol=1e-3)
    np.testing.assert_allclose(xy[inside], 0.0, atol=1e-9)
    np.testing.assert_allclose(area_change(gradient)[inside],
                               stretch * stretch - 1.0, rtol=1e-3)

    turned = math.degrees(omega * steps * cfg.dt)
    np.testing.assert_allclose(rotation_degrees(gradient)[inside], turned, rtol=0.02)


# --- deformations with a closed-form answer ----------------------------------------------
def test_simple_shear_matches_the_closed_form():
    """For engineering shear g: E_xx = 0, E_yy = g^2/2, 2E_xy = g."""
    cfg, grid, p, nodes = _case()
    rate = 0.08
    _, y = _node_xy(grid)
    nodes.velocity[:, 0] = rate * y
    usable = _run(cfg, grid, p, nodes, 1)
    inside = _interior(p, usable)

    gamma = rate * cfg.dt
    xx, yy, xy = green_lagrange(deformation_gradient(p, usable))
    np.testing.assert_allclose(xy[inside], gamma, rtol=1e-9)
    np.testing.assert_allclose(xx[inside], 0.0, atol=1e-12)
    np.testing.assert_allclose(yy[inside], gamma * gamma / 2.0, rtol=1e-9)


def test_a_uniform_stretch_changes_the_area_by_the_right_amount():
    """u = rate*x for one step stretches by (1 + rate*dt), and the area with it."""
    cfg, grid, p, nodes = _case()
    rate = 0.06
    x, _ = _node_xy(grid)
    nodes.velocity[:, 0] = rate * x
    usable = _run(cfg, grid, p, nodes, 1)
    inside = _interior(p, usable)

    stretch = 1.0 + rate * cfg.dt
    np.testing.assert_allclose(area_change(deformation_gradient(p, usable))[inside],
                               stretch - 1.0, rtol=1e-9)
    assert np.abs(rotation_degrees(deformation_gradient(p, usable))[inside]).max() < 1e-9


# --- against the incremental measure on the real synthetic case --------------------------
def test_the_rotation_case_loses_its_artifact(workdir: Path):
    """``rotation_1P`` turns 50 degrees without deforming; the incremental measure disagrees."""
    work = workdir / "rotation"
    shutil.copytree(SYNTHETIC / "rotation_1P", work, dirs_exist_ok=True,
                    ignore=shutil.ignore_patterns("expected"))
    sim = Simulation.from_directory(work)
    sim.contour = make_contour_correction(3)
    sim.run()

    p = sim.particles
    usable = output_mask(p, sim.grid, sim.config.total_steps) & (p.nan_initial == ACTIVE)
    gradient = deformation_gradient(p, usable)
    finite = equivalent_shear(green_lagrange(gradient))

    # what the incremental measure claims, and what is actually there
    assert p.eq_strain[usable].mean() > 0.004
    assert np.abs(p.vol_strain[usable].mean()) > 0.01
    np.testing.assert_allclose(finite[usable], 0.0, atol=1e-9)
    np.testing.assert_allclose(area_change(gradient)[usable], 0.0, atol=1e-9)
    np.testing.assert_allclose(np.abs(rotation_degrees(gradient)[usable]), 50.0, atol=1e-3)


# --- what gets written -------------------------------------------------------------------
def test_the_four_blocks_are_published_and_kept_out_of_compatibility_mode():
    new = {"Finite_strain", "Fin_equi_strain", "Finite_rotation", "Finite_area"}
    assert new <= {spec.name for spec in result_specs()}
    assert not (new & {spec.name for spec in result_specs(legacy_compat=True)})


def test_they_reach_the_results_file(workdir: Path):
    work = workdir / "case"
    shutil.copytree(SYNTHETIC / "shear_1P", work, dirs_exist_ok=True,
                    ignore=shutil.ignore_patterns("expected"))
    run_case(work)
    text = (work / "zapatak.POST.RES").read_text()
    for name in ("Finite_strain", "Fin_equi_strain", "Finite_rotation", "Finite_area"):
        assert name in text, name
