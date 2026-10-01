"""Vorticity and accumulated rotation.

The curl of the velocity field is what the strains cannot tell you: a rigid rotation does
not deform the material, yet accumulating linear strain increments reports strain for it
(see ``tests/test_synthetic.py::test_rigid_body_rotation``). These two results separate the
two things, so the tests here check them against fields whose answer is known by hand.
"""

from __future__ import annotations

import math
import shutil
from pathlib import Path

import numpy as np
import pytest

from pivnp.config import parse_par
from pivnp.contour import make_contour_correction
from pivnp.gid_writer import result_specs
from pivnp.mesh import particle_grid
from pivnp.nodal import compute_nodal_momentum_v1
from pivnp.particles import create_particles
from pivnp.restart import read_restart
from pivnp.simulation import RunOptions, Simulation, run_case
from pivnp.solver import (
    advance_particles,
    output_mask,
    rotation_angle,
    update_strains,
    vorticity_number,
)
from pivnp.state import Nodes

SYNTHETIC = Path(__file__).parent / "data" / "synthetic"

# 3x2 cells of 1 m, 2x2 particles per cell, dt = 0.5
PAR = """title
block 2
6 12 2 2 1.0 1.0
block 3
0.5 {steps} 1 0 1 1 0 0 0
block 4
2000 0.4
"""


def _case(steps: int = 1):
    cfg = parse_par(PAR.format(steps=steps))
    grid = particle_grid(cfg)
    return cfg, grid, create_particles(cfg, grid), Nodes.zeros(cfg.n_nodes, grid.n_nodes)


def _step(cfg, grid, particles, nodes, step):
    compute_nodal_momentum_v1(nodes)
    advance_particles(particles, nodes, grid, cfg, step)
    update_strains(particles, nodes, grid, cfg, step)


def _node_xy(grid):
    j, i = np.divmod(np.arange(grid.n_nodes), grid.n_cols + 1)
    return i * grid.dx, j * grid.dy


# --- fields whose curl is known by hand --------------------------------------------------
def test_a_uniform_field_does_not_rotate():
    cfg, grid, p, nodes = _case()
    nodes.velocity[:] = [0.2, -0.1]
    _step(cfg, grid, p, nodes, 1)
    assert not p.vorticity.any()
    assert not p.rotation.any()


def test_pure_stretching_does_not_rotate():
    """u = rate*x deforms without turning: the curl stays zero."""
    cfg, grid, p, nodes = _case()
    x, _ = _node_xy(grid)
    nodes.velocity[:, 0] = 0.01 * x
    _step(cfg, grid, p, nodes, 1)
    assert p.strain[:, 0].any()          # it does deform
    assert not p.vorticity.any()         # but it does not turn


def test_simple_shear_turns_at_half_its_shear_rate():
    """Simple shear is pure shear plus rotation, so its curl is not zero."""
    cfg, grid, p, nodes = _case()
    rate = 0.04
    _, y = _node_xy(grid)
    nodes.velocity[:, 0] = rate * y      # u grows with height
    _step(cfg, grid, p, nodes, 1)

    located = output_mask(p, grid, 1)
    # curl = dv/dx - du/dy = -rate
    np.testing.assert_allclose(p.vorticity[located], -rate, atol=1e-12)
    expected = math.degrees(-rate / 2 * cfg.dt)
    np.testing.assert_allclose(p.rotation[located], expected, atol=1e-12)
    np.testing.assert_allclose(p.strain[located, 2], rate * cfg.dt, atol=1e-12)


def test_rigid_rotation_turns_without_deforming():
    """u = -W*y, v = W*x: curl = 2W, and no strain at all."""
    cfg, grid, p, nodes = _case(steps=2)
    omega = 0.03                                   # rad/s
    x, y = _node_xy(grid)
    nodes.velocity[:, 0] = -omega * y
    nodes.velocity[:, 1] = omega * x
    for step in (1, 2):
        _step(cfg, grid, p, nodes, step)

    located = output_mask(p, grid, 2)
    np.testing.assert_allclose(p.vorticity[located], 2 * omega, atol=1e-12)
    turned = math.degrees(omega * 2 * cfg.dt)      # W * total time
    np.testing.assert_allclose(p.rotation[located], turned, atol=1e-12)
    assert np.abs(p.strain[located, :3]).max() < 1e-15


def test_the_rotation_accumulates_over_the_steps():
    cfg, grid, p, nodes = _case(steps=3)
    x, y = _node_xy(grid)
    nodes.velocity[:, 0] = -0.02 * y
    nodes.velocity[:, 1] = 0.02 * x
    seen = []
    for step in (1, 2, 3):
        _step(cfg, grid, p, nodes, step)
        seen.append(p.rotation[0])
    np.testing.assert_allclose(seen, [seen[0], 2 * seen[0], 3 * seen[0]], rtol=1e-12)


# --- the real case, where the answer is known --------------------------------------------
def test_a_finite_rotation_is_tracked_to_a_thousandth_of_a_degree(workdir: Path):
    """``rotation_1P`` turns 1 deg/s for 50 s, so the answer is 50 degrees.

    This is the case where the strains lie: the same analysis reports an equivalent shear
    strain of about 0.005 where the material is not deforming at all. The rotation, with the
    boundary rebuilt, lands within 0.01 deg of the truth and is the same on every particle.
    """
    work = workdir / "rotation"
    shutil.copytree(SYNTHETIC / "rotation_1P", work, dirs_exist_ok=True,
                    ignore=shutil.ignore_patterns("expected"))
    sim = Simulation.from_directory(work)
    sim.contour = make_contour_correction(3)       # rebuild the boundary by extrapolation
    sim.run()

    p = sim.particles
    alive = output_mask(p, sim.grid, sim.config.total_steps) & (p.nan_initial == 0)
    turned = np.abs(p.rotation[alive])             # the case turns clockwise
    expected = float(sim.config.total_steps - 1)   # the first instant does not rotate

    assert turned.mean() == pytest.approx(expected, abs=0.01)
    assert turned.max() - turned.min() < 1e-6      # the same on every particle
    # and the strain, on the very same particles, claims a deformation that is not there
    assert p.eq_strain[alive].mean() > 0.004


def test_the_rotation_survives_a_restart(workdir: Path):
    """Without carrying it in the .REC, a restarted analysis would turn from zero again.

    Four steps and then four more, with the second batch renumbered from 1 as a second
    stage of a test would be, has to give what eight steps in a row give.
    """
    def build(directory: Path, steps: int, restart: int, first: int = 1) -> Path:
        source = SYNTHETIC / "rotation_1P"
        (directory / "pivlab").mkdir(parents=True, exist_ok=True)
        shutil.copy(source / "PIV-NP.TXT", directory)
        for k in range(steps):
            shutil.copy(source / "pivlab" / f"datos ({first + k}).txt",
                        directory / "pivlab" / f"datos ({k + 1}).txt")
        par = (source / "zapatak.PAR").read_text().splitlines()
        par[4] = f"         1. {steps} 1 0 1 1 0 {restart} 0"
        (directory / "zapatak.PAR").write_text("\n".join(par) + "\n")
        return directory

    whole = build(workdir / "whole", steps=8, restart=0)
    run_case(whole)

    part = build(workdir / "part", steps=4, restart=0)
    run_case(part)
    build(part, steps=4, restart=1, first=5)      # instants 5..8, renumbered to 1..4
    run_case(part)

    restarted = read_restart(part / "zapatak.REC")
    assert restarted.rotation is not None
    assert np.abs(restarted.rotation).max() > 1.0  # it really did turn
    np.testing.assert_allclose(restarted.rotation,
                               read_restart(whole / "zapatak.REC").rotation, rtol=1e-12)


# --- how the deformation divides between shear and rotation ------------------------------
def test_pure_shear_is_zero_and_zero_degrees():
    """Stretching along one axis and shortening along the other: it deforms, it does not turn."""
    cfg, grid, p, nodes = _case()
    x, y = _node_xy(grid)
    nodes.velocity[:, 0] = 0.01 * x
    nodes.velocity[:, 1] = -0.01 * y
    _step(cfg, grid, p, nodes, 1)

    located = output_mask(p, grid, 1)
    np.testing.assert_allclose(vorticity_number(p)[located], 0.0, atol=1e-12)
    np.testing.assert_allclose(rotation_angle(p)[located], 0.0, atol=1e-9)


def test_simple_shear_is_one_and_forty_five_degrees():
    """The reference point of the scale: a shear band."""
    cfg, grid, p, nodes = _case()
    _, y = _node_xy(grid)
    nodes.velocity[:, 0] = 0.04 * y
    _step(cfg, grid, p, nodes, 1)

    located = output_mask(p, grid, 1)
    np.testing.assert_allclose(vorticity_number(p)[located], 1.0, rtol=1e-12)
    np.testing.assert_allclose(rotation_angle(p)[located], 45.0, rtol=1e-12)


def test_rigid_rotation_is_undefined_and_ninety_degrees():
    """No deviatoric strain to divide by, which is exactly what makes it a rigid rotation."""
    cfg, grid, p, nodes = _case()
    x, y = _node_xy(grid)
    nodes.velocity[:, 0] = -0.03 * y
    nodes.velocity[:, 1] = 0.03 * x
    _step(cfg, grid, p, nodes, 1)

    located = output_mask(p, grid, 1)
    assert np.isnan(vorticity_number(p)[located]).all()      # the ratio has no value here
    np.testing.assert_allclose(rotation_angle(p)[located], 90.0, rtol=1e-12)


def test_a_particle_that_did_nothing_describes_nothing():
    cfg, grid, p, nodes = _case()
    _step(cfg, grid, p, nodes, 1)                            # the field is zero everywhere
    located = output_mask(p, grid, 1)
    assert np.isnan(vorticity_number(p)[located]).all()
    assert np.isnan(rotation_angle(p)[located]).all()


def test_the_angle_is_the_arctangent_of_the_number():
    """They are one quantity in two shapes, so Wm = tan(angle) has to hold."""
    cfg, grid, p, nodes = _case()
    _, y = _node_xy(grid)
    nodes.velocity[:, 0] = 0.04 * y
    nodes.velocity[:, 1] = 0.01 * y                          # something less tidy
    _step(cfg, grid, p, nodes, 1)

    located = output_mask(p, grid, 1)
    number = vorticity_number(p)[located]
    angle = rotation_angle(p)[located]
    np.testing.assert_allclose(np.tan(np.radians(angle)), number, rtol=1e-10)


def test_the_synthetic_rotation_case_is_called_a_rotation(workdir: Path):
    """The case we know is a rigid rotation has to be reported as one, not as shear."""
    work = workdir / "rotation"
    shutil.copytree(SYNTHETIC / "rotation_1P", work, dirs_exist_ok=True,
                    ignore=shutil.ignore_patterns("expected"))
    sim = Simulation.from_directory(work)
    sim.contour = make_contour_correction(3)
    sim.run()

    p = sim.particles
    alive = output_mask(p, sim.grid, sim.config.total_steps) & (p.nan_initial == 0)
    np.testing.assert_allclose(rotation_angle(p)[alive], 90.0, atol=1e-6)
    # the strain it reports is an apparent isotropic contraction, with no shear in it
    np.testing.assert_allclose(p.strain[alive, 0], p.strain[alive, 1], atol=1e-12)
    np.testing.assert_allclose(p.strain[alive, 2], 0.0, atol=1e-12)
    assert p.vol_strain[alive].mean() < -0.01      # and it claims a volume loss


# --- what gets written -------------------------------------------------------------------
def test_the_two_results_are_published():
    names = [spec.name for spec in result_specs()]
    assert "Vorticity" in names and "Rotation" in names
    assert "Vorticity_num" in names and "Rot_angle" in names


def test_compatibility_mode_leaves_them_out():
    """The original Fortran wrote no such block, and the regression suite compares bytes."""
    names = [spec.name for spec in result_specs(legacy_compat=True)]
    for new in ("Vorticity", "Rotation", "Vorticity_num", "Rot_angle"):
        assert new not in names


def test_they_reach_the_results_file(workdir: Path):
    work = workdir / "case"
    shutil.copytree(SYNTHETIC / "shear_1P", work, dirs_exist_ok=True,
                    ignore=shutil.ignore_patterns("expected"))
    run_case(work)
    text = (work / "zapatak.POST.RES").read_text()
    assert "Vorticity" in text and "Rotation" in text

    run_case(work, options=RunOptions(legacy_compat=True))
    legacy = (work / "zapatak.POST.RES").read_text()
    assert "Vorticity" not in legacy and "Rotation" not in legacy
