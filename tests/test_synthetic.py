"""Validation with synthetic cases of known solution.

They are velocity fields built by hand on small meshes: displacement without strain,
uniform and variable horizontal strain, shear strain and rigid-body rotation. Each one has
an analytic solution, so they check the computation itself, not just that it has not
changed with respect to earlier versions.

The data and the historical results come from the group's analyses (`cases.json` records
the original name of each case) and their `.PAR` files are in the old format, so these
tests also check that those are still read.
"""

from __future__ import annotations

import gzip
import math
import shutil
from pathlib import Path

import numpy as np
import pytest

from pivnp.contour import make_contour_correction
from pivnp.particles import seed_positions
from pivnp.pivlab_io import FrameSource, pivlab_to_node
from pivnp.simulation import Simulation
from pivnp.solver import output_mask
from pivnp.vtk_export import iter_time_steps

SYNTHETIC = Path(__file__).parent / "data" / "synthetic"
CASES = sorted(d.name for d in SYNTHETIC.iterdir() if d.is_dir())
#: In the rotation every particle changes cell, so the cell-by-cell check does not apply;
#: that case has a test of its own.
CASES_WITHOUT_ROTATION = [c for c in CASES if not c.startswith("rotation")]


def run(name: str, destination: Path) -> Simulation:
    shutil.copytree(SYNTHETIC / name, destination, dirs_exist_ok=True,
                    ignore=shutil.ignore_patterns("expected"))
    sim = Simulation.from_directory(destination)
    sim.run()
    return sim


def active(sim: Simulation) -> np.ndarray:
    p = sim.particles
    return output_mask(p, sim.grid, sim.config.total_steps) & (p.nan_initial == 0)


def expected_strain(sim: Simulation) -> np.ndarray:
    """Accumulated strain of each cell, by finite differences of the velocities.

    This is a computation path independent of the one the program uses (derivatives of the
    shape functions), but equivalent for a bilinear field.
    """
    cfg, grid = sim.config, sim.grid
    conn = pivlab_to_node(cfg.n_cols, cfg.n_rows)
    source = FrameSource(sim.case_dir, cfg.n_nodes, cfg.pivlab_format, prefetch=0)
    cell_nodes = grid.cell_nodes(np.arange(grid.n_cells))  # (cells, 4)
    total = np.zeros((grid.n_cells, 3))
    for step in range(1, cfg.total_steps + 1):
        frame = source.read(step)
        u = np.zeros(cfg.n_nodes)
        v = np.zeros(cfg.n_nodes)
        u[conn] = np.nan_to_num(frame.u)
        v[conn] = -np.nan_to_num(frame.v)  # the y axis of the image points downwards
        u1, u2, u3, u4 = (u[cell_nodes[:, k]] for k in range(4))
        v1, v2, v3, v4 = (v[cell_nodes[:, k]] for k in range(4))
        dudx = ((u2 - u1) + (u4 - u3)) / (2 * grid.dx)
        dvdy = ((v3 - v1) + (v4 - v2)) / (2 * grid.dy)
        dudy = ((u3 - u1) + (u4 - u2)) / (2 * grid.dy)
        dvdx = ((v2 - v1) + (v4 - v3)) / (2 * grid.dx)
        total += np.stack([dudx, dvdy, dudy + dvdx], axis=1) * cfg.dt
    return total


@pytest.mark.parametrize("name", CASES_WITHOUT_ROTATION)
def test_strain_matches_the_analytic_solution(name: str, workdir: Path):
    """The particles that do not change cell accumulate the strain of that cell."""
    sim = run(name, workdir / name)
    p, grid = sim.particles, sim.grid
    start = seed_positions(grid, sim.config.particles_per_side)
    first_cell = grid.locate(start)
    last_cell = grid.locate(p.position)
    stayed = active(sim) & (first_cell == last_cell) & (first_cell >= 0)
    assert stayed.sum() >= 4

    expected = expected_strain(sim)[first_cell[stayed]]
    np.testing.assert_allclose(p.strain[stayed, :3], expected, atol=1e-12)


def test_displacement_without_strain(workdir: Path):
    """20 cm in x and −20 cm in y, the same on every particle and without deforming the solid."""
    for name in ("displacement_1P", "displacement_4P"):
        sim = run(name, workdir / name)
        p = sim.particles
        alive = active(sim)
        np.testing.assert_allclose(p.displacement[alive, 0], 0.20, atol=1e-12)
        np.testing.assert_allclose(p.displacement[alive, 1], -0.20, atol=1e-12)
        assert not p.strain[alive].any()
        assert not p.eq_strain[alive].any()


def test_uniform_shear_strain(workdir: Path):
    """Pure shear: γxy = 0.1 on every particle, with no normal or volumetric strain."""
    for name in ("shear_1P", "shear_4P"):
        sim = run(name, workdir / name)
        p = sim.particles
        alive = active(sim)
        np.testing.assert_allclose(p.strain[alive, 2], 0.10, atol=1e-12)
        assert not p.strain[alive, 0].any() and not p.strain[alive, 1].any()
        assert not p.vol_strain[alive].any()
        np.testing.assert_allclose(p.eq_strain[alive], 0.10 / math.sqrt(3.0), rtol=1e-12)


def test_uniform_horizontal_strain(workdir: Path):
    """Uniform stretching: εxx = 0.24 everywhere, εyy = γxy = 0, and εvol = εxx."""
    for name in ("horizontal_uniform_1P", "horizontal_uniform_4P"):
        sim = run(name, workdir / name)
        p = sim.particles
        alive = active(sim)
        np.testing.assert_allclose(p.strain[alive, 0], 0.24, atol=1e-12)
        assert not p.strain[alive, 1].any() and not p.strain[alive, 2].any()
        np.testing.assert_allclose(p.vol_strain[alive], 0.24, atol=1e-12)


def test_variable_horizontal_strain(workdir: Path):
    """Stretching that grows to the right: each column of cells has its own strain."""
    for name in ("horizontal_1P", "horizontal_4P"):
        sim = run(name, workdir / name)
        p = sim.particles
        alive = active(sim)
        values = np.unique(np.round(p.strain[alive, 0], 9))
        np.testing.assert_allclose(values[:3], [0.03, 0.06, 0.18], atol=1e-12)
        assert not p.strain[alive, 1].any() and not p.strain[alive, 2].any()


@pytest.mark.parametrize("name", ["rotation_1P", "rotation_4P"])
def test_rigid_body_rotation(name: str, workdir: Path):
    """A rigid rotation should not deform, but the incremental formulation does deform it.

    This is a known limitation of computing the strain from linear increments: the test
    pins the current error down so that any change making it worse is caught. Fixing it
    would need a finite-strain measure (the deformation gradient), not a local change.
    """
    sim = run(name, workdir / name)
    p = sim.particles
    alive = active(sim)
    assert np.abs(p.strain[alive, :3]).max() == pytest.approx(0.0439, abs=5e-4)
    assert np.abs(p.vol_strain[alive]).max() < 0.045
    # the rotation is symmetric about the centre, and so are the displacements
    np.testing.assert_allclose(p.displacement[alive, 0].min(), -p.displacement[alive, 0].max(),
                               rtol=1e-9)


def test_rotation_with_contour_correction_leaves_only_the_formulation_error(workdir: Path):
    """With the boundary corrected, exactly the theoretical formulation error is left.

    The rotation case has instants with 24 of its 49 nodes without data, so most of the
    apparent strain comes from the boundary. Rebuilding those nodes by extrapolation leaves
    every particle with the same apparent strain, which matches what the theory predicts:
    accumulating linear increments through a finite rotation leaves
    εxx = εyy = n·(cos Δθ − 1).
    """
    sim = run("rotation_1P", workdir / "rotation")
    sim_corrected = Simulation.from_directory(workdir / "rotation")
    sim_corrected.contour = make_contour_correction(3)
    sim_corrected.run()

    steps = sim.config.total_steps - 1  # the first instant does not rotate
    eps = steps * (math.cos(math.radians(1.0)) - 1.0)
    mean = 2 * eps / 3
    j2 = (2 * (eps - mean) ** 2 + mean**2) / 2
    theoretical = 2 * math.sqrt(3 * j2) / 3

    alive = active(sim_corrected)
    eq = sim_corrected.particles.eq_strain[alive]
    assert eq.max() - eq.min() < 1e-6  # uniform across the solid
    assert eq.mean() == pytest.approx(theoretical, rel=0.02)
    # uncorrected, the boundary multiplies the error by three
    assert sim.particles.eq_strain[active(sim)].mean() > 3 * theoretical


def test_the_strain_does_not_depend_on_the_particles_per_cell(workdir: Path):
    """With 1 and with 4 particles per cell, each cell gives the same strain."""
    for base in ("shear", "horizontal", "horizontal_uniform", "rotation"):
        values = []
        for suffix in ("1P", "4P"):
            sim = run(f"{base}_{suffix}", workdir / f"{base}_{suffix}")
            alive = active(sim)
            values.append(np.unique(np.round(sim.particles.strain[alive, :3], 9), axis=0))
        common = min(len(values[0]), len(values[1]))
        assert common > 0
        for row in values[0][:common]:
            assert any(np.allclose(row, other, atol=1e-9) for other in values[1]), base


@pytest.mark.parametrize("name", CASES)
def test_reproduces_the_historical_results(name: str, workdir: Path):
    """The analyses the group ran at the time are reproduced particle by particle."""
    sim = run(name, workdir / name)
    historical = workdir / f"{name}.POST.RES"
    with gzip.open(SYNTHETIC / name / "expected" / "zapatak.POST.RES.gz", "rb") as fin:
        historical.write_bytes(fin.read())

    old = {round(t, 3): b for t, b in iter_time_steps(historical)}
    new = {round(t, 3): b for t, b in iter_time_steps(sim.case_dir / "zapatak.POST.RES")}
    assert set(old) == set(new)

    last_time = max(old)
    for field in ("Displacement", "Total_strain", "Equi_strain"):
        a, b = old[last_time][field], new[last_time][field]
        common = np.intersect1d(a.ids, b.ids)
        assert common.size >= 16
        va = a.values[np.searchsorted(a.ids, common)]
        vb = b.values[np.searchsorted(b.ids, common)]
        np.testing.assert_allclose(vb, va, atol=1e-13, err_msg=f"{name}: {field}")
