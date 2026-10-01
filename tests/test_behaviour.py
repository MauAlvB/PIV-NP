"""Expected behaviour of the analysis, checked against known solutions.

Each test checks a property that must hold (acceleration of a field with constant
acceleration, strain of a linear field, equivalence between a single run and a restarted
one...) and, where it makes sense, that the ``legacy_compat`` mode still reproduces the
behaviour of the original Fortran.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
import pytest

from pivnp.config import ConfigError, parse_par
from pivnp.gid_writer import GidWriter
from pivnp.mesh import particle_grid
from pivnp.nodal import load_measurements
from pivnp.particles import create_particles, local_coordinates
from pivnp.pivlab_io import Frame, pivlab_to_node
from pivnp.restart import read_restart, write_restart
from pivnp.simulation import RunOptions, Simulation, run_case
from pivnp.solver import output_mask
from pivnp.state import Nodes
from pivnp.vtk_export import iter_time_steps, read_gid_mesh

from .test_config import PAR

HEADER = "PIVlab\nFRAME\nx,y,u,v\n"


def _par(npc=2, version=1, steps=4, print_every=1, restart=0, n_cols=3, n_rows=2, size=1.0):
    return (f"case\nb2\n{n_cols * n_rows} {(n_cols + 1) * (n_rows + 1)} {npc} {n_rows} "
            f"{size} {size}\nb3\n0.5 {steps} {print_every} 0 {version} 1 0 {restart} 0\n"
            "b4\n2000 0.4\n")


def _write_frames(directory: Path, velocities, n_cols=3, n_rows=2, size=1.0):
    """Write PIVlab files with the velocity returned by ``velocities(step, x, y)``."""
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "PIV-NP.TXT").write_text("case\n")
    for step, uv in enumerate(velocities, start=1):
        rows = []
        for col in range(n_cols + 1):  # PIVlab: by columns, from top to bottom
            for row_from_top in range(n_rows + 1):
                x = col * size
                y = (n_rows - row_from_top) * size
                u, v = uv(x, y)
                rows.append(f"{x},{y},{u!r},{-v!r}")  # the v of PIVlab has the y axis downwards
        (directory / f"datos ({step}).txt").write_text(HEADER + "\n".join(rows) + "\n")


# --- acceleration ----------------------------------------------------------------------
def test_previous_velocity_is_the_whole_field():
    conn = pivlab_to_node(2, 2)
    nodes = Nodes.zeros(9, 9)
    rng = np.random.default_rng(0)
    first, second = rng.normal(size=9), rng.normal(size=9)
    zeros = np.zeros(9)
    for values in (first, second):
        load_measurements(Frame(1, None, values, zeros, zeros, zeros), conn, nodes)
    np.testing.assert_array_equal(nodes.velocity[conn, 0], second)
    np.testing.assert_array_equal(nodes.previous_velocity[conn, 0], first)


def test_acceleration_is_exact_for_every_node(workdir: Path):
    # u = a·t with a = 0.2 m/s²: the acceleration must be 0.2 on every particle.
    dt, accel = 0.5, 0.2
    _write_frames(workdir, [lambda x, y, s=s: (accel * s * dt, 0.0) for s in range(1, 5)],
                  n_cols=6, n_rows=6)
    (workdir / "case.PAR").write_text(_par(steps=4, n_cols=6, n_rows=6))
    sim = Simulation.from_directory(workdir)
    sim.run()
    located = output_mask(sim.particles, sim.grid, 4)
    np.testing.assert_allclose(sim.particles.acceleration[located, 0], accel)

    legacy = Simulation.from_directory(workdir, options=RunOptions(legacy_compat=True))
    legacy.run()
    assert not np.allclose(legacy.particles.acceleration[located, 0], accel)


# --- strains ------------------------------------------------------------
@pytest.mark.parametrize("npc", [2, 3])
@pytest.mark.parametrize("version", [1, 2])
def test_strain_matches_analytic_field(npc, version, workdir: Path):
    # u = rate·x  ->  eps_xx = rate·dt on each step, with any mesh and any NPC.
    rate, dt, steps = 0.01, 0.5, 2
    _write_frames(workdir, [lambda x, y: (rate * x, 0.0)] * steps, n_cols=6, n_rows=6)
    (workdir / "case.PAR").write_text(
        _par(npc=npc, version=version, steps=steps, n_cols=6, n_rows=6))
    sim = Simulation.from_directory(workdir)
    sim.run()
    p = sim.particles
    # Interior region: with IVERSION=2 the nodes on the outer boundary of the staggered
    # mesh receive fewer contributions, so the velocity there is not exact without a
    # contour correction.
    x, y = p.position[:, 0], p.position[:, 1]
    inner = (output_mask(p, sim.grid, steps)
             & (x > 1.0) & (x < 5.0) & (y > 1.0) & (y < 5.0))
    assert inner.sum() > 10
    np.testing.assert_allclose(p.strain[inner, 0], steps * rate * dt, rtol=1e-9)

    legacy = Simulation.from_directory(workdir, options=RunOptions(legacy_compat=True))
    legacy.run()
    # The original divides by the nodal mass: with IVERSION=1 it is 1 and changes nothing;
    # with IVERSION=2 it is roughly NPC² (not exactly, because the particles move).
    ratio = legacy.particles.strain[inner, 0] / p.strain[inner, 0]
    expected_ratio = 1.0 if version == 1 else 1 / npc**2
    np.testing.assert_allclose(ratio, expected_ratio, rtol=1e-9 if version == 1 else 0.05)


# --- rejected inputs ------------------------------------------------------------
@pytest.mark.parametrize("npc", [0, 7, 11])
def test_particles_per_side_limited_to_six(npc):
    with pytest.raises(ConfigError, match="NPC"):
        parse_par(PAR.replace("2006  2100    3", f"2006  2100    {npc}"))
    with pytest.raises(ValueError):
        local_coordinates(npc)


def test_six_particles_per_side_is_accepted():
    assert parse_par(PAR.replace("2006  2100    3", "2006  2100    6")).particles_per_side == 6


def test_tracking_particles_are_rejected():
    text = PAR.replace("0\t0\t0\n", "0\t0\t1\n") + "BLOQUE 5\n1 1 1 1 2 2 3 3\n"
    with pytest.raises(ConfigError, match="PTV"):
        parse_par(text)


# --- GiD mesh with the initial positions -----------------------------------------------
def test_mesh_plus_displacement_is_the_current_position(workdir: Path):
    _write_frames(workdir, [lambda x, y: (0.4, 0.2)] * 3)
    (workdir / "case.PAR").write_text(_par(steps=3, print_every=1))
    sim = Simulation.from_directory(workdir)
    sim.run()

    mesh = read_gid_mesh(workdir / "case.POST.MSH")
    steps = list(iter_time_steps(workdir / "case.POST.RES"))
    for _time, blocks in steps:
        disp = blocks["Displacement"]
        positions = mesh.coords[disp.ids - 1] + disp.values
        assert np.isfinite(positions).all()
    # the last instant must match the final position of the particles
    last = steps[-1][1]["Displacement"]
    np.testing.assert_allclose(mesh.coords[last.ids - 1] + last.values,
                               sim.particles.position[last.ids - 1], atol=1e-6)
    # and the mesh, the initial position (particles at half a cell, unmoved)
    expected_start = create_particles(sim.config, particle_grid(sim.config)).position
    np.testing.assert_allclose(mesh.coords, expected_start, atol=1e-6)


def test_legacy_mode_keeps_the_old_mesh(workdir: Path):
    _write_frames(workdir, [lambda x, y: (0.4, 0.2)] * 2)
    (workdir / "case.PAR").write_text(_par(steps=2))
    sim = Simulation.from_directory(workdir, options=RunOptions(legacy_compat=True))
    sim.run()
    mesh = read_gid_mesh(workdir / "case.POST.MSH")
    start = create_particles(sim.config, particle_grid(sim.config)).position
    np.testing.assert_allclose(mesh.coords, start + [0.2, 0.1], atol=1e-6)  # moved by one step


# --- restart with the staggered mesh ------------------------------------------------------------
@pytest.mark.parametrize("version", [1, 2])
def test_restart_continues_the_analysis(version, workdir: Path):
    """4 steps in a row == 2 steps + a restart with the next 2."""
    def field(step):
        return lambda x, y, s=step: (0.05 * s * (1 + 0.1 * y), 0.02 * s * x)

    full = workdir / "full"
    _write_frames(full, [field(s) for s in range(1, 5)], n_cols=4, n_rows=4)
    (full / "case.PAR").write_text(_par(version=version, steps=4, n_cols=4, n_rows=4))
    run_case(full)

    part = workdir / "partial"
    _write_frames(part, [field(s) for s in (1, 2)], n_cols=4, n_rows=4)
    (part / "case.PAR").write_text(_par(version=version, steps=2, n_cols=4, n_rows=4))
    run_case(part)
    _write_frames(part, [field(s) for s in (3, 4)], n_cols=4, n_rows=4)  # renumbered to 1 and 2
    (part / "case.PAR").write_text(_par(version=version, steps=2, restart=1, n_cols=4, n_rows=4))
    run_case(part)

    expected = read_restart(full / "case.REC")
    restarted = read_restart(part / "case.REC")
    np.testing.assert_allclose(restarted.position, expected.position, rtol=1e-12)
    np.testing.assert_allclose(restarted.displacement, expected.displacement, rtol=1e-12)
    np.testing.assert_allclose(restarted.strain, expected.strain, rtol=1e-12)


def test_legacy_restart_sends_everything_to_the_first_cell(workdir: Path):
    _write_frames(workdir, [lambda x, y: (0.1, 0.0)] * 2, n_cols=4, n_rows=4)
    (workdir / "case.PAR").write_text(_par(version=2, steps=2, n_cols=4, n_rows=4))
    options = RunOptions(legacy_compat=True)
    run_case(workdir, options=options)
    (workdir / "case.PAR").write_text(_par(version=2, steps=2, restart=1, n_cols=4, n_rows=4))
    sim = Simulation.from_directory(workdir, options=options)
    sim.run()
    moving = np.abs(sim.nodes.momentum[:, 0]) > 0
    assert moving.sum() == 4  # only the 4 nodes of cell 1

    sim_fixed = Simulation.from_directory(workdir)
    sim_fixed.run()
    assert (np.abs(sim_fixed.nodes.momentum[:, 0]) > 0).sum() > 4


# --- continuous restart -----------------------------------------------------------------
def _mini_case(directory: Path, steps: int, restart: bool, first_frame: int = 1) -> Path:
    """Case with the real data cut down, starting at instant ``first_frame``."""
    regression = Path(__file__).parent / "data" / "regression"
    (directory / "pivlab").mkdir(parents=True, exist_ok=True)
    for k in range(steps):
        shutil.copy(regression / "frames" / "pivlab" / f"datos ({first_frame + k}).txt",
                    directory / "pivlab" / f"datos ({k + 1}).txt")
    (directory / "PIV-NP.TXT").write_text("mini\n")
    par = (regression / "v1_npc3" / "mini.PAR").read_text().splitlines()
    par[4] = f"0.8 {steps} 1 0 1 1 0 {int(restart)} 0"
    (directory / "mini.PAR").write_text("\n".join(par) + "\n")
    return directory


def test_restart_produces_the_same_results_file_as_a_single_run(workdir: Path):
    full = _mini_case(workdir / "full", steps=8, restart=False)
    run_case(full)

    part = _mini_case(workdir / "partial", steps=4, restart=False)
    run_case(part)
    _mini_case(part, steps=4, restart=True, first_frame=5)  # instants 5..8 renumbered
    run_case(part)

    assert (part / "mini.POST.RES").read_bytes() == (full / "mini.POST.RES").read_bytes()
    assert (part / "mini.POST.MSH").read_bytes() == (full / "mini.POST.MSH").read_bytes()
    assert (part / "mini.REC").read_bytes() == (full / "mini.REC").read_bytes()


def test_instant_displacement_after_restart_is_only_the_step(workdir: Path):
    part = _mini_case(workdir / "case", steps=2, restart=False)
    run_case(part)
    _mini_case(part, steps=1, restart=True, first_frame=3)
    sim = Simulation.from_directory(part)
    sim.run()
    p = sim.particles
    active = output_mask(p, sim.grid, sim.config.total_steps) & (p.nan_initial == 0)
    assert np.abs(p.step_displacement[active]).max() < np.abs(p.displacement[active]).max()


def test_restart_file_keeps_the_original_fortran_records(workdir: Path):
    case = _mini_case(workdir / "case", steps=2, restart=False)
    run_case(case)
    extended = read_restart(case / "mini.REC")
    assert extended.step == 2 and extended.time == pytest.approx(1.6)
    assert extended.nodes is not None and extended.initial_position is not None

    # the extended file starts exactly with the 7 records the Fortran reads
    write_restart(workdir / "legacy_only.REC", extended, extended=False)
    write_restart(workdir / "extended.REC", extended, extended=True)
    legacy_bytes = (workdir / "legacy_only.REC").read_bytes()
    assert (workdir / "extended.REC").read_bytes().startswith(legacy_bytes)

    # and reading it without the extension gives the particle state, without the nodal one
    plain = read_restart(workdir / "legacy_only.REC")
    assert plain.nodes is None and plain.step == 0 and plain.time == 0.0
    np.testing.assert_array_equal(plain.position, extended.position)


# --- the "NaNs" result -------------------------------------------------------------------
@pytest.mark.parametrize(("version", "maximum"), [(1, 4), (2, 1)])
def test_missing_data_counts_what_is_missing(version, maximum, workdir: Path):
    regression = Path(__file__).parent / "data" / "regression"
    shutil.copytree(regression / "frames", workdir, dirs_exist_ok=True)
    par = (regression / "v1_npc3" / "mini.PAR").read_text().splitlines()
    par[4] = f"0.8 2 1 0 {version} 1 0 0 0"
    (workdir / "mini.PAR").write_text("\n".join(par) + "\n")
    (workdir / "PIV-NP.TXT").write_text("mini\n")

    sim = Simulation.from_directory(workdir)
    sim.run()
    counts = sim.particles.missing_data
    located = output_mask(sim.particles, sim.grid, 2)
    assert counts[located].max() == maximum
    assert counts[located].min() == 0
    # the particles marked from the first step are the ones with no data at all
    if version == 1:
        assert (counts[sim.particles.nan_initial == 1] == 4).all()


# --- kinetic energy --------------------------------------------------------------------
def test_kinetic_energy_is_a_single_scalar(workdir: Path):
    case = _mini_case(workdir / "case", steps=3, restart=False)
    sim = Simulation.from_directory(case)
    sim.run()

    blocks = list(iter_time_steps(case / "mini.POST.RES"))[-1][1]
    kinetic = blocks["E_kinetic"]
    assert kinetic.values.shape[1] == 1  # one value per line, as the header says
    p = sim.particles
    expected = p.kinetic_energy.sum(axis=1)[kinetic.ids - 1]
    np.testing.assert_allclose(kinetic.values[:, 0], expected, rtol=1e-5)
    # and now the sum of energies adds up
    total = blocks["E_total"].values[:, 0]
    potential = blocks["E_potential"].values[:, 0]
    np.testing.assert_allclose(total, potential + kinetic.values[:, 0], rtol=1e-5)


def test_legacy_mode_keeps_both_components(workdir: Path):
    case = _mini_case(workdir / "case", steps=2, restart=False)
    run_case(case, options=RunOptions(legacy_compat=True))
    kinetic = dict(list(iter_time_steps(case / "mini.POST.RES"))[-1][1])["E_kinetic"]
    assert kinetic.values.shape[1] == 2


# --- cross-check -------------------------------------------------------------------------
def test_fixed_and_legacy_modes_differ_only_where_expected(workdir: Path):
    """With IVERSION=1 and no restart only the acceleration, the mesh and the rounding change.

    The positions and strains differ only in the last digits, because of the constants the
    original kept in single precision.
    """
    regression = Path(__file__).parent / "data" / "regression"
    for name in ("PIV-NP.TXT", "mini.PAR"):
        shutil.copy(regression / "v1_npc3" / name, workdir)
    shutil.copytree(regression / "frames", workdir, dirs_exist_ok=True)

    fixed = Simulation.from_directory(workdir)
    fixed.run()
    legacy = Simulation.from_directory(workdir, options=RunOptions(legacy_compat=True))
    legacy.run()

    for field in ("position", "displacement", "strain", "eq_strain", "velocity"):
        np.testing.assert_allclose(getattr(fixed.particles, field),
                                   getattr(legacy.particles, field), rtol=1e-4, atol=1e-9)
    assert not np.array_equal(fixed.particles.acceleration, legacy.particles.acceleration)


def test_writer_still_produces_two_materials(workdir: Path):
    nan_initial = np.array([0, 1, 0], dtype=np.int8)
    writer = GidWriter(workdir, "case", eol="\n")
    writer.write_mesh(np.zeros((3, 2)), nan_initial)
    writer.close()
    text = (workdir / "case.POST.MSH").read_text()
    assert text.strip().endswith("End Elements")
    assert [line.split()[-1] for line in text.splitlines()[-4:-1]] == ["1", "2", "1"]
