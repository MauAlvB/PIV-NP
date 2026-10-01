"""The displacement source is interchangeable.

The point of these tests is the last one: a source written from outside PIV-NP, reading no
files at all, drives a complete analysis. If that keeps working, adding a real source (a PIV
analysis built in, another package, a simulation) is writing one class.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pivnp.pivlab_io import PivlabSource
from pivnp.restart import read_restart
from pivnp.simulation import RunOptions, Simulation
from pivnp.solver import output_mask
from pivnp.sources import (
    SOURCES,
    DisplacementSource,
    Frame,
    available_sources,
    build_source,
    register_source,
)

# 3x2 cells of 1 m, 2x2 particles per cell, dt = 0.5, 2 steps
PAR = """a case with no files
block 2
6 12 2 2 1.0 1.0
block 3
0.5 2 1 0 1 1 0 0 0
block 4
2000 0.4
"""

U = 0.2  # m/s along x, the same on every point and every step


class ConstantSource:
    """A source with no files behind it: every point moves at the same speed."""

    def __init__(self, n_nodes: int) -> None:
        self.n_nodes = n_nodes
        self.moisture = False
        self.images = None
        self.asked: list[int] = []

    def frames(self, steps: range):
        for step in steps:
            self.asked.append(step)
            zeros = np.zeros(self.n_nodes)
            yield Frame(step, f"constant field {step}", np.full(self.n_nodes, U), zeros,
                        zeros, zeros)

    def mesh_in_metres(self, step: int = 1):
        x = np.zeros(self.n_nodes)
        return x, x, 0.001, np.ones(self.n_nodes, dtype=bool)

    def frame_interval(self) -> float | None:
        return None


@register_source("constant-for-tests")
def _build_constant(case_dir: Path, config, prefetch: int = 4) -> ConstantSource:
    return ConstantSource(config.n_nodes)


def _case(directory: Path) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "PIV-NP.TXT").write_text("nofiles\n")
    (directory / "nofiles.PAR").write_text(PAR)
    return directory


# --- the contract ----------------------------------------------------------------------
def test_the_shipped_source_satisfies_the_protocol(workdir: Path):
    # the protocol carries attributes (moisture, images), so it is checked on instances
    assert isinstance(PivlabSource(workdir, n_nodes=1), DisplacementSource)
    assert isinstance(ConstantSource(3), DisplacementSource)


def test_pivlab_is_available_by_name():
    assert "pivlab" in available_sources()
    assert available_sources() == sorted(SOURCES)


def test_an_unknown_source_says_which_ones_there_are():
    with pytest.raises(ValueError, match="unknown displacement source 'nope'"):
        build_source("nope", Path("."), None)


def test_frame_label_describes_any_origin():
    zeros = np.zeros(1)
    args = (zeros, zeros, zeros, zeros)
    assert Frame(1, Path("dir") / "datos (7).txt", *args).label == "datos (7).txt"
    assert Frame(2, "my own data", *args).label == "my own data"
    assert Frame(3, None, *args).label == "step 3"


def test_a_source_without_an_interval_is_not_checked(workdir: Path):
    """``frame_interval`` returning None has to silence the DT check, not crash it."""
    sim = Simulation.from_directory(_case(workdir / "case"),
                                    options=RunOptions(source="constant-for-tests"))
    assert sim.check_frame_interval() is None


def test_pivlab_interval_is_none_without_files(workdir: Path):
    source = PivlabSource(workdir, n_nodes=1)
    assert source.frame_interval() is None


# --- the point -------------------------------------------------------------------------
def test_a_source_from_outside_drives_a_whole_analysis(workdir: Path):
    """No PIVlab file exists here: every velocity comes from ConstantSource."""
    case = _case(workdir / "case")
    sim = Simulation.from_directory(case, options=RunOptions(source="constant-for-tests"))
    assert isinstance(sim.frames, ConstantSource)

    summary = sim.run()

    p = sim.particles
    located = output_mask(p, sim.grid, sim.config.total_steps)
    expected = U * sim.config.dt * sim.config.total_steps  # 0.2 m/s * 0.5 s * 2 steps
    np.testing.assert_allclose(p.displacement[located, 0], expected)
    np.testing.assert_allclose(p.displacement[located, 1], 0.0)
    assert not p.strain[located].any()  # a uniform field does not deform

    # the steps were asked for in order, and the results were written
    assert sim.frames.asked == [1, 2]
    assert summary.steps == 2
    assert (case / "nofiles.POST.RES").exists()
    assert (case / "nofiles.POST.MSH").exists()
    assert read_restart(case / "nofiles.REC").step == 2
