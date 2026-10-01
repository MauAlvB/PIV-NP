"""The shipped example runs, and gives the numbers the guide promises.

``docs/GUIDE.md`` tells a first-time reader which values to expect after running
``examples/shear-block``. These tests pin those values down, so that the guide cannot go
stale without something failing here.
"""

from __future__ import annotations

import math
import shutil
from pathlib import Path

import numpy as np
import pytest

from pivnp.par_migrate import ANALYSIS_NAMES, GEOMETRY_NAMES, LEGEND, SOIL_NAMES, block
from pivnp.simulation import Simulation
from pivnp.solver import output_mask

EXAMPLE = Path(__file__).parents[1] / "examples" / "shear-block"
CASE = "shearblock"

#: Imposed field: u = SHEAR_RATE * height above the base, for TOTAL_TIME seconds.
SHEAR_RATE, TOTAL_TIME = 0.05, 2.0
SHEAR_STRAIN = SHEAR_RATE * TOTAL_TIME                      # 0.10
EQ_STRAIN = SHEAR_STRAIN / math.sqrt(3.0)                   # 0.057735


def _run(directory: Path, moisture: int = 0) -> Simulation:
    """Copy the example, set the moisture mode in its .PAR, and run it."""
    shutil.copytree(EXAMPLE, directory, dirs_exist_ok=True)
    lines = ["Shear block: a synthetic example with a known solution"]
    lines += block("BLOCK 2:", GEOMETRY_NAMES, ["96", "117", "2", "8", "0.01", "0.01"])
    lines += block("BLOCK 3:", ANALYSIS_NAMES,
                   ["0.1", "20", "1", str(moisture), "1", "1", "0", "0", "0"])
    lines += block("BLOCK 4:", SOIL_NAMES, ["2650.0", "0.4"])
    (directory / f"{CASE}.PAR").write_text("\n".join(lines) + "\n" + LEGEND, encoding="ascii")

    sim = Simulation.from_directory(directory)
    sim.run()
    return sim


def _published(sim: Simulation) -> np.ndarray:
    return output_mask(sim.particles, sim.grid, sim.config.total_steps)


def test_the_example_ships_every_file_it_needs():
    names = {p.name for p in EXAMPLE.iterdir()}
    assert {"PIV-NP.TXT", f"{CASE}.PAR", f"{CASE}.HUM", "calibration.csv",
            "make_example.py", "README.md"} <= names
    assert (EXAMPLE / "datos (1).txt").exists() and (EXAMPLE / "datos (20).txt").exists()
    assert len(list((EXAMPLE / "images").glob("wet_*.png"))) == 20


def test_the_example_gives_the_exact_solution(workdir: Path):
    """Simple shear has an analytic answer, and the example reproduces it exactly."""
    sim = _run(workdir / "case")
    p = sim.particles
    published = _published(sim)

    assert p.lost.size == 384                       # 96 cells x 2x2 particles
    assert published.sum() == 372                   # 12 leave the grid through the right
    assert not (p.nan_initial != 0).any()           # every point was measured

    np.testing.assert_allclose(p.strain[published, 2], SHEAR_STRAIN, atol=1e-12)
    np.testing.assert_allclose(p.eq_strain[published], EQ_STRAIN, atol=1e-12)
    assert not p.strain[published, 0].any()         # no horizontal strain
    assert not p.strain[published, 1].any()         # no vertical strain
    assert not p.vol_strain[published].any()        # shear does not change the volume
    assert not p.displacement[published, 1].any()   # nothing moves vertically

    # displacement = SHEAR_RATE * height * time, so it grows with height in 16 steps
    dx = p.displacement[published, 0]
    assert dx.min() == pytest.approx(0.00025)       # 0.25 mm, the lowest row of particles
    assert dx.max() == pytest.approx(0.00775)       # 7.75 mm, the highest one
    assert len(np.unique(np.round(dx, 12))) == 16   # 8 rows of cells x 2 particle rows


def test_turning_the_moisture_on_adds_fields_without_touching_the_mechanics(workdir: Path):
    """``moisture = 2`` in the .PAR is the whole switch, and it changes nothing else."""
    plain = _run(workdir / "plain", moisture=0)
    wet = _run(workdir / "wet", moisture=2)

    for field in ("position", "displacement", "strain", "eq_strain", "velocity"):
        np.testing.assert_array_equal(getattr(plain.particles, field),
                                      getattr(wet.particles, field))

    published = _published(wet)
    assert not plain.particles.saturation.any()     # without it, the fields stay at zero
    saturation = wet.particles.saturation[published]
    moisture = wet.particles.moisture[published]
    # the wetting front has climbed the block: the base is saturated, the top only partly
    assert saturation.max() == pytest.approx(1.0, abs=1e-9)
    assert 0.5 < saturation.min() < 0.6
    assert moisture.max() == pytest.approx(25.0, abs=1e-9)   # the wet end of the table
    # the cubic can land one ulp past 1.0 where the exact answer is 1.0
    assert (saturation >= 0).all() and (saturation <= 1.0 + 1e-12).all()


def test_every_moisture_value_of_the_example_is_a_measurement(workdir: Path):
    """The guide shows a quality line of "100 % measured", so it has to stay that way.

    Its calibration reaches one row past each end of the range the test produces, which is
    what keeps a node sitting exactly at 0 or 100 from being reported as a mere bound.
    """
    wet = _run(workdir / "wet", moisture=2)
    summary = wet.frames.images.quality_summary()
    assert "100 %) measured" in summary, summary
    assert "limit" not in summary, summary


def test_the_example_needs_no_boundary_correction(workdir: Path):
    """Its data has no gaps, so contour 0 and contour 1 give the same thing."""
    sim = _run(workdir / "case")
    assert not sim.nodes.filled.any()
