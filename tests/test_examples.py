"""The shipped examples run, and give the numbers their documentation promises.

``docs/GUIDE.md`` and each example's own ``README.md`` tell a first-time reader which values
to expect. These tests pin those values down, so that the documentation cannot go stale
without something failing here.

Two examples, and the difference between them is the point:

* ``examples/shear-block`` is driven by PIVlab files and has an **exact** answer, so it is
  checked to twelve decimals;
* ``examples/piv-from-images`` is driven by the photographs, through the built-in PIV, and
  has a *known* answer rather than an exact one -- a correlation measures to about a tenth
  of a pixel. It is checked to the accuracy that was measured, which is about one per cent
  on the shear.
"""

from __future__ import annotations

import math
import shutil
from pathlib import Path

import numpy as np
import pytest

from pivnp.par_migrate import ANALYSIS_NAMES, GEOMETRY_NAMES, LEGEND, SOIL_NAMES, block
from pivnp.simulation import RunOptions, Simulation
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
    # the inputs live in their own subfolders, so the case root stays readable
    assert len(list((EXAMPLE / "pivlab").glob("datos (*).txt"))) == 20
    assert len(list((EXAMPLE / "images").glob("wet_*.png"))) == 20
    assert not list(EXAMPLE.glob("datos (*).txt"))


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


# --- the example that needs no PIVlab -------------------------------------------------
PHOTOS = Path(__file__).parents[1] / "examples" / "piv-from-images"
PHOTO_CASE = "shearphotos"

#: What ``make_example.py`` put into the photographs: the top of a 128 px block slides 2 px
#: further per step, over 10 steps.
PHOTO_SHEAR = 10 * 2.0 / 128.0                              # 0.15625
PHOTO_EQ = PHOTO_SHEAR / math.sqrt(3.0)                     # 0.090211


def _run_from_photos(directory: Path, smoothing: float | None = None,
                     subpixel_offset: bool = False) -> Simulation:
    """Copy the example and measure its displacements from its photographs."""
    shutil.copytree(PHOTOS, directory, dirs_exist_ok=True)
    settings = directory / f"{PHOTO_CASE}.PIV"
    if smoothing is not None:
        text = settings.read_text(encoding="latin-1")
        settings.write_text(text.replace("SMOOTH = 0.6", f"SMOOTH = {smoothing}"),
                            encoding="latin-1")
    if subpixel_offset:
        with settings.open("a", encoding="latin-1") as handle:
            handle.write("\nSUBPIXEL_OFFSET = 1\n")
    sim = Simulation.from_directory(directory,
                                    options=RunOptions(source="images", prefetch=0))
    sim.run()
    return sim


def test_the_photograph_example_ships_every_file_it_needs():
    names = {p.name for p in PHOTOS.iterdir()}
    assert {"PIV-NP.TXT", f"{PHOTO_CASE}.PAR", f"{PHOTO_CASE}.PIV",
            "make_example.py", "README.md"} <= names
    assert len(list((PHOTOS / "images").glob("shear_*.png"))) == 11   # one more than steps
    # the whole point of this example: there is no PIVlab export anywhere in it
    assert not (PHOTOS / "pivlab").exists()
    assert not list(PHOTOS.glob("datos (*).txt"))


def test_the_photograph_example_recovers_the_shear_that_was_put_in(workdir: Path):
    """The headline claim of the built-in PIV, on a case whose answer is known.

    The tolerances are the measured accuracy, not an aspiration: the shear comes back 1.1 %
    low and the equivalent strain 3.2 % high, so the limits are set a little outside that.
    A change that made the correlation worse would fail here before anyone noticed it on
    real data.
    """
    sim = _run_from_photos(workdir / "photos")
    p = sim.particles
    published = output_mask(p, sim.grid, sim.config.total_steps) & (p.nan_initial == 0)

    assert p.lost.size == 1500                      # 375 cells x 2x2 particles
    assert published.sum() == 1440                  # 60 leave the grid as the block shears

    assert np.median(p.strain[published, 2]) == pytest.approx(PHOTO_SHEAR, rel=0.03)
    assert np.median(p.eq_strain[published]) == pytest.approx(PHOTO_EQ, rel=0.05)

    # a simple shear changes no length and nothing moves vertically, so these are the
    # artifact, and they are small rather than absent -- which is what the README says
    assert abs(np.median(p.strain[published, 0])) < 0.01
    assert abs(np.median(p.strain[published, 1])) < 0.01
    assert np.median(np.abs(p.displacement[published, 1])) < 2e-4     # under 0.2 mm

    # the top of the block slides and the base does not
    height = p.position[published, 1]
    slide = p.displacement[published, 0]
    top = height > np.percentile(height, 90)
    base = height < np.percentile(height, 10)
    assert np.median(slide[top]) > 0.008            # over 8 mm of the 10 imposed
    assert abs(np.median(slide[base])) < 0.002      # under 2 mm at the base


def test_reading_between_pixels_still_accumulates_worse_than_rounding(workdir: Path):
    """The open question, pinned so that solving it cannot go unnoticed.

    Reading each window between the pixels of the photograph gives a better field by every
    measure taken of a single step: four to five times less bias and scatter on known
    shifts, and peak locking gone. It also makes the *accumulated* shear of this example
    worse -- about 8 % low against about 1 % -- while making the artifacts that should be
    zero several times smaller. Both estimators read the gradient of each single step to
    within 1.3 %, so the loss happens over the ten steps of accumulation and nobody has
    found where.

    That is why ``SUBPIXEL_OFFSET`` is off by default. This test fails the day the
    accumulation stops losing it, which is exactly when the default should change -- so read
    the failure as the answer arriving, not as a break.
    """
    rounded = _run_from_photos(workdir / "rounded")
    between = _run_from_photos(workdir / "between", subpixel_offset=True)

    def shear(sim: Simulation) -> float:
        p = sim.particles
        kept = output_mask(p, sim.grid, sim.config.total_steps) & (p.nan_initial == 0)
        return float(np.median(p.strain[kept, 2]))

    off_rounded = abs(shear(rounded) - PHOTO_SHEAR) / PHOTO_SHEAR
    off_between = abs(shear(between) - PHOTO_SHEAR) / PHOTO_SHEAR
    assert off_rounded < 0.03, "the shipped default should stay within 3 % of the truth"
    assert off_between > 0.05, (
        "reading between pixels now accumulates within 5 % of the truth, which was the "
        "thing stopping it being the default -- re-measure and consider turning it on"
    )


def test_smoothing_quietens_the_strain_without_moving_the_displacement(workdir: Path):
    """Why ``SMOOTH`` defaults to 0.6, measured on the example rather than argued.

    The displacement is what the smoothing must not touch, and the strain is what it is
    for. Both are checked in one test because either alone would be easy to satisfy: a
    smoother that did nothing would pass the first, and one that flattened everything would
    pass the second.
    """
    raw = _run_from_photos(workdir / "raw", smoothing=0.0)
    smoothed = _run_from_photos(workdir / "smoothed", smoothing=0.6)

    def summary(sim: Simulation) -> tuple[float, float]:
        p = sim.particles
        kept = output_mask(p, sim.grid, sim.config.total_steps) & (p.nan_initial == 0)
        shear = p.strain[kept, 2]
        height = p.position[kept, 1]
        top = height > np.percentile(height, 90)
        spread = np.percentile(shear, 90) - np.percentile(shear, 10)
        return float(np.median(p.displacement[kept, 0][top])), float(spread)

    raw_slide, raw_spread = summary(raw)
    smooth_slide, smooth_spread = summary(smoothed)

    # the displacement is left alone: measured 9.275 mm against 9.274
    assert smooth_slide == pytest.approx(raw_slide, abs=5e-5)
    # and the scatter of the strain is cut by a third: measured 0.096 to 0.063
    assert smooth_spread < 0.8 * raw_spread
