"""Moisture computed from the images inside a full analysis (MOISTER=2).

A minimal case is set up: a one-cell grid, four nodes, three steps and uniform images that
get progressively darker, which is what the method reads as soil getting wet.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pivnp.config import ConfigError
from pivnp.moisture.source import Mesh, source_for_case
from pivnp.simulation import Simulation, write_moisture_files

HEADER = ("PIVlab\nConversion factor (px -> m): 0.001, (px/frame -> m/s): 0.001\n"
          "x [m],y [m],u [m/s],v [m/s]\n")
#: Nodes in metres. With 0.001 m/px they land on pixels (1,1), (1,4), (4,1) and (4,4).
NODES = [(0.002, 0.002), (0.002, 0.005), (0.005, 0.002), (0.005, 0.005)]
REFERENCE_GRAY = 200
#: Linear calibration table: normalized gray 0 (saturated) to 100 (dry).
CALIBRATION = "gray,saturation,moisture\n0,1,30\n100,0,0\n"
PAR = """\
Moisture-from-images test case
BLOQUE 2: N_cel N_nod N_part_celda N_fil Ancho Alto
          1     4     2            1     0.003 0.003
BLOQUE 3: del_t total_steps impresion moister version pivlab contour rec track
          1.0   3           1         {moister}       1       1      0       0   0
BLOQUE 4: s_density porosity
          2650.0    0.4
"""
HUM = """\
! moisture measurement settings
IMAGES = img_{n}.png
CHANNEL = gray
SIGMA = 1
DRY_REFERENCE = ref.png
CALIBRATION = soil.csv
"""


def build_case(directory: Path, moister: int = 2, grays=(199, 197, 194)) -> Path:
    """Write every file of the case and return its directory."""
    pytest.importorskip("PIL")
    from PIL import Image

    directory.mkdir(parents=True, exist_ok=True)
    (directory / "PIV-NP.TXT").write_text("'test'\n")
    (directory / "test.PAR").write_text(PAR.format(moister=moister))
    (directory / "test.HUM").write_text(HUM)
    (directory / "soil.csv").write_text(CALIBRATION)

    Image.fromarray(np.full((8, 8), REFERENCE_GRAY, dtype=np.uint8)).save(directory / "ref.png")
    for step, gray in enumerate(grays, 1):
        rows = "".join(f"{x},{y},0.0,0.0\n" for x, y in NODES)
        (directory / f"datos ({step}).txt").write_text(HEADER + rows)
        Image.fromarray(np.full((8, 8), gray, dtype=np.uint8)).save(directory / f"img_{step}.png")
    return directory


def normalized(gray: int) -> float:
    """Normalized gray in the dry-saturated band, which here is 11 levels wide."""
    return (gray - (REFERENCE_GRAY - 6)) / 11 * 100


def test_source_for_case_reads_the_settings_and_the_calibration(workdir: Path):
    case = build_case(workdir / "case")
    mesh = Mesh(np.array([x for x, _ in NODES]), np.array([y for _, y in NODES]), 0.001)
    source = source_for_case(case, "test", mesh, np.ones(4, dtype=bool))

    assert source.nodes_outside == 0
    assert source.settings.sigma == 1 and source.settings.channel == 0
    state = source.at_step(1, np.ones(4, dtype=bool))
    np.testing.assert_allclose(state.saturation, 1 - normalized(199) / 100, atol=1e-9)
    np.testing.assert_allclose(state.moisture, 30 * (1 - normalized(199) / 100), atol=1e-9)


def test_the_saturation_does_not_go_down_and_grows_as_it_darkens(workdir: Path):
    case = build_case(workdir / "case", grays=(199, 197, 199))
    (case / "test.HUM").write_text(HUM + "INCREMENTAL = 1\n")
    mesh = Mesh(np.array([x for x, _ in NODES]), np.array([y for _, y in NODES]), 0.001)
    source = source_for_case(case, "test", mesh, np.ones(4, dtype=bool))
    has_data = np.ones(4, dtype=bool)

    saturations = [source.at_step(step, has_data).saturation[0] for step in (1, 2, 3)]
    assert saturations[1] > saturations[0]
    assert saturations[2] == saturations[1]  # the third image brightens, but it does not undo


def test_full_analysis_with_moister_2(workdir: Path):
    case = build_case(workdir / "case")
    simulation = Simulation.from_directory(case)

    assert simulation.config.moisture and simulation.config.moisture_from_images
    assert simulation.frames.images is not None
    assert not simulation.frames.moisture  # no Moist_<n>.TXT are read: they do not exist

    simulation.run()
    assert (simulation.particles.moisture > 0).all()
    assert (simulation.nodes.saturation_measured > 0).all()
    assert (case / "test.POST.RES").exists()


def test_moister_2_without_the_hum_file(workdir: Path):
    case = build_case(workdir / "case")
    (case / "test.HUM").unlink()
    with pytest.raises(FileNotFoundError, match="HUM"):
        Simulation.from_directory(case)


def test_moister_2_without_the_calibration(workdir: Path):
    case = build_case(workdir / "case")
    (case / "soil.csv").unlink()
    with pytest.raises(FileNotFoundError):
        Simulation.from_directory(case)


def test_moister_1_still_reads_the_files(workdir: Path):
    """With MOISTER=1 no image is touched: moisture comes from the Moist_<n>.TXT."""
    case = build_case(workdir / "case", moister=1)
    for step in (1, 2, 3):
        rows = "".join(f"{x},{y},{step / 10},0.5\n" for x, y in NODES)
        (case / f"Moist_{step}.TXT").write_text("x,y,moisture,saturation\n" + rows)
    simulation = Simulation.from_directory(case)

    assert simulation.frames.images is None and simulation.frames.moisture
    simulation.run()
    assert (simulation.nodes.saturation_measured == 0.5).all()


def test_writing_the_moisture_files(workdir: Path):
    case = build_case(workdir / "case")
    assert write_moisture_files(case) == 3

    lines = (case / "Moist_2.TXT").read_text().splitlines()
    assert lines[0] == "x_m,y_m,moisture,saturation_degree"
    assert len(lines) == 1 + len(NODES)
    x, y, moisture, saturation = (float(v) for v in lines[1].split(","))
    assert (x, y) == NODES[0]
    assert saturation == pytest.approx(1 - normalized(197) / 100, abs=1e-9)
    assert moisture == pytest.approx(30 * saturation, abs=1e-9)

    # the written files then serve as input of an analysis with MOISTER=1
    (case / "test.PAR").write_text(PAR.format(moister=1))
    simulation = Simulation.from_directory(case)
    simulation.run()
    assert (simulation.particles.moisture > 0).all()


def test_the_analysis_says_how_much_of_it_is_a_measurement(workdir: Path):
    """With the gray below the band, the saturation is a bound and that has to be visible."""
    case = build_case(workdir / "case", grays=(199, 190, 185))  # the last two, past the band
    mesh = Mesh(np.array([x for x, _ in NODES]), np.array([y for _, y in NODES]), 0.001)
    source = source_for_case(case, "test", mesh, np.ones(4, dtype=bool))
    for step in (1, 2, 3):
        source.at_step(step, np.ones(4, dtype=bool))

    summary = source.quality_summary()
    assert "of 12 values with data" in summary
    assert "4 (33 %) measured" in summary
    assert "8 (67 %) at the wet limit" in summary


def test_unknown_moister(workdir: Path):
    case = build_case(workdir / "case", moister=5)
    with pytest.raises(ConfigError, match="MOISTER"):
        Simulation.from_directory(case)
