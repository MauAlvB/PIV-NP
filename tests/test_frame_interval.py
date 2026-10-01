"""Warning when the DT of the .PAR does not match the interval it was exported from PIVlab."""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from pivnp.pivlab_io import frame_interval_in_header
from pivnp.simulation import Simulation

HEADER = ("PIVlab by W.Th. & E.J.S., ASCII chart output\n"
          "FRAME: 1, filenames: A: a.jpg & B: b.jpg, conversion factor xy (px -> m): "
          "{xy}, conversion factor uv (px/frame -> m/s): {uv}\n"
          "x [m],y [m],u [m/s],v [m/s]\n")

PAR = """case
b2
6 12 1 2 1.0 1.0
b3
{dt} 2 1 0 1 1 0 0 0
b4
2000 0.4
"""


def build(directory: Path, dt: float, xy: float, uv: float) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "PIV-NP.TXT").write_text("case\n")
    (directory / "case.PAR").write_text(PAR.format(dt=dt))
    for step in (1, 2):
        rows = "\n".join("0,0,0.001,0.0" for _ in range(12))
        (directory / f"datos ({step}).txt").write_text(
            HEADER.format(xy=xy, uv=uv) + rows + "\n")
    return directory


def test_reads_the_interval_from_the_header(workdir: Path):
    case = build(workdir, dt=0.02, xy=0.00030773, uv=0.015387)
    assert frame_interval_in_header(case / "datos (1).txt") == pytest.approx(0.02, rel=1e-4)


def test_without_factors_in_the_header(workdir: Path):
    (workdir / "d.txt").write_text("title\nno factors\nx,y,u,v\n")
    assert frame_interval_in_header(workdir / "d.txt") is None


def test_it_warns_when_the_dt_does_not_match(workdir: Path, caplog):
    # exported assuming 1 s between images, but the .PAR says 0.8 s
    case = build(workdir, dt=0.8, xy=0.0042423, uv=0.0042423)
    sim = Simulation.from_directory(case)
    with caplog.at_level(logging.WARNING, logger="pivnp"):
        interval = sim.check_frame_interval()
    assert interval == pytest.approx(1.0)
    assert "DT=0.8" in caplog.text and "0.8" in caplog.text


def test_it_does_not_warn_when_they_match(workdir: Path, caplog):
    case = build(workdir, dt=0.02, xy=0.00030773, uv=0.015387)
    sim = Simulation.from_directory(case)
    with caplog.at_level(logging.WARNING, logger="pivnp"):
        sim.check_frame_interval()
    assert not caplog.text
