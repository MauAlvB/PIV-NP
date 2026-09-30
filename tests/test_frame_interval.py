"""Aviso cuando el DT del .PAR no coincide con el intervalo con que se exportó de PIVlab."""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from pivnp.pivlab_io import frame_interval_in_header
from pivnp.simulation import Simulation

CABECERA = ("PIVlab by W.Th. & E.J.S., ASCII chart output\n"
            "FRAME: 1, filenames: A: a.jpg & B: b.jpg, conversion factor xy (px -> m): "
            "{xy}, conversion factor uv (px/frame -> m/s): {uv}\n"
            "x [m],y [m],u [m/s],v [m/s]\n")

PAR = """caso
b2
6 12 1 2 1.0 1.0
b3
{dt} 2 1 0 1 1 0 0 0
b4
2000 0.4
"""


def preparar(directorio: Path, dt: float, xy: float, uv: float) -> Path:
    directorio.mkdir(parents=True, exist_ok=True)
    (directorio / "PIV-NP.TXT").write_text("caso\n")
    (directorio / "caso.PAR").write_text(PAR.format(dt=dt))
    for paso in (1, 2):
        filas = "\n".join("0,0,0.001,0.0" for _ in range(12))
        (directorio / f"datos ({paso}).txt").write_text(
            CABECERA.format(xy=xy, uv=uv) + filas + "\n")
    return directorio


def test_lee_el_intervalo_de_la_cabecera(workdir: Path):
    caso = preparar(workdir, dt=0.02, xy=0.00030773, uv=0.015387)
    assert frame_interval_in_header(caso / "datos (1).txt") == pytest.approx(0.02, rel=1e-4)


def test_sin_factores_en_la_cabecera(workdir: Path):
    (workdir / "d.txt").write_text("titulo\nsin factores\nx,y,u,v\n")
    assert frame_interval_in_header(workdir / "d.txt") is None


def test_avisa_si_el_dt_no_coincide(workdir: Path, caplog):
    # exportado suponiendo 1 s entre imágenes, pero el .PAR dice 0.8 s
    caso = preparar(workdir, dt=0.8, xy=0.0042423, uv=0.0042423)
    sim = Simulation.from_directory(caso)
    with caplog.at_level(logging.WARNING, logger="pivnp"):
        intervalo = sim.check_frame_interval()
    assert intervalo == pytest.approx(1.0)
    assert "DT=0.8" in caplog.text and "0.8" in caplog.text


def test_no_avisa_si_coinciden(workdir: Path, caplog):
    caso = preparar(workdir, dt=0.02, xy=0.00030773, uv=0.015387)
    sim = Simulation.from_directory(caso)
    with caplog.at_level(logging.WARNING, logger="pivnp"):
        sim.check_frame_interval()
    assert not caplog.text
