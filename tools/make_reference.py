"""Genera los resultados de referencia del código Fortran original para las pruebas.

1. Compila ``legacy/*.for`` con gfortran (si no existe ``legacy/pivnp_legacy.exe``).
2. Recorta los datos PIVlab de un caso real a una ventana pequeña
   (``tests/data/regression/frames*``) y sintetiza archivos de humedad.
3. Para cada escenario de ``tests/data/regression/scenarios.json`` escribe el ``.PAR``,
   ejecuta el Fortran y guarda sus salidas comprimidas en ``<escenario>/expected``.

Uso (desde la raíz del repositorio, con gfortran en el PATH)::

    python tools/make_reference.py --source "../caso paper centrifuga"
"""

from __future__ import annotations

import argparse
import gzip
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
LEGACY = ROOT / "legacy"
REGRESSION = ROOT / "tests" / "data" / "regression"
CASE = "mini"
GFORTRAN_FLAGS = ["-O2", "-finit-local-zero", "-ffixed-line-length-none", "-Wno-tabs", "-static"]


def build_legacy(exe: Path) -> Path:
    if exe.exists():
        return exe
    sources = [LEGACY / "MainCodePIV-NP.for", LEGACY / "contour_stub.for"]
    cmd = ["gfortran", *GFORTRAN_FLAGS, "-o", str(exe), *map(str, sources)]
    print("Compilando:", " ".join(cmd))
    subprocess.run(cmd, check=True, cwd=LEGACY)
    return exe


def crop_frames(source: Path, window: dict, dest: Path, dest5: Path) -> None:
    """Recorta cada ``datos (n).txt`` a la ventana y crea la variante de 5 columnas."""
    n_cols_src, n_rows_src = window["source_grid"]
    c0, c1 = window["cols"]
    r0, r1 = window["rows_from_top"]
    keep = [c * n_rows_src + r for c in range(c0, c1) for r in range(r0, r1)]
    dest.mkdir(parents=True, exist_ok=True)
    dest5.mkdir(parents=True, exist_ok=True)
    for step in range(1, window["n_frames"] + 1):
        lines = (source / f"datos ({step}).txt").read_text(encoding="latin-1").splitlines()
        assert len(lines) - 3 == n_cols_src * n_rows_src
        header, body = lines[:3], lines[3:]
        rows = [body[k] for k in keep]
        (dest / f"datos ({step}).txt").write_text("\n".join(header + rows) + "\n")
        rows5 = [f"{row},1" for row in rows]
        (dest5 / f"datos ({step}).txt").write_text("\n".join(header + rows5) + "\n")
        _write_moisture(dest / f"Moist_{step}.TXT", rows, step)
        shutil.copy(dest / f"Moist_{step}.TXT", dest5 / f"Moist_{step}.TXT")


def _write_moisture(path: Path, rows: list[str], step: int) -> None:
    """Humedad sintética: NaN donde no hay velocidad y algún valor negativo."""
    out = ["Moisture synthetic data"]
    for k, row in enumerate(rows):
        x, y, u, _ = row.split(",")[:4]
        if u.strip().lower() == "nan":
            moist, sat = "NaN", "NaN"
        else:
            moist = f"{0.05 + 0.01 * np.sin(k + step) - (0.08 if k % 17 == 0 else 0):.6f}"
            sat = f"{0.5 + 0.3 * np.cos(0.3 * k - step):.6f}"
        out.append(f"{x},{y},{moist},{sat}")
    path.write_text("\n".join(out) + "\n")


def par_text(spec: dict, n_cols: int, n_rows: int, size: float, restart: bool) -> str:
    """El ``.PAR`` en el formato único, el mismo que escribe ``pivnp --convert-par``."""
    from pivnp.par_migrate import ANALYSIS_NAMES, GEOMETRY_NAMES, SOIL_NAMES, bloque

    geometria = [str(n_cols * n_rows), str((n_cols + 1) * (n_rows + 1)), str(spec["npc"]),
                 str(n_rows), str(size), str(size)]
    analisis = ["0.8", str(spec["restart_steps"] if restart else spec["steps"]),
                str(spec["print_every"]), str(int(spec.get("moisture", False))),
                str(spec["version"]), str(spec.get("pivlab_format", 1)), "0",
                str(int(restart)), "0"]
    lineas = ["Caso de regresion PIV-NP"]
    lineas += bloque("BLOQUE 2:", GEOMETRY_NAMES, geometria)
    lineas += bloque("BLOQUE 3:", ANALYSIS_NAMES, analisis)
    lineas += bloque("BLOQUE 4:", SOIL_NAMES, ["2650.0", "0.4"])
    return "\n".join(lineas) + "\n"


def run_legacy(exe: Path, workdir: Path) -> None:
    subprocess.run([str(exe)], check=True, cwd=workdir, stdout=subprocess.DEVNULL)


def gzip_copy(src: Path, dest: Path) -> None:
    with open(src, "rb") as fin, gzip.GzipFile(dest, "wb", mtime=0) as fout:
        shutil.copyfileobj(fin, fout)


def make_scenario(name: str, spec: dict, exe: Path, window: dict, size: float) -> None:
    n_cols = window["cols"][1] - window["cols"][0] - 1
    n_rows = window["rows_from_top"][1] - window["rows_from_top"][0] - 1
    scenario_dir = REGRESSION / name
    expected = scenario_dir / "expected"
    shutil.rmtree(scenario_dir, ignore_errors=True)
    expected.mkdir(parents=True)
    (scenario_dir / "PIV-NP.TXT").write_text(f"{CASE}\n")
    (scenario_dir / f"{CASE}.PAR").write_text(par_text(spec, n_cols, n_rows, size, False))
    if "restart_steps" in spec:
        (scenario_dir / "restart.PAR").write_text(par_text(spec, n_cols, n_rows, size, True))

    frames = REGRESSION / ("frames5" if spec.get("pivlab_format", 1) != 1 else "frames")
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        shutil.copytree(frames, work, dirs_exist_ok=True)
        shutil.copy(scenario_dir / "PIV-NP.TXT", work)
        shutil.copy(scenario_dir / f"{CASE}.PAR", work)
        run_legacy(exe, work)
        if "restart_steps" in spec:
            shutil.copy(work / f"{CASE}.REC", expected / "phase1.REC")
            shutil.copy(scenario_dir / "restart.PAR", work / f"{CASE}.PAR")
            run_legacy(exe, work)
        for suffix in (".POST.RES", ".POST.MSH"):
            gzip_copy(work / f"{CASE}{suffix}", expected / f"{CASE}{suffix}.gz")
        shutil.copy(work / f"{CASE}.REC", expected / f"{CASE}.REC")
    print(f"  {name}: ok")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--source", type=Path, required=True,
                        help="directorio con los 'datos (n).txt' del caso real")
    parser.add_argument("--exe", type=Path, default=LEGACY / "pivnp_legacy.exe")
    parser.add_argument("--only", nargs="*", help="escenarios a regenerar")
    args = parser.parse_args(argv)

    spec = json.loads((REGRESSION / "scenarios.json").read_text(encoding="utf-8"))
    exe = build_legacy(args.exe.resolve())
    crop_frames(args.source, spec["window"], REGRESSION / "frames", REGRESSION / "frames5")
    print("Generando referencias:")
    for name, scenario in spec["scenarios"].items():
        if not args.only or name in args.only:
            make_scenario(name, scenario, exe, spec["window"], spec["cell_size"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
