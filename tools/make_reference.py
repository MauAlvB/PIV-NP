"""Generate the reference results of the original Fortran code for the tests.

1. Builds ``legacy/*.for`` with gfortran (unless ``legacy/pivnp_legacy.exe`` is there).
2. Crops the PIVlab data of a real case down to a small window
   (``tests/data/regression/frames*``) and synthesises moisture files.
3. For every scenario in ``tests/data/regression/scenarios.json`` it writes the ``.PAR``,
   runs the Fortran and stores its outputs, compressed, in ``<scenario>/expected``.

Usage (from the root of the repository, with gfortran on the PATH)::

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
    print("Building:", " ".join(cmd))
    subprocess.run(cmd, check=True, cwd=LEGACY)
    return exe


def crop_frames(source: Path, window: dict, dest: Path, dest5: Path) -> None:
    """Crop every ``datos (n).txt`` to the window and build the five-column variant.

    The files go into the ``pivlab/`` and ``moisture/`` subfolders of each frame set, which
    is the layout of a case. ``stage_flat`` undoes that when the Fortran has to read them.
    """
    n_cols_src, n_rows_src = window["source_grid"]
    c0, c1 = window["cols"]
    r0, r1 = window["rows_from_top"]
    keep = [c * n_rows_src + r for c in range(c0, c1) for r in range(r0, r1)]
    for folder in (dest, dest5):
        (folder / "pivlab").mkdir(parents=True, exist_ok=True)
        (folder / "moisture").mkdir(parents=True, exist_ok=True)
    for step in range(1, window["n_frames"] + 1):
        lines = (source / f"datos ({step}).txt").read_text(encoding="latin-1").splitlines()
        assert len(lines) - 3 == n_cols_src * n_rows_src
        header, body = lines[:3], lines[3:]
        rows = [body[k] for k in keep]
        (dest / "pivlab" / f"datos ({step}).txt").write_text("\n".join(header + rows) + "\n")
        rows5 = [f"{row},1" for row in rows]
        (dest5 / "pivlab" / f"datos ({step}).txt").write_text("\n".join(header + rows5) + "\n")
        moist = dest / "moisture" / f"Moist_{step}.TXT"
        _write_moisture(moist, rows, step)
        shutil.copy(moist, dest5 / "moisture" / f"Moist_{step}.TXT")


def stage_flat(frames: Path, work: Path) -> None:
    """Copy a frame set into ``work`` with every file in one directory.

    The original Fortran reads ``datos (n).txt`` from its own working directory and knows
    nothing about subfolders, so the files have to be flattened before running it. The
    Python version reads either layout; this is only for the executable.
    """
    for path in frames.rglob("*"):
        if path.is_file():
            shutil.copy(path, work / path.name)


def _write_moisture(path: Path, rows: list[str], step: int) -> None:
    """Synthetic moisture: NaN wherever there is no velocity, plus a few negative values."""
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
    """The ``.PAR`` in the single format, the same one ``pivnp --convert-par`` writes."""
    from pivnp.par_migrate import ANALYSIS_NAMES, GEOMETRY_NAMES, SOIL_NAMES, block

    geometry = [str(n_cols * n_rows), str((n_cols + 1) * (n_rows + 1)), str(spec["npc"]),
                str(n_rows), str(size), str(size)]
    analysis = ["0.8", str(spec["restart_steps"] if restart else spec["steps"]),
                str(spec["print_every"]), str(int(spec.get("moisture", False))),
                str(spec["version"]), str(spec.get("pivlab_format", 1)), "0",
                str(int(restart)), "0"]
    lines = ["PIV-NP regression case"]
    lines += block("BLOCK 2:", GEOMETRY_NAMES, geometry)
    lines += block("BLOCK 3:", ANALYSIS_NAMES, analysis)
    lines += block("BLOCK 4:", SOIL_NAMES, ["2650.0", "0.4"])
    return "\n".join(lines) + "\n"


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
        stage_flat(frames, work)  # the Fortran cannot read the subfolders
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
                        help="directory holding the 'datos (n).txt' of the real case")
    parser.add_argument("--exe", type=Path, default=LEGACY / "pivnp_legacy.exe")
    parser.add_argument("--only", nargs="*", help="scenarios to regenerate")
    args = parser.parse_args(argv)

    spec = json.loads((REGRESSION / "scenarios.json").read_text(encoding="utf-8"))
    exe = build_legacy(args.exe.resolve())
    crop_frames(args.source, spec["window"], REGRESSION / "frames", REGRESSION / "frames5")
    print("Generating references:")
    for name, scenario in spec["scenarios"].items():
        if not args.only or name in args.only:
            make_scenario(name, scenario, exe, spec["window"], spec["cell_size"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
