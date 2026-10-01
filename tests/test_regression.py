"""Equivalence with the original Fortran code.

Every scenario in ``tests/data/regression`` was run with the original Fortran
(``tools/make_reference.py``). Here the Python version is run and the outputs
(``.POST.RES``, ``.POST.MSH`` and ``.REC``) are required to be identical byte for byte.
"""

from __future__ import annotations

import gzip
import json
import shutil
from pathlib import Path

import pytest

from pivnp.compare import compare_files
from pivnp.simulation import RunOptions, run_case

REGRESSION = Path(__file__).parent / "data" / "regression"
SCENARIOS = json.loads((REGRESSION / "scenarios.json").read_text(encoding="utf-8"))["scenarios"]
CASE = "mini"
# Compatibility mode: reproduces the behaviour of the original Fortran.
OPTIONS = RunOptions(eol="\r\n", prefetch=2, legacy_compat=True)


def _prepare(name: str, spec: dict, directory: Path) -> Path:
    frames = REGRESSION / ("frames5" if spec.get("pivlab_format", 1) != 1 else "frames")
    shutil.copytree(frames, directory, dirs_exist_ok=True)
    for file in ("PIV-NP.TXT", f"{CASE}.PAR"):
        shutil.copy(REGRESSION / name / file, directory)
    return directory


def _assert_same_text(expected_gz: Path, actual: Path, directory: Path) -> None:
    expected = directory / ("expected_" + actual.name)
    with gzip.open(expected_gz, "rb") as fin, open(expected, "wb") as fout:
        shutil.copyfileobj(fin, fout)
    if expected.read_bytes() != actual.read_bytes():
        report = compare_files(expected, actual)
        details = "\n".join(f"  line {n}:\n    ref: {a}\n    new: {b}"
                            for n, a, b in report.examples)
        pytest.fail(f"{actual.name} differs from the Fortran: {report.summary()}\n{details}")


@pytest.mark.parametrize("name", sorted(SCENARIOS))
def test_matches_legacy_fortran(name: str, workdir: Path) -> None:
    spec = SCENARIOS[name]
    expected = REGRESSION / name / "expected"
    work = _prepare(name, spec, workdir / "run")

    run_case(work, options=OPTIONS)
    if "restart_steps" in spec:
        assert (work / f"{CASE}.REC").read_bytes() == (expected / "phase1.REC").read_bytes()
        shutil.copy(REGRESSION / name / "restart.PAR", work / f"{CASE}.PAR")
        run_case(work, options=OPTIONS)

    for suffix in (".POST.MSH", ".POST.RES"):
        _assert_same_text(expected / f"{CASE}{suffix}.gz", work / f"{CASE}{suffix}", workdir)
    assert (work / f"{CASE}.REC").read_bytes() == (expected / f"{CASE}.REC").read_bytes()
