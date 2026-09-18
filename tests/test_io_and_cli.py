import gzip
import shutil
from pathlib import Path

import numpy as np
import pytest

from pivnp.cli import main
from pivnp.compare import compare_files
from pivnp.gid_writer import RESULTS, result_header
from pivnp.restart import RestartData, read_restart, write_restart

REGRESSION = Path(__file__).parent / "data" / "regression"


def test_restart_roundtrip(workdir: Path):
    rng = np.random.default_rng(0)
    data = RestartData(1, rng.normal(size=(5, 2)), rng.normal(size=(5, 2)),
                       rng.normal(size=(5, 4)), rng.normal(size=5),
                       np.array([0, 1, 0, 0, 1], dtype=np.int8))
    write_restart(workdir / "c.REC", data)
    back = read_restart(workdir / "c.REC")
    assert back.mesh_version == 1 and back.n_particles == 5
    np.testing.assert_array_equal(back.position, data.position)
    np.testing.assert_array_equal(back.strain, data.strain)
    np.testing.assert_array_equal(back.nan_initial, data.nan_initial)


def test_reads_restart_written_by_fortran():
    data = read_restart(REGRESSION / "v1_restart" / "expected" / "phase1.REC")
    assert data.mesh_version == 1
    assert data.position.shape == (99 * 4, 2)
    assert set(np.unique(data.nan_initial)) <= {0, 1}


def test_truncated_restart_raises(workdir: Path):
    path = workdir / "c.REC"
    path.write_bytes(b"\x04\x00\x00\x00\x05\x00\x00\x00\x04\x00\x00\x00")
    with pytest.raises(EOFError):
        read_restart(path)


def test_result_header_uses_a15_fields():
    assert result_header("Inst_displacement", "Vector", 0.8) == (
        "Result Inst_displaceme Isochrones   0.800000E+00 Vector         OnNodes")
    assert result_header("NaNs", "Scalar", 1.6).startswith("Result            NaNs")


def test_result_table_order_matches_legacy():
    names = [spec.name for spec in RESULTS]
    assert names[:4] == ["Displacement", "Inst_displacement", "NaNs", "Velocity"]
    assert names[-2:] == ["Moisture", "Saturation"]


def test_compare_tolerates_last_digit(workdir: Path):
    a, b = workdir / "a.RES", workdir / "b.RES"
    a.write_text("Values\n  1  0.123456E+00\n  2  0.100000E+01\n")
    b.write_text("Values\n  1  0.123457E+00\n  2  0.100002E+01\n")
    report = compare_files(a, b)
    assert (report.tolerated_lines, report.mismatched_lines) == (1, 1)
    assert not report.equivalent


def test_cli_runs_a_case(workdir: Path, capsys):
    shutil.copytree(REGRESSION / "frames", workdir, dirs_exist_ok=True)
    for name in ("PIV-NP.TXT", "mini.PAR"):
        shutil.copy(REGRESSION / "v1_npc3" / name, workdir)
    assert main([str(workdir), "--threads", "2", "-q"]) == 0
    expected = gzip.decompress((REGRESSION / "v1_npc3/expected/mini.POST.MSH.gz").read_bytes())
    produced = (workdir / "mini.POST.MSH").read_bytes()
    assert produced.replace(b"\r\n", b"\n") == expected.replace(b"\r\n", b"\n")


def test_cli_reports_bad_input(workdir: Path):
    (workdir / "PIV-NP.TXT").write_text("nada\n")
    assert main([str(workdir), "-q"]) == 1
