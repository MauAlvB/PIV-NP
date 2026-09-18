import base64
import gzip
import xml.etree.ElementTree as ET
import zlib
from pathlib import Path

import numpy as np
import pytest

from pivnp.vtk_export import export_vtk, iter_time_steps, read_gid_mesh

EXPECTED = Path(__file__).parent / "data" / "regression" / "v1_npc3" / "expected"


def _decode(element) -> np.ndarray:
    """Decodifica un DataArray binario comprimido (UInt32, un bloque zlib)."""
    text = element.text.strip()
    header = np.frombuffer(base64.b64decode(text[:24]), "<u4")  # 16 bytes -> 24 chars
    raw = zlib.decompress(base64.b64decode(text[24:]))
    assert header[0] == 1 and header[1] == len(raw)
    dtype = {"Float32": "<f4", "Int32": "<i4", "UInt8": "u1"}[element.get("type")]
    data = np.frombuffer(raw, dtype)
    n_comp = int(element.get("NumberOfComponents"))
    return data.reshape(-1, n_comp) if n_comp > 1 else data


@pytest.fixture
def gid_case(workdir: Path) -> Path:
    for suffix in (".POST.RES", ".POST.MSH"):
        (workdir / f"mini{suffix}").write_bytes(gzip.decompress(
            (EXPECTED / f"mini{suffix}.gz").read_bytes()))
    return workdir / "mini.POST.RES"


def test_export_writes_one_file_per_time(gid_case: Path):
    pvd = export_vtk(gid_case)
    root = ET.parse(pvd).getroot()
    steps = [(float(d.get("timestep")), d.get("file")) for d in root.iter("DataSet")]
    assert [t for t, _ in steps] == [t for t, _ in iter_time_steps(gid_case)]
    assert all((pvd.parent / f).exists() for _, f in steps)


def test_positions_and_values_match_gid_files(gid_case: Path):
    pvd = export_vtk(gid_case)
    mesh = read_gid_mesh(gid_case.with_name("mini.POST.MSH"))
    steps = list(iter_time_steps(gid_case))
    first_disp = dict(zip(steps[0][1]["Displacement"].ids,
                          steps[0][1]["Displacement"].values, strict=True))

    last_file = [d.get("file") for d in ET.parse(pvd).getroot().iter("DataSet")][-1]
    piece = ET.parse(pvd.parent / last_file).getroot().find(".//Piece")
    arrays = {a.get("Name"): _decode(a) for a in piece.iter("DataArray")}
    blocks = steps[-1][1]
    ids = blocks["Displacement"].ids
    np.testing.assert_array_equal(arrays["id"], ids)

    disp = blocks["Displacement"].values
    expected = mesh.coords[ids - 1] + disp - np.array([first_disp[i] for i in ids])
    np.testing.assert_allclose(arrays["Points"][:, :2], expected, rtol=1e-6, atol=1e-6)
    assert not arrays["Points"][:, 2].any()
    np.testing.assert_allclose(arrays["Displacement"][:, :2], disp, rtol=1e-6)
    assert arrays["Displacement"].shape[1] == 3  # vector completado con z = 0
    np.testing.assert_allclose(arrays["Equi_strain"], blocks["Equi_strain"].values[:, 0],
                               rtol=1e-6)
    assert arrays["Total_strain"].shape[1] == 3
    assert set(np.unique(arrays["material"])) == {1}


def test_every_option(gid_case: Path, workdir: Path):
    pvd = export_vtk(gid_case, out_dir=workdir / "vtk", every=2)
    n_steps = len(list(iter_time_steps(gid_case)))
    assert len(list(ET.parse(pvd).getroot().iter("DataSet"))) == (n_steps + 1) // 2


def test_readable_by_vtk(gid_case: Path):
    pv = pytest.importorskip("pyvista")
    reader = pv.get_reader(str(export_vtk(gid_case)))
    reader.set_active_time_point(len(reader.time_values) - 1)
    grid = reader.read()[0]
    assert grid.n_points > 0 and "Equi_strain" in grid.point_data
