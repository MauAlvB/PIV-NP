"""Conversion of GiD results (``.POST.MSH`` + ``.POST.RES``) to VTK for ParaView.

It writes one ``.vtu`` file per step and a ``.pvd`` collection that ParaView opens as an
animation. It works both with the results of this version and with those of the original
Fortran executable.

Usage::

    pivnp-vtk path/case.POST.RES            # looks for case.POST.MSH next to it
    pivnp-vtk path/case.POST.RES -o vtk --every 5

The coordinates of every step are the **current** positions of the particles:
``mesh + Displacement``. With ``--legacy-msh`` the displacement of the first step is also
subtracted, which is needed for results of the original Fortran, where the mesh was written
after the first step.
"""

from __future__ import annotations

import argparse
import base64
import sys
import zlib
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path

import numpy as np

#: Results that are physical vectors (padded with z = 0 so that arrows can be used).
VECTORS = {"Displacement", "Inst_displaceme", "Inst_displacement", "Velocity", "Acceleration"}
COMPONENT_NAMES = {
    "Total_strain": ("xx", "yy", "xy"),
    "Inc_strain": ("xx", "yy", "xy"),
    "E_kinetic": ("x", "y"),
}
RENAMES = {"Inst_displaceme": "Inst_displacement"}  # nombre recortado por el formato A15


@dataclass(frozen=True)
class GidMesh:
    ids: np.ndarray  # (n,)
    coords: np.ndarray  # (n, 2)
    material: np.ndarray  # (n,)


@dataclass(frozen=True)
class ResultBlock:
    name: str
    time: float
    ids: np.ndarray
    values: np.ndarray  # (n, k)


# --- reading GiD -----------------------------------------------------------------------------
def read_gid_mesh(path: Path) -> GidMesh:
    lines = Path(path).read_text(encoding="latin-1").splitlines()
    start, end = lines.index("Coordinates") + 2, lines.index("End Coordinates")
    coords = np.loadtxt(lines[start:end], ndmin=2)
    start, end = lines.index("Elements") + 2, lines.index("End Elements")
    elements = np.loadtxt(lines[start:end], dtype=np.int64, ndmin=2)
    material = np.zeros(coords.shape[0], dtype=np.int64)
    material[elements[:, 0] - 1] = elements[:, 2]
    return GidMesh(coords[:, 0].astype(np.int64), coords[:, 1:3], material)


def iter_gid_results(path: Path) -> Iterator[ResultBlock]:
    """Walk the ``Result`` blocks of a ``.POST.RES`` without loading it all into memory."""
    with open(path, encoding="latin-1") as f:
        header = None
        for line in f:
            if line.startswith("Result"):
                tokens = line.split()
                header = (RENAMES.get(tokens[1], tokens[1]), float(tokens[3]))
            elif line.startswith("Values") and header:
                rows = []
                for row in f:
                    if row.startswith("End values"):
                        break
                    rows.append(row)
                data = np.loadtxt(rows, ndmin=2) if rows else np.zeros((0, 2))
                yield ResultBlock(header[0], header[1], data[:, 0].astype(np.int64), data[:, 1:])
                header = None


def iter_time_steps(path: Path) -> Iterator[tuple[float, dict[str, ResultBlock]]]:
    """Group the blocks by step."""
    current: dict[str, ResultBlock] = {}
    time = None
    for block in iter_gid_results(path):
        if time is not None and block.time != time:
            yield time, current
            current = {}
        time = block.time
        current[block.name] = block
    if current:
        yield time, current


# --- writing VTK -----------------------------------------------------------------------------
def _encode(array: np.ndarray) -> str:
    """Datos binarios comprimidos con zlib en base64 (formato 'binary' de VTK XML)."""
    raw = np.ascontiguousarray(array).tobytes()
    compressed = zlib.compress(raw, 6)
    header = np.array([1, len(raw), len(raw), len(compressed)], dtype="<u4").tobytes()
    return base64.b64encode(header).decode() + base64.b64encode(compressed).decode()


def _data_array(name: str, array: np.ndarray, components: tuple[str, ...] = ()) -> str:
    vtk_type = {"f": "Float32", "i": "Int32", "u": "UInt8"}[array.dtype.kind]
    n_comp = 1 if array.ndim == 1 else array.shape[1]
    names = "".join(f' ComponentName{i}="{c}"' for i, c in enumerate(components))
    return (f'<DataArray type="{vtk_type}" Name="{name}" NumberOfComponents="{n_comp}"'
            f'{names} format="binary">{_encode(array)}</DataArray>')


def write_vtu(path: Path, points: np.ndarray, point_data: dict[str, np.ndarray]) -> None:
    """Point cloud (VTK_VERTEX cells) with per-point data."""
    n = points.shape[0]
    xyz = np.zeros((n, 3), dtype=np.float32)
    xyz[:, : points.shape[1]] = points
    arrays = []
    for name, values in point_data.items():
        values = np.asarray(values)
        if values.dtype.kind == "f":
            values = values.astype(np.float32)
        elif values.dtype.kind in "iu":
            values = values.astype(np.int32)
        if values.ndim == 2 and values.shape[1] == 1:
            values = values[:, 0]
        if name in VECTORS and values.ndim == 2 and values.shape[1] == 2:
            values = np.column_stack([values, np.zeros(n, np.float32)])
        arrays.append(_data_array(name, values, COMPONENT_NAMES.get(name, ())))
    idx = np.arange(n, dtype=np.int32)
    Path(path).write_text(
        '<?xml version="1.0"?>\n'
        '<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian" '
        'header_type="UInt32" compressor="vtkZLibDataCompressor">\n'
        f'<UnstructuredGrid><Piece NumberOfPoints="{n}" NumberOfCells="{n}">\n'
        f'<PointData>{"".join(arrays)}</PointData>\n'
        f'<Points>{_data_array("Points", xyz)}</Points>\n'
        f'<Cells>{_data_array("connectivity", idx)}{_data_array("offsets", idx + 1)}'
        f'{_data_array("types", np.ones(n, dtype=np.uint8))}</Cells>\n'
        '</Piece></UnstructuredGrid>\n</VTKFile>\n',
        encoding="ascii",
    )


def write_pvd(path: Path, entries: list[tuple[float, str]]) -> None:
    datasets = "\n".join(f'  <DataSet timestep="{t!r}" part="0" file="{f}"/>'
                         for t, f in entries)
    Path(path).write_text(
        '<?xml version="1.0"?>\n<VTKFile type="Collection" version="1.0">\n<Collection>\n'
        f"{datasets}\n</Collection>\n</VTKFile>\n", encoding="ascii")


# --- conversion ------------------------------------------------------------------------------
def _align(block: ResultBlock, ids: np.ndarray, n_total: int) -> np.ndarray:
    """Values of the block in the order of ``ids`` (NaN when a particle does not appear)."""
    if np.array_equal(block.ids, ids):
        return block.values
    dense = np.full((n_total + 1, block.values.shape[1]), np.nan)
    dense[block.ids] = block.values
    return dense[ids]


def export_vtk(res_path: Path, msh_path: Path | None = None, out_dir: Path | None = None,
               every: int = 1, legacy_mesh: bool = False) -> Path:
    """Convert a GiD case to VTK. Returns the path of the ``.pvd``.

    ``legacy_mesh``: the mesh holds the positions of the first step instead of the initial
    ones (results of the original Fortran, or computed with ``--legacy-compat``).
    """
    res_path = Path(res_path)
    case = res_path.name.removesuffix(".POST.RES").removesuffix(".post.res")
    msh_path = Path(msh_path) if msh_path else res_path.with_name(f"{case}.POST.MSH")
    out_dir = Path(out_dir) if out_dir else res_path.parent / f"{case}_vtk"
    out_dir.mkdir(parents=True, exist_ok=True)

    mesh = read_gid_mesh(msh_path)
    n_total = int(mesh.ids.max())
    coords = np.zeros((n_total + 1, 2))
    coords[mesh.ids] = mesh.coords
    material = np.zeros(n_total + 1, dtype=np.int64)
    material[mesh.ids] = mesh.material
    first_disp = None

    entries = []
    for k, (time, blocks) in enumerate(iter_time_steps(res_path)):
        disp_block = blocks["Displacement"]
        if first_disp is None:
            first_disp = np.zeros((n_total + 1, 2))
            if legacy_mesh:
                first_disp[disp_block.ids] = disp_block.values
        if k % every:
            continue
        ids = disp_block.ids
        positions = coords[ids] + disp_block.values - first_disp[ids]
        data = {"id": ids, "material": material[ids]}
        data.update({name: _align(b, ids, n_total) for name, b in blocks.items()})
        file_name = f"{case}_{k:04d}.vtu"
        write_vtu(out_dir / file_name, positions, data)
        entries.append((time, file_name))

    pvd = out_dir / f"{case}.pvd"
    write_pvd(pvd, entries)
    return pvd


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="pivnp-vtk", description="Convert PIV-NP GiD results to VTK (ParaView).")
    parser.add_argument("res", type=Path, help="the <case>.POST.RES file")
    parser.add_argument("--msh", type=Path, help=".POST.MSH file (next to the .RES by default)")
    parser.add_argument("-o", "--out", type=Path, help="output folder (<case>_vtk by default)")
    parser.add_argument("--every", type=int, default=1, help="export 1 out of every N steps")
    parser.add_argument("--legacy-msh", action="store_true",
                        help="the mesh holds the positions of the first step (results of the "
                             "original Fortran, or computed with --legacy-compat)")
    args = parser.parse_args(argv)
    pvd = export_vtk(args.res, args.msh, args.out, max(1, args.every), args.legacy_msh)
    print(f"Open in ParaView: {pvd}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
