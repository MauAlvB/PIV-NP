from pathlib import Path

import numpy as np
import pytest

from pivnp.pivlab_io import FrameSource, pivlab_to_node, read_moisture_file, read_velocity_file

from .legacy_reference import iconectividad

HEADER = "PIVlab\nFRAME: 1\nx [m],y [m],u [m/s],v [m/s]\n"


@pytest.mark.parametrize(("nch", "nfil"), [(1, 1), (3, 2), (59, 34)])
def test_pivlab_to_node_matches_iconectividad(nch, nfil):
    conn = iconectividad(nch, nfil)
    ours = pivlab_to_node(nch, nfil)
    assert [ours[p - 1] + 1 for p in sorted(conn)] == [conn[p] for p in sorted(conn)]
    assert sorted(ours.tolist()) == list(range((nch + 1) * (nfil + 1)))  # permutación


def test_pivlab_order_is_column_major_from_top():
    # 2x2 celdas -> 3x3 puntos. El primer punto PIVlab (columna 0, arriba) es el nodo 6.
    assert pivlab_to_node(2, 2).tolist() == [6, 3, 0, 7, 4, 1, 8, 5, 2]


def test_read_four_column_file_with_nan(workdir: Path):
    path = workdir / "datos (1).txt"
    path.write_text(HEADER + "0.1,0.2,NaN,NaN\n0.1,0.4,1.5,-2e-3\n")
    x, y, u, v = read_velocity_file(path, 2, pivlab_format=1)
    assert np.isnan(u[0]) and np.isnan(v[0])
    assert (u[1], v[1]) == (1.5, -2e-3)


def test_four_column_read_is_free_format(workdir: Path):
    path = workdir / "d.txt"
    path.write_text(HEADER + "0.1 0.2 1\n2 0.1\n0.4 3 4\nsobra\n")
    _, _, u, v = read_velocity_file(path, 2, pivlab_format=1)
    assert u.tolist() == [1.0, 3.0] and v.tolist() == [2.0, 4.0]


def test_read_five_column_file(workdir: Path):
    path = workdir / "d.txt"
    path.write_text(HEADER + "0.1,0.2,1.0,2.0,1\n0.1,0.4,3.0,4.0,0,extra\n")
    _, _, u, v = read_velocity_file(path, 2, pivlab_format=2)
    assert u.tolist() == [1.0, 3.0] and v.tolist() == [2.0, 4.0]


def test_missing_values_raise(workdir: Path):
    path = workdir / "d.txt"
    path.write_text(HEADER + "0.1,0.2,1.0,2.0\n")
    with pytest.raises(ValueError, match="se esperaban 8"):
        read_velocity_file(path, 2, pivlab_format=1)


def test_read_moisture(workdir: Path):
    path = workdir / "Moist_1.TXT"
    path.write_text("titulo\n0,0,0.1,0.5\n0,1,NaN,NaN\n")
    moisture, saturation = read_moisture_file(path, 2)
    assert moisture[0] == 0.1 and saturation[0] == 0.5 and np.isnan(moisture[1])


def test_frame_source_order_prefetch_and_case(workdir: Path):
    for step in (1, 2, 3):
        (workdir / f"DATOS ({step}).TXT").write_text(HEADER + f"0,0,{step},0\n")
        (workdir / f"moist_{step}.txt").write_text(f"t\n0,0,{step / 10},1\n")
    source = FrameSource(workdir, n_nodes=1, moisture=True, prefetch=2)
    frames = list(source.frames(range(1, 4)))
    assert [f.u[0] for f in frames] == [1.0, 2.0, 3.0]
    assert [f.moisture[0] for f in frames] == [0.1, 0.2, 0.3]
    with pytest.raises(FileNotFoundError):
        source.read(4)


def test_frame_source_without_prefetch_or_moisture(workdir: Path):
    (workdir / "datos (1).txt").write_text(HEADER + "0,0,1,0\n")
    frame = next(FrameSource(workdir, n_nodes=1, prefetch=0).frames(range(1, 2)))
    assert frame.moisture.tolist() == [0.0] and frame.saturation.tolist() == [0.0]
