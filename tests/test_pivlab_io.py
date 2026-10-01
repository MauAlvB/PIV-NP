from pathlib import Path

import numpy as np
import pytest

from pivnp.config import ConfigError
from pivnp.pivlab_io import FrameSource, pivlab_to_node, read_moisture_file, read_velocity_file

from .legacy_reference import iconectividad

HEADER = "PIVlab\nFRAME: 1\nx [m],y [m],u [m/s],v [m/s]\n"


# --- where the input files live ----------------------------------------------------------
def _write_step(directory: Path, step: int, u: float, moisture: float | None = None) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    (directory / f"datos ({step}).txt").write_text(HEADER + f"0,0,{u},0\n")
    if moisture is not None:
        (directory / f"Moist_{step}.TXT").write_text(f"t\n0,0,{moisture},1\n")


def test_inputs_are_read_from_the_case_root(workdir: Path):
    """How every case was laid out before the subfolders: everything beside the .PAR."""
    _write_step(workdir, 1, 1.5, moisture=0.3)
    frame = FrameSource(workdir, n_nodes=1, moisture=True, prefetch=0).read(1)
    assert frame.u[0] == 1.5 and frame.moisture[0] == 0.3


def test_inputs_are_read_from_their_subfolders(workdir: Path):
    _write_step(workdir / "pivlab", 1, 2.5)
    _write_step(workdir / "moisture", 1, 0.0, moisture=0.7)
    (workdir / "moisture" / "datos (1).txt").unlink()  # only the Moist file belongs there
    source = FrameSource(workdir, n_nodes=1, moisture=True, prefetch=0)
    frame = source.read(1)
    assert frame.u[0] == 2.5 and frame.moisture[0] == 0.7
    assert source.velocity_path(1).parent.name == "pivlab"


def test_the_two_layouts_can_be_mixed(workdir: Path):
    """Velocities in their subfolder and moisture in the root, or the other way round."""
    _write_step(workdir / "pivlab", 1, 3.5)
    (workdir / "Moist_1.TXT").write_text("t\n0,0,0.9,1\n")
    frame = FrameSource(workdir, n_nodes=1, moisture=True, prefetch=0).read(1)
    assert frame.u[0] == 3.5 and frame.moisture[0] == 0.9


def test_the_same_file_in_both_places_is_refused(workdir: Path):
    """Choosing one in silence is how someone edits a file that is not the one being read."""
    _write_step(workdir, 1, 1.0)
    _write_step(workdir / "pivlab", 1, 9.0)
    with pytest.raises(ConfigError, match="pivlab/ subfolder"):
        FrameSource(workdir, n_nodes=1, prefetch=0)


def test_a_missing_step_names_the_file(workdir: Path):
    _write_step(workdir / "pivlab", 1, 1.0)
    source = FrameSource(workdir, n_nodes=1, prefetch=0)
    with pytest.raises(FileNotFoundError, match=r"datos \(2\)"):
        source.read(2)


@pytest.mark.parametrize(("nch", "nfil"), [(1, 1), (3, 2), (59, 34)])
def test_pivlab_to_node_matches_iconectividad(nch, nfil):
    conn = iconectividad(nch, nfil)
    ours = pivlab_to_node(nch, nfil)
    assert [ours[p - 1] + 1 for p in sorted(conn)] == [conn[p] for p in sorted(conn)]
    assert sorted(ours.tolist()) == list(range((nch + 1) * (nfil + 1)))  # a permutation


def test_pivlab_order_is_column_major_from_top():
    # 2x2 cells -> 3x3 points. The first PIVlab point (column 0, top) is node 6.
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


def test_five_columns_are_detected_even_if_the_par_says_four(workdir: Path):
    """A .PAR with IPIVLAB=1 and a five-column file: the file wins.

    Reading five columns as if they were four shifts every value out of place, and that is
    what happens with the old .PAR files, which carry no IPIVLAB and default to 1.
    """
    path = workdir / "d.txt"
    path.write_text(HEADER + "0.1,0.2,1.0,2.0,1\n0.1,0.4,3.0,4.0,0\n")
    x, y, u, v = read_velocity_file(path, 2, pivlab_format=1)
    assert x.tolist() == [0.1, 0.1] and y.tolist() == [0.2, 0.4]
    assert u.tolist() == [1.0, 3.0] and v.tolist() == [2.0, 4.0]


def test_missing_values_raise(workdir: Path):
    path = workdir / "d.txt"
    path.write_text(HEADER + "0.1,0.2,1.0,2.0\n")
    with pytest.raises(ValueError, match="expected 8 values"):
        read_velocity_file(path, 2, pivlab_format=1)


def test_read_moisture(workdir: Path):
    path = workdir / "Moist_1.TXT"
    path.write_text("title\n0,0,0.1,0.5\n0,1,NaN,NaN\n")
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


HEADER_WITH_FACTORS = ("PIVlab\nconversion factor xy (px -> m): 0.004, "
                       "conversion factor uv (px/frame -> m/s): 0.004\nx,y,u,v\n")


def test_mesh_in_metres(workdir: Path):
    (workdir / "datos (1).txt").write_text(HEADER_WITH_FACTORS + "0.1,0.2,1,0\n0.3,0.4,NaN,0\n")
    x, y, factor, has_data = FrameSource(workdir, n_nodes=2).mesh_in_metres()
    assert x.tolist() == [0.1, 0.3] and y.tolist() == [0.2, 0.4]
    assert factor == 0.004 and has_data.tolist() == [True, False]


def test_header_without_conversion_factor(workdir: Path):
    (workdir / "datos (1).txt").write_text(HEADER + "0,0,1,0\n")
    with pytest.raises(ValueError, match="pixels-to-metres"):
        FrameSource(workdir, n_nodes=1).mesh_in_metres()


class _FakeMoistureSource:
    """Stand-in for MoistureSource: records the order in which it is asked for things."""

    def __init__(self) -> None:
        self.grays: list[int] = []
        self.model: list[int] = []

    def normalized_gray(self, step, has_data):
        self.grays.append(step)
        return np.array([10.0 * step])

    def from_gray(self, step, normalized):
        self.model.append(step)

        class State:
            moisture = normalized / 100
            saturation = normalized / 10

        return State()


@pytest.mark.parametrize("prefetch", [0, 3])
def test_moisture_from_images_is_applied_in_order(workdir: Path, prefetch):
    """The model carries memory, so it has to see the steps in order even though the images
    were read ahead of time on several threads."""
    for step in (1, 2, 3):
        (workdir / f"datos ({step}).txt").write_text(HEADER + f"0,0,{step},0\n")
    moisture = _FakeMoistureSource()
    source = FrameSource(workdir, n_nodes=1, moisture=True, prefetch=prefetch, images=moisture)
    frames = list(source.frames(range(1, 4)))

    assert moisture.model == [1, 2, 3]
    assert sorted(moisture.grays) == [1, 2, 3]
    assert [f.saturation[0] for f in frames] == [1.0, 2.0, 3.0]
    assert [f.moisture[0] for f in frames] == [0.1, 0.2, 0.3]
    assert not source.moisture  # the Moist_<n>.TXT are not read, and do not exist here
