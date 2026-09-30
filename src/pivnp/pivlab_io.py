"""Reading the PIVlab files (``datos (n).txt``) and the moisture ones (``Moist_n.TXT``).

PIVlab stores the points by columns: constant x and y growing *downwards* (the image axis).
PIV-NP numbers the nodes by rows from the bottom left and uses the y axis pointing up, which
is why the nodes are reordered and the sign of the vertical velocity is flipped.
"""

from __future__ import annotations

import logging
import re
from collections.abc import Iterator
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, replace
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:  # types only: the moisture package does not depend on this module
    from .moisture.source import MoistureSource

VELOCITY_PATTERN = "datos ({step}).TXT"
MOISTURE_PATTERN = "Moist_{step}.TXT"

log = logging.getLogger("pivnp")

#: Separators of the PIVlab files: commas, spaces or tabs.
_TOKENS = re.compile(r"[,\s]+")

#: Conversion factors PIVlab writes on the second line of every file.
_XY_FACTOR = re.compile(r"px\s*->\s*m\)\s*:\s*([0-9.eE+-]+)")
_UV_FACTOR = re.compile(r"px/frame\s*->\s*m/s\)\s*:\s*([0-9.eE+-]+)")


@dataclass(frozen=True)
class Frame:
    """Data of one step, in the order of the points of the PIVlab file."""

    step: int
    source: Path
    u: np.ndarray  # x velocity (NaN when PIVlab did not measure it)
    v: np.ndarray  # y velocity, image axis (downwards)
    moisture: np.ndarray
    saturation: np.ndarray
    #: Gray of the test image, already normalized, when the moisture is computed from the
    #: images. It is the intermediate step: the moisture comes from applying the model to it,
    #: and the model carries memory of the previous steps, so it cannot be computed here.
    normalized_gray: np.ndarray | None = None


def pivlab_to_node(n_cols: int, n_rows: int) -> np.ndarray:
    """PIV-NP node of every PIVlab point (``ICONECTIVIDAD`` in the original).

    Point ``p = col * (n_rows + 1) + row_from_top`` corresponds to node
    ``(n_rows - row_from_top) * (n_cols + 1) + col``.
    """
    col, row_from_top = np.divmod(np.arange((n_cols + 1) * (n_rows + 1)), n_rows + 1)
    return (n_rows - row_from_top) * (n_cols + 1) + col


def xy_factor_in_header(path: Path) -> float:
    """Metres per pixel PIVlab exported the file with (second line)."""
    lines = Path(path).read_text(encoding="latin-1").splitlines()
    found = _XY_FACTOR.search(lines[1]) if len(lines) >= 2 else None
    if not found:
        raise ValueError(f"{path}: the header does not carry the pixels-to-metres factor, "
                         "which is needed to place the nodes on the image")
    return float(found.group(1))


def frame_interval_in_header(path: Path) -> float | None:
    """Interval between images the file was exported with, or ``None`` if it is not stated.

    PIVlab writes two conversion factors: one from pixels to metres and another from pixels
    per frame to metres per second. Their ratio is the time between images it was told about,
    which should match the DT of the ``.PAR``: if it does not, displacements come out scaled.
    """
    lines = Path(path).read_text(encoding="latin-1").splitlines()
    if len(lines) < 2:
        return None
    xy, uv = _XY_FACTOR.search(lines[1]), _UV_FACTOR.search(lines[1])
    if not xy or not uv:
        return None
    try:
        divisor = float(uv.group(1))
        return float(xy.group(1)) / divisor if divisor else None
    except (ValueError, ZeroDivisionError):
        return None


def _parse_values(lines: list[str], count: int, path: Path) -> np.ndarray:
    tokens = " ".join(lines).replace(",", " ").split()
    if len(tokens) < count:
        raise ValueError(f"{path}: expected {count} values, found {len(tokens)}")
    return np.array(tokens[:count], dtype=np.float64)


def read_velocity_file(path: Path, n_nodes: int = 0,
                       pivlab_format: int = 0) -> tuple[np.ndarray, ...]:
    """Read ``x, y, u, v`` from a PIVlab file (3 header lines).

    The format is deduced from the file itself and not from the ``IPIVLAB`` of the ``.PAR``:
    if every node takes one line with four values or more, the first four of each are taken;
    otherwise 4·NN values are read in a row, ignoring line breaks. Both paths give the same
    thing with four-column files, and it is the only way not to get five-column ones wrong:
    reading them as if they had four shifts every value from the first one on.

    ``pivlab_format`` is only used to warn when the ``.PAR`` says otherwise.
    """
    lines = [line for line in Path(path).read_text(encoding="latin-1").splitlines()[3:]
             if line.strip()]
    by_lines = (len(lines) >= n_nodes > 0
                and all(len(_TOKENS.split(line.strip())) >= 4 for line in lines[:n_nodes]))
    if by_lines:
        if pivlab_format == 1:
            log.warning("%s: the .PAR says IPIVLAB=1 (a single list of values), but the file "
                        "carries one node per line; it is read by lines", path.name)
        data = np.array([_TOKENS.split(line.strip())[:4] for line in lines[:n_nodes]],
                        dtype=np.float64)
    else:
        data = _parse_values(lines, 4 * n_nodes, path).reshape(n_nodes, 4)
    return data[:, 0], data[:, 1], data[:, 2], data[:, 3]


def read_moisture_file(path: Path, n_nodes: int) -> tuple[np.ndarray, np.ndarray]:
    """Read moisture and saturation from ``Moist_n.TXT`` (1 header line, 4 columns)."""
    lines = Path(path).read_text(encoding="latin-1").splitlines()[1:]
    data = _parse_values(lines, 4 * n_nodes, path).reshape(n_nodes, 4)
    return data[:, 2].copy(), data[:, 3].copy()


class FrameSource:
    """Provider of PIVlab steps, reading ahead on threads."""

    def __init__(
        self,
        directory: Path,
        n_nodes: int,
        pivlab_format: int = 1,
        moisture: bool = False,
        prefetch: int = 4,
        images: MoistureSource | None = None,
    ) -> None:
        self.directory = Path(directory)
        self.n_nodes = n_nodes
        self.pivlab_format = pivlab_format
        #: Read the moisture from the ``Moist_<n>.TXT``; incompatible with computing it.
        self.moisture = moisture and images is None
        self.prefetch = max(0, prefetch)
        self.images = images
        self._index = {p.name.lower(): p for p in self.directory.iterdir()}

    def _path(self, pattern: str, step: int) -> Path:
        name = pattern.format(step=step)
        try:
            return self._index[name.lower()]
        except KeyError:
            raise FileNotFoundError(self.directory / name) from None

    def velocity_path(self, step: int) -> Path:
        """PIVlab file of step ``step``."""
        return self._path(VELOCITY_PATTERN, step)

    def mesh_in_metres(self, step: int = 1) -> tuple[np.ndarray, np.ndarray, float, np.ndarray]:
        """Nodes in metres, metres per pixel and which nodes PIVlab measured, from a file."""
        path = self._path(VELOCITY_PATTERN, step)
        x, y, u, _ = read_velocity_file(path, self.n_nodes, self.pivlab_format)
        return x, y, xy_factor_in_header(path), np.isfinite(u)

    def read(self, step: int) -> Frame:
        """Read one step. It does not apply the moisture model: that goes in order, in
        ``frames``."""
        path = self._path(VELOCITY_PATTERN, step)
        _, _, u, v = read_velocity_file(path, self.n_nodes, self.pivlab_format)
        if self.moisture:
            moisture, saturation = read_moisture_file(
                self._path(MOISTURE_PATTERN, step), self.n_nodes
            )
        else:
            moisture = saturation = np.zeros(self.n_nodes)
        gray = self.images.normalized_gray(step, np.isfinite(u)) if self.images else None
        return Frame(step, path, u, v, moisture, saturation, gray)

    def _with_moisture(self, frame: Frame) -> Frame:
        """Apply the moisture model, which needs the steps in order."""
        if self.images is None or frame.normalized_gray is None:
            return frame
        state = self.images.from_gray(frame.step, frame.normalized_gray)
        return replace(frame, moisture=state.moisture, saturation=state.saturation)

    def frames(self, steps: range) -> Iterator[Frame]:
        """Yield the steps in order, reading the following ones in the background."""
        if self.prefetch == 0:
            yield from (self._with_moisture(self.read(step)) for step in steps)
            return
        with ThreadPoolExecutor(max_workers=self.prefetch) as pool:
            pending = [pool.submit(self.read, s) for s in steps[: self.prefetch]]
            for k in range(len(steps)):
                frame = pending[k].result()
                ahead = k + self.prefetch
                if ahead < len(steps):
                    pending.append(pool.submit(self.read, steps[ahead]))
                pending[k] = None  # libera memoria
                yield self._with_moisture(frame)
