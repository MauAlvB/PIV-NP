"""The displacement source that measures its own velocities from the photographs.

Selected with ``pivnp <case> --source images``. Everything downstream is unchanged: the
analysis asks for a :class:`~pivnp.sources.Frame` per step and never learns where it came
from, which is what the source contract was built for.

What this has to get right, beyond the correlation itself, is the bookkeeping that PIVlab
normally does: the grid of interrogation windows has to be the grid the ``.PAR`` describes,
the displacement in pixels has to become a velocity in metres per second, and the vertical
axis has to end up the way the rest of PIV-NP expects it.
"""

from __future__ import annotations

import logging
from collections.abc import Iterator
from dataclasses import replace
from pathlib import Path

import numpy as np

from ..config import CaseConfig, ConfigError, find_file
from ..moisture.images import read_image
from ..sources import Frame, register_source
from .correlation import analyse, window_centres
from .settings import PivSettings, read

log = logging.getLogger("pivnp")


class ImageSource:
    """Measures the velocity field from the photographs of the test, step by step."""

    def __init__(self, directory: Path, settings: PivSettings, config: CaseConfig) -> None:
        self.directory = Path(directory)
        self.settings = settings
        self.config = config
        self.moisture = False
        self.images = None
        self._shape: tuple[int, int] | None = None
        self._mask: np.ndarray | None = None
        #: The photograph read last, kept because step *n* reads images *n* and *n+1* and
        #: step *n+1* then reads *n+1* again. On a 3840x2160 JPEG the decoding is most of
        #: the time the analysis takes, so remembering one image halves the run. Only one:
        #: these are tens of megabytes each once decoded.
        self._last: tuple[int, np.ndarray] | None = None

    # --- the photographs -------------------------------------------------------------
    def _read(self, step: int) -> np.ndarray:
        if self._last is not None and self._last[0] == step:
            return self._last[1]
        path = self.settings.image_path(self.directory, step)
        if not path.exists():
            raise FileNotFoundError(
                f"{path}: the photograph for step {step} is not there. The pattern is "
                f"{self.settings.image_pattern!r} and the first image is numbered "
                f"{self.settings.first_image}."
            )
        image = read_image(path, self.settings.channel).astype(np.float64)
        if self._shape is None:
            self._shape = image.shape
            self._mask = self._build_mask(image)
        elif image.shape != self._shape:
            raise ConfigError(f"{path.name} is {image.shape[1]}x{image.shape[0]} px and the "
                              f"first photograph is {self._shape[1]}x{self._shape[0]}; every "
                              "image of a sequence has to be the same size")
        self._last = (step, image)
        return image

    def _build_mask(self, first: np.ndarray) -> np.ndarray | None:
        """Where there is material to measure, if the settings say how to tell."""
        if self.settings.mask_image:
            path = find_file(self.directory, self.settings.mask_image)
            drawn = read_image(path, 0).astype(np.float64)
            if drawn.shape != first.shape:
                raise ConfigError(f"{path.name} is not the size of the photographs")
            return drawn > 0
        if self.settings.mask_below is not None:
            return first > self.settings.mask_below
        return None

    # --- the grid --------------------------------------------------------------------
    def grid_shape(self) -> tuple[int, int]:
        """Rows and columns of interrogation windows, which is the grid of PIV points."""
        if self._shape is None:
            self._read(1)
        rows, cols, _ = window_centres(self._shape, self.settings.window,
                                       self.settings.overlap, self.settings.region)
        return rows.size, cols.size

    def describe_grid(self) -> str:
        """The block 2 the ``.PAR`` would need for these photographs and settings."""
        rows, cols = self.grid_shape()
        cell = self.settings.window * (1.0 - self.settings.overlap) * self.settings.scale
        return (f"n_cells {(cols - 1) * (rows - 1)}  n_nodes {cols * rows}  "
                f"n_rows {rows - 1}  width {cell:.6g}  height {cell:.6g}")

    def check_against_the_par(self) -> None:
        """Refuse to run on a grid the ``.PAR`` does not describe, and say what to put."""
        rows, cols = self.grid_shape()
        if cols * rows != self.config.n_nodes or rows - 1 != self.config.n_rows:
            raise ConfigError(
                f"the photographs and the .PIV settings give a grid of {cols} x {rows} "
                f"points, and the .PAR describes {self.config.n_cols + 1} x "
                f"{self.config.n_rows + 1}. Block 2 of the .PAR should read:\n"
                f"    {self.describe_grid()}"
            )

    # --- the contract ----------------------------------------------------------------
    def mesh_in_metres(self, step: int = 1) -> tuple[np.ndarray, np.ndarray, float,
                                                     np.ndarray]:
        """Where the grid points sit, in metres, in the order PIVlab writes them."""
        self._read(step)
        rows, cols, _ = window_centres(self._shape, self.settings.window,
                                       self.settings.overlap, self.settings.region)
        # by columns, y downwards: the order the rest of the code expects
        x = np.repeat(cols, rows.size) * self.settings.scale
        y = np.tile(rows, cols.size) * self.settings.scale
        field = self._field(step if step < self.config.total_steps else 1)
        return x, y, self.settings.scale, self._flatten(field.measured).astype(bool)

    def frame_interval(self) -> float | None:
        """Unknown: the photographs do not say how far apart in time they were taken.

        The ``DT`` of the ``.PAR`` is the only statement of it, so there is nothing to check
        it against and the analysis skips the comparison.
        """
        return None

    def _flatten(self, grid: np.ndarray) -> np.ndarray:
        """From the (rows, cols) grid to the PIVlab point order: by columns, y downwards."""
        return grid.T.ravel()

    def _field(self, step: int):
        first = self._read(step)
        second = self._read(step + 1)
        return analyse(first, second, window=self.settings.window,
                       overlap=self.settings.overlap, passes=self.settings.passes,
                       mask=self._mask, threshold=self.settings.outlier_threshold,
                       region=self.settings.region, smoothing=self.settings.smoothing,
                       between_pixels=self.settings.between_pixels)

    def read(self, step: int) -> Frame:
        """Measure one step and hand it over as the analysis expects it."""
        field = self._field(step)
        seconds = self.config.dt
        # pixels per step to metres per second; v keeps the image axis, pointing down,
        # because that is what load_measurements expects to flip
        u = self._flatten(field.u) * self.settings.scale / seconds
        v = self._flatten(field.v) * self.settings.scale / seconds
        zeros = np.zeros(u.size)
        return Frame(step, f"{self.settings.image_path(self.directory, step).name} -> "
                           f"{self.settings.image_path(self.directory, step + 1).name}",
                     u, v, zeros, zeros)

    def frames(self, steps: range) -> Iterator[Frame]:
        """Yield the steps in order. The correlation is the slow part, not the reading."""
        for step in steps:
            frame = self.read(step)
            if self.images is not None:
                gray = self.images.normalized_gray(step, np.isfinite(frame.u))
                state = self.images.from_gray(step, gray)
                frame = replace(frame, moisture=state.moisture, saturation=state.saturation,
                                normalized_gray=gray)
            yield frame


@register_source("images")
def build_image_source(case_dir: Path, config: CaseConfig, prefetch: int = 4) -> ImageSource:
    """Build the built-in PIV source of a case, as ``--source images`` does."""
    case_dir = Path(case_dir)
    name = config.case_name if hasattr(config, "case_name") else None
    candidates = [case_dir / f"{name}.PIV"] if name else []
    candidates += sorted(case_dir.glob("*.PIV")) + sorted(case_dir.glob("*.piv"))
    if not candidates:
        raise ConfigError(
            f"{case_dir}: --source images needs a .PIV file saying where the photographs are "
            "and how many metres a pixel is worth."
        )
    settings = read(candidates[0])
    if settings.unknown_keys:
        log.warning("%s: keys that are not recognized and are ignored: %s", settings.source,
                    ", ".join(settings.unknown_keys))
    source = ImageSource(case_dir, settings, config)
    source.check_against_the_par()
    log.info("PIV from the photographs: %s, window %d px, overlap %g, %d passes, "
             "%g m/px", settings.image_pattern, settings.window, settings.overlap,
             settings.passes, settings.scale)
    return source
