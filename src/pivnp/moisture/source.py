"""Orchestration: from the test images to the water content and saturation of every node.

It puts the previous pieces together and delivers, step by step, the two fields the analysis
needs. It can also write the ``Moist_<n>.TXT`` files, the format the MATLAB code worked with,
which is handy for comparing results.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .calibration import Calibration
from .images import gaussian_blur, read_image
from .model import (
    DRY_LIMIT,
    MEASURED,
    NO_DATA,
    WET_LIMIT,
    MoistureModel,
    References,
    State,
    global_references,
    normalize,
    references_from_offsets,
    warn_if_band_is_narrow,
)
from .sampling import outside_the_image, pixel_coordinates, sample
from .settings import MoistureSettings, read

MOIST_HEADER = "x_m,y_m,moisture,saturation_degree"

log = logging.getLogger("pivnp")


@dataclass(frozen=True)
class Mesh:
    """Nodes of the PIV grid, in metres, and the conversion factor to pixels."""

    x_m: np.ndarray
    y_m: np.ndarray
    metres_per_pixel: float

    def __post_init__(self) -> None:
        if self.x_m.shape != self.y_m.shape:
            raise ValueError("x and y must have the same number of nodes")
        if self.metres_per_pixel <= 0:
            raise ValueError("the conversion factor must be positive")


class MoistureSource:
    """Compute water content and saturation from the images, step by step."""

    def __init__(self, directory: Path, settings: MoistureSettings,
                 calibration: Calibration, mesh: Mesh) -> None:
        self.directory = Path(directory)
        self.settings = settings
        self.mesh = mesh
        self.model = MoistureModel(calibration, settings.saturation_threshold,
                                   settings.incremental)
        self.column, self.row = pixel_coordinates(mesh.x_m, mesh.y_m, mesh.metres_per_pixel,
                                                  settings.registration)
        self.references: References | None = None
        self.nodes_outside = 0
        #: How many values of each quality mark have been published in the whole analysis.
        self.quality: dict[int, int] = dict.fromkeys((MEASURED, NO_DATA, WET_LIMIT,
                                                      DRY_LIMIT), 0)

    # --- preparation ----------------------------------------------------------------------
    def prepare(self, has_data: np.ndarray) -> References:
        """Fix the dry and saturated references, which do not change during the test.

        With ``DRY_BAND`` and ``SATURATED_BAND`` in the ``.HUM`` they are two intensities for
        the whole image and no reference image is needed. Otherwise they are taken from the
        reference image node by node, which is what the RGB flow does.
        """
        settings = self.settings
        if settings.has_global_band:
            self.references = global_references(self.column.shape, settings.dry_band,
                                                settings.saturated_band)
        else:
            filtered = self._filtered_image(self.directory / settings.dry_reference)
            self.nodes_outside = outside_the_image(filtered, self.column, self.row)
            gray = sample(filtered, self.column, self.row, has_data)
            self.references = references_from_offsets(gray, settings.dry_offset,
                                                      settings.saturated_offset)
        warn_if_band_is_narrow(self.references, settings.source)
        return self.references

    # --- computation ----------------------------------------------------------------------
    def normalized_gray(self, step: int, has_data: np.ndarray) -> np.ndarray:
        """Gray of the step, on the 0 (saturated) - 100 (dry) scale."""
        if self.references is None:
            raise RuntimeError("prepare() must be called before the first step")
        filtered = self._filtered_image(self.settings.image_path(self.directory, step))
        gray = sample(filtered, self.column, self.row, has_data)
        return normalize(gray, self.references)

    def _filtered_image(self, path: Path) -> np.ndarray:
        image = read_image(path, self.settings.channel)
        return gaussian_blur(image, self.settings.sigma, self.settings.legacy_rounding)

    def quality_summary(self) -> str:
        """How much of what is published is a measurement and how much a bound."""
        total = sum(self.quality.values())
        if not total:
            return "no step has been computed"
        with_data = total - self.quality[NO_DATA]
        if not with_data:
            return "no node with data"
        parts = [f"{self.quality[mark]} ({100 * self.quality[mark] / with_data:.0f} %) {name}"
                 for mark, name in ((MEASURED, "measured"),
                                    (WET_LIMIT, "at the wet limit (an 'at least')"),
                                    (DRY_LIMIT, "at the dry limit (an 'at most')"))
                 if self.quality[mark]]
        return f"of {with_data} values with data: " + ", ".join(parts)

    def from_gray(self, step: int, normalized: np.ndarray) -> State:
        """Water content and saturation from the already normalized gray.

        Kept apart from reading the image because this part carries memory: the saturation
        does not go down, so the steps have to come through here in order. What comes before
        carries none and can be computed ahead of time on another thread.
        """
        state = self.model.evaluate(normalized)
        for mark in self.quality:
            self.quality[mark] += int((state.quality == mark).sum())
        if step == 1 and self.settings.first_step == "legacy":
            # MATLAB left the water content at zero on the first step.
            state = State(state.saturation, np.zeros_like(state.moisture), state.quality)
        return state

    def at_step(self, step: int, has_data: np.ndarray) -> State:
        """Water content and saturation of one step."""
        return self.from_gray(step, self.normalized_gray(step, has_data))

    def reset(self) -> None:
        self.model.reset()


def source_for_case(case_dir: Path, case_name: str, mesh: Mesh,
                    has_data: np.ndarray) -> MoistureSource:
    """Prepare the moisture source of a case from its ``<case>.HUM``.

    The calibration file is looked up next to the case file, unless an absolute path is given.
    """
    from ..config import find_file  # imported here to avoid a circular dependency

    case_dir = Path(case_dir)
    settings = read(find_file(case_dir, f"{case_name}.HUM"))
    if settings.unknown_keys:
        log.warning("%s: keys that are not recognized and are ignored: %s", settings.source,
                    ", ".join(settings.unknown_keys))
    calibration_path = Path(settings.calibration)
    if not calibration_path.is_absolute():
        calibration_path = find_file(case_dir, settings.calibration)
    source = MoistureSource(case_dir, settings, Calibration.from_csv(calibration_path), mesh)
    source.prepare(has_data)
    if source.nodes_outside:
        log.warning("%d grid nodes fall outside the moisture image: check the registration "
                    "(SCALE_X, OFFSET_X, SCALE_Y, OFFSET_Y) in the .HUM", source.nodes_outside)
    log.info("Moisture from the images: %s, channel %d, sigma %g, calibration %s",
             settings.image_pattern, settings.channel, settings.sigma,
             Path(calibration_path).name)
    return source


def write_moist(path: Path, mesh: Mesh, moisture: np.ndarray,
                saturation: np.ndarray) -> None:
    """Write a ``Moist_<n>.TXT`` file in the format of the MATLAB code."""
    lines = [MOIST_HEADER]
    for x, y, m, s in zip(mesh.x_m, mesh.y_m, moisture, saturation, strict=True):
        lines.append(f"{_number(x, 10)},{_number(y, 10)},"
                     f"{_number(m, 17)},{_number(s, 15)}")
    Path(path).write_text("\n".join(lines) + "\n", encoding="ascii")


def _number(value: float, digits: int) -> str:
    """MATLAB's format: ``NaN`` for what is missing and short notation for the rest."""
    if not np.isfinite(value):
        return "NaN"
    return f"{float(value):.{digits}g}"
