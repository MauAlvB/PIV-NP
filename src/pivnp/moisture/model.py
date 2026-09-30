"""From the measured gray value to degree of saturation and water content.

Three steps, the same ones the original MATLAB code had:

1. **References**: from the gray value of a reference image, the gray of the dry soil and of
   the saturated soil is defined node by node, by adding and subtracting two offsets.
2. **Normalization**: the gray of each step is taken to a scale from 0 (saturated) to 100
   (dry), clamping negatives to 0.
3. **Calibration**: that normalized value becomes saturation and water content through the
   curve of the soil.

An **incremental policy** may act on top of those three steps: the gray of a node never goes
back up, and once it passes the threshold the node is taken as saturated. It models a wetting
front that advances, and that is why it is **off by default**: it is a hypothesis about the
test, not a measurement, and with it on you cannot measure drying.

When enabled it acts on the normalized gray and not on the saturation, so that both fields
always come out of the same value. The original code applied it to the saturation only and
recomputed the water content from scratch at every step, which let the two fields contradict
each other: on the reference case it ended up with 634 nodes taken as saturated whose water
content was down to almost zero.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np

from .calibration import Calibration

log = logging.getLogger("pivnp")


@dataclass(frozen=True)
class References:
    """Gray value of the dry and of the saturated soil at every node."""

    dry: np.ndarray
    saturated: np.ndarray

    def __post_init__(self) -> None:
        if self.dry.shape != self.saturated.shape:
            raise ValueError("both references must have the same size")


def references_from_offsets(reference_gray: np.ndarray, dry_offset: float,
                            saturated_offset: float) -> References:
    """References from a single image, by adding and subtracting a fixed value.

    This is what MATLAB did (+5 and −6 gray levels). Those two numbers set the whole
    saturation scale, so they are worth checking against laboratory data.
    """
    if dry_offset <= saturated_offset:
        raise ValueError("the dry offset must be greater than the saturated one")
    gray = np.asarray(reference_gray, dtype=np.float64)
    return References(gray + dry_offset, gray + saturated_offset)


#: Below this width, one gray level weighs more than 5 % of the saturation scale. With the
#: 11 levels of the RGB flow it weighs 9 %, and there the method becomes very fragile.
MIN_BAND_WIDTH = 20.0


def warn_if_band_is_narrow(references: References, source: str = "") -> float:
    """Check the width of the band, the most sensitive number in the whole method.

    Returns the median width. The narrower it is, the more each gray level weighs: measured
    on the paper case, going from 40 levels to 11 moves the saturation by 0.105 on average,
    an order of magnitude more than any other decision in the computation.
    """
    width = references.dry - references.saturated
    finite = width[np.isfinite(width)]
    if finite.size == 0:
        return float("nan")
    median = float(np.median(finite))
    if median < MIN_BAND_WIDTH:
        log.warning("%sthe band between dry and saturated soil is %.0f gray levels, so one "
                    "level is %.0f %% of the saturation scale. With a band that narrow the "
                    "result depends heavily on image noise; it is worth revisiting it with "
                    "the calibration test of the soil.",
                    f"{source}: " if source else "", median, 100.0 / median)
    return median


def global_references(shape, dry_gray: float, saturated_gray: float) -> References:
    """References that are the same at every node, measured on the soil of the test.

    This is what the SWIR flow does: instead of taking the reference from the image node by
    node, two intensities are fixed for the whole material, that of the dry soil and that of
    the saturated one. The band then stops depending on the texture of each point.
    """
    if dry_gray <= saturated_gray:
        raise ValueError(f"the gray of the dry soil ({dry_gray}) must be greater than that "
                         f"of the saturated one ({saturated_gray})")
    return References(np.full(shape, float(dry_gray)), np.full(shape, float(saturated_gray)))


def normalize(gray: np.ndarray, references: References) -> np.ndarray:
    """Take the gray to the 0 (saturated) - 100 (dry) scale, clamping negatives."""
    span = references.dry - references.saturated
    with np.errstate(divide="ignore", invalid="ignore"):
        normalized = (np.asarray(gray, dtype=np.float64) - references.saturated) / span
    normalized = np.where(np.isfinite(normalized), normalized * 100.0, np.nan)
    return np.where(normalized < 0.0, 0.0, normalized)


#: Quality marks of every node, so that results can be read knowing what is a measurement.
MEASURED = 0  #: the gray falls inside the band: the value is a measurement
NO_DATA = 1  #: the node has no data (outside the image, or PIVlab did not measure it)
WET_LIMIT = 2  #: the gray falls below the saturated end: the value is an "at least"
DRY_LIMIT = 3  #: the gray falls above the dry end: the value is an "at most"


@dataclass
class State:
    """Result of one step."""

    saturation: np.ndarray
    moisture: np.ndarray
    #: Mark per node (:data:`MEASURED`, :data:`NO_DATA`, :data:`WET_LIMIT`,
    #: :data:`DRY_LIMIT`). A node at a limit is not a measurement, it is a bound: worth
    #: knowing when reading a field, because they tend to be many.
    quality: np.ndarray

    @property
    def clamped(self) -> int:
        """Nodes falling outside the band, whose value is therefore a bound."""
        return int(np.isin(self.quality, (WET_LIMIT, DRY_LIMIT)).sum())

    @property
    def measured(self) -> int:
        return int((self.quality == MEASURED).sum())


class MoistureModel:
    """Turn the gray of each step into degree of saturation and water content.

    The incremental policy acts **on the normalized gray**, not on the saturation: that way
    both fields come out of the same value and cannot contradict each other. The original
    code applied it to the saturation only and left the water content free, so a node could
    end up marked as saturated with its water content down to almost zero.
    """

    def __init__(self, calibration: Calibration, saturation_threshold: float = 0.95,
                 incremental: bool = False) -> None:
        if not 0 < saturation_threshold <= 1:
            raise ValueError(f"the threshold must be in (0, 1] and it is {saturation_threshold}")
        self.calibration = calibration
        self.threshold = float(saturation_threshold)
        self.incremental = bool(incremental)
        #: Lowest gray reached by each node: the soil does not dry, so it never goes back up.
        self.reached: np.ndarray | None = None

    def reset(self) -> None:
        self.reached = None

    def evaluate(self, normalized_gray: np.ndarray) -> State:
        """Saturation and water content of the step, applying the incremental policy."""
        gray = np.asarray(normalized_gray, dtype=np.float64)
        measured = np.isfinite(gray)
        if self.incremental:
            gray = self._apply_ratchet(gray, measured)
        result = self.calibration.evaluate(gray)
        return State(result.saturation, result.moisture, self._quality(gray, measured))

    # --- incremental policy ---------------------------------------------------------------
    def _apply_ratchet(self, gray: np.ndarray, measured: np.ndarray) -> np.ndarray:
        if self.reached is None:
            self.reached = np.full(gray.shape, np.inf)
        elif self.reached.shape != gray.shape:
            raise ValueError("the number of nodes changed between steps")

        gray = np.where(measured, np.minimum(gray, self.reached), gray)
        if self.threshold < 1.0:
            # Past the threshold the node is taken as saturated and never goes back.
            saturated = self.calibration.evaluate(gray).saturation >= self.threshold
            gray = np.where(measured & saturated, float(self.calibration.limits[0]), gray)
        self.reached = np.where(measured, gray, self.reached)
        return gray

    def _quality(self, gray: np.ndarray, measured: np.ndarray) -> np.ndarray:
        lowest, highest = self.calibration.limits
        quality = np.full(gray.shape, NO_DATA, dtype=np.int8)
        quality[measured] = MEASURED
        # normalize() already clamps from below, hence the <=
        quality[measured & (gray <= lowest)] = WET_LIMIT
        quality[measured & (gray >= highest)] = DRY_LIMIT
        return quality
