"""Settings file of the moisture measurement (``<case>.HUM``).

Every value is identified by its name, the order does not matter and comments start with
``!`` as in Fortran::

    ! moisture measurement settings
    IMAGES      = vis_{n}.jpg      ! {n} is replaced by the step number
    CHANNEL     = gray             ! 1 red, 2 green, 3 blue, 0 or "gray" for grayscale
    SIGMA       = 40               ! averaging radius, in pixels
    DRY_REFERENCE = ref2.jpg
    CALIBRATION = calibration_slope_rgb.csv

Keys that do not appear take their default value, which is the one of the original MATLAB
code except where a different criterion was agreed: ``FIRST_STEP``, ``LEGACY_ROUNDING`` and
``INCREMENTAL`` with its ``SATURATION_THRESHOLD``, which are off because they are a
hypothesis about the test rather than a measurement. The old behaviour comes back by setting
them to ``legacy`` and to ``1``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from .images import channel_number
from .sampling import Registration

#: Default value of every key.
DEFAULTS: dict[str, str] = {
    "IMAGES": "vis_{n}.jpg",
    "CHANNEL": "gray",
    "SIGMA": "40",
    "DRY_REFERENCE": "ref2.jpg",
    "DRY_OFFSET": "5",
    "SATURATED_OFFSET": "-6",
    "DRY_BAND": "",
    "SATURATED_BAND": "",
    "CALIBRATION": "",
    "SATURATION_THRESHOLD": "0.95",
    "INCREMENTAL": "0",
    "FIRST_STEP": "same",
    "LEGACY_ROUNDING": "0",
    "SCALE_X": "1.0",
    "OFFSET_X": "0.0",
    "SCALE_Y": "1.0",
    "OFFSET_Y": "0.0",
    "SHEAR_XY": "0.0",
    "SHEAR_YX": "0.0",
    "PERSPECTIVE_X": "0.0",
    "PERSPECTIVE_Y": "0.0",
}

#: What to do with the first step: like the rest, or like MATLAB (water content 0 and
#: constant saturation).
FIRST_STEP = ("same", "legacy")


class SettingsError(ValueError):
    """The moisture settings file is missing or inconsistent."""


@dataclass(frozen=True)
class MoistureSettings:
    """Parameters of the moisture measurement of a test."""

    image_pattern: str
    channel: int
    sigma: float
    dry_reference: str
    dry_offset: float
    saturated_offset: float
    #: Intensities of the dry and of the saturated soil, the same for the whole image. When
    #: given they replace the per-node reference; this is what the SWIR flow does.
    dry_band: float | None
    saturated_band: float | None
    calibration: str
    saturation_threshold: float
    incremental: bool
    first_step: str
    legacy_rounding: bool
    scale_x: float
    offset_x: float
    scale_y: float
    offset_y: float
    shear_xy: float = 0.0
    shear_yx: float = 0.0
    perspective_x: float = 0.0
    perspective_y: float = 0.0
    source: str = "<memory>"
    unknown_keys: tuple[str, ...] = field(default_factory=tuple)

    def image_path(self, directory: Path, step: int) -> Path:
        return Path(directory) / self.image_pattern.format(n=step)

    @property
    def has_global_band(self) -> bool:
        """Whether the dry-saturated band is fixed with two intensities for the whole image."""
        return self.dry_band is not None and self.saturated_band is not None

    @property
    def registration(self) -> Registration:
        """Transform from the PIV grid to the moisture image."""
        return Registration(self.scale_x, self.offset_x, self.scale_y, self.offset_y,
                            self.shear_xy, self.shear_yx,
                            self.perspective_x, self.perspective_y)

    @property
    def registration_is_identity(self) -> bool:
        return self.registration.is_identity

    def validate(self) -> None:
        errors = []
        if "{n}" not in self.image_pattern:
            errors.append("IMAGES must contain {n}, which is replaced by the step number")
        if self.sigma <= 0:
            errors.append(f"SIGMA={self.sigma} must be positive")
        if not 0 < self.saturation_threshold <= 1:
            errors.append(f"SATURATION_THRESHOLD={self.saturation_threshold} must be in (0, 1]")
        if self.dry_offset <= self.saturated_offset:
            errors.append("DRY_OFFSET must be greater than SATURATED_OFFSET")
        if (self.dry_band is None) != (self.saturated_band is None):
            errors.append("DRY_BAND and SATURATED_BAND go together: either both or neither")
        elif self.has_global_band and self.dry_band <= self.saturated_band:
            errors.append(f"DRY_BAND={self.dry_band} must be greater than "
                          f"SATURATED_BAND={self.saturated_band}")
        if self.first_step not in FIRST_STEP:
            errors.append(f"FIRST_STEP={self.first_step!r} must be "
                          f"{' or '.join(FIRST_STEP)}")
        if self.scale_x == 0 or self.scale_y == 0:
            errors.append("SCALE_X and SCALE_Y cannot be zero")
        if not self.calibration:
            errors.append("CALIBRATION is missing, the file with the curve of the soil")
        if errors:
            raise SettingsError(f"{self.source}: " + "; ".join(errors))


def parse(text: str, source: str = "<memory>") -> MoistureSettings:
    """Interpret the contents of a ``.HUM`` file."""
    values = dict(DEFAULTS)
    unknown = []
    for number, line in enumerate(text.splitlines(), 1):
        clean = line.split("!", 1)[0].strip()
        if not clean:
            continue
        if "=" not in clean:
            raise SettingsError(f"{source}:{number}: expected 'KEY = value' and read "
                                f"{line.strip()!r}")
        key, value = (part.strip() for part in clean.split("=", 1))
        key = key.upper()
        if key not in DEFAULTS:
            unknown.append(key)
        values[key] = value

    def real_number(key: str) -> float:
        try:
            return float(values[key].replace(",", "."))
        except ValueError:
            raise SettingsError(f"{source}: {key} must be a number and it is "
                                f"{values[key]!r}") from None

    def affirmative(key: str) -> bool:
        return values[key].strip().lower() in ("1", "yes", "true")

    def number_or_nothing(key: str) -> float | None:
        return real_number(key) if values[key].strip() else None

    try:
        channel = channel_number(values["CHANNEL"] if not values["CHANNEL"].lstrip("-").isdigit()
                                 else int(values["CHANNEL"]))
    except ValueError as error:
        raise SettingsError(f"{source}: {error}") from None

    settings = MoistureSettings(
        image_pattern=values["IMAGES"],
        channel=channel,
        sigma=real_number("SIGMA"),
        dry_reference=values["DRY_REFERENCE"],
        dry_offset=real_number("DRY_OFFSET"),
        saturated_offset=real_number("SATURATED_OFFSET"),
        dry_band=number_or_nothing("DRY_BAND"),
        saturated_band=number_or_nothing("SATURATED_BAND"),
        calibration=values["CALIBRATION"],
        saturation_threshold=real_number("SATURATION_THRESHOLD"),
        incremental=affirmative("INCREMENTAL"),
        first_step=values["FIRST_STEP"].strip().lower(),
        legacy_rounding=affirmative("LEGACY_ROUNDING"),
        scale_x=real_number("SCALE_X"),
        offset_x=real_number("OFFSET_X"),
        scale_y=real_number("SCALE_Y"),
        offset_y=real_number("OFFSET_Y"),
        shear_xy=real_number("SHEAR_XY"),
        shear_yx=real_number("SHEAR_YX"),
        perspective_x=real_number("PERSPECTIVE_X"),
        perspective_y=real_number("PERSPECTIVE_Y"),
        source=source,
        unknown_keys=tuple(unknown),
    )
    settings.validate()
    return settings


def read(path: Path) -> MoistureSettings:
    """Read the ``.HUM`` file of a case."""
    path = Path(path)
    if not path.exists():
        raise SettingsError(f"{path}: the moisture settings file does not exist")
    return parse(path.read_text(encoding="latin-1"), str(path))
