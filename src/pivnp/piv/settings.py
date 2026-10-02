"""Settings of the built-in PIV (``<case>.PIV``).

Same shape as the ``<case>.HUM`` of the moisture: every value is named, the order does not
matter, and ``!`` starts a comment::

    ! PIV settings
    IMAGES   = images/test_{n}.jpg   ! {n} is replaced by the image number
    SCALE    = 0.00054               ! metres per pixel
    WINDOW   = 32                    ! interrogation window, in pixels
    OVERLAP  = 0.5
    PASSES   = 2
    SMOOTH   = 0.6                   ! 0 keeps the raw field

Only ``IMAGES`` and ``SCALE`` have no sensible default. The scale is what turns a
displacement in pixels into one in metres, and there is no way to guess it: measure
something of known length in a photograph, or take it from the calibration of the test.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from ..moisture.images import channel_number

#: Default of every key. ``IMAGES`` and ``SCALE`` are required and have none.
DEFAULTS: dict[str, str] = {
    "IMAGES": "",
    "SCALE": "",
    "CHANNEL": "gray",
    "WINDOW": "32",
    "OVERLAP": "0.5",
    "PASSES": "2",
    "FIRST_IMAGE": "1",
    "OUTLIER_THRESHOLD": "2.0",
    "SUBPIXEL_OFFSET": "0",
    "SMOOTH": "0.6",
    "MASK_BELOW": "",
    "MASK_IMAGE": "",
    "REGION": "",
}


class PivSettingsError(ValueError):
    """The PIV settings file is missing or inconsistent."""


@dataclass(frozen=True)
class PivSettings:
    """How to measure the displacements of a test from its photographs."""

    image_pattern: str
    #: Metres per pixel. Everything downstream is in metres because of this one number.
    scale: float
    channel: int
    window: int
    overlap: float
    passes: int
    #: Number of the first photograph, for sequences that do not start at 1.
    first_image: int
    outlier_threshold: float
    #: Whether the passes after the first read their windows *between* the pixels of the
    #: photograph, applying the offset they were given in full instead of rounding it to a
    #: whole pixel. Off by default: it gives a four to five times better field per step and
    #: removes peak locking, but on the example whose accumulated answer is known it makes
    #: the accumulated shear worse, for a reason not yet found. See the correlation module.
    between_pixels: bool
    #: Width in grid points of the Gaussian smoothing of the finished field. The strain is a
    #: difference between neighbouring vectors, so the sub-pixel scatter that hardly shows in
    #: the displacement dominates it; 0 keeps the raw field.
    smoothing: float
    #: Pixels darker than this are not material. Useful when the background is masked out.
    mask_below: float | None
    #: An image marking the material instead, where anything non-black is material.
    mask_image: str
    #: ``left, top, right, bottom`` in pixels: the part of the photograph the grid covers.
    #: Without it the grid spans the whole image, which puts points on background that will
    #: never move and makes the grid bigger than it needs to be.
    region: tuple[int, int, int, int] | None
    source: str = "<memory>"
    unknown_keys: tuple[str, ...] = field(default_factory=tuple)

    def image_path(self, directory: Path, step: int) -> Path:
        """The photograph for a step. Step 1 is the first image, step 2 the second..."""
        return Path(directory) / self.image_pattern.format(n=self.first_image + step - 1)

    def validate(self) -> None:
        errors = []
        if not self.image_pattern:
            errors.append("IMAGES is missing, the name of the photographs")
        else:
            try:
                first = self.image_pattern.format(n=1)
                second = self.image_pattern.format(n=2)
            except (KeyError, IndexError, ValueError):
                errors.append("IMAGES must contain {n}, which is replaced by the image "
                              "number; {n:03d} pads it with zeros, as in image_001.jpg")
            else:
                if first == second:
                    errors.append("IMAGES must contain {n}, which is replaced by the image "
                                  "number; without it every step reads the same photograph")
        if self.scale <= 0:
            errors.append("SCALE is missing or not positive; it is the metres per pixel of "
                          "the photographs, and without it a displacement in pixels cannot "
                          "become one in metres")
        if self.window < 8:
            errors.append(f"WINDOW={self.window} is too small to correlate; 16 or 32 is usual")
        if self.window % 2:
            errors.append(f"WINDOW={self.window} should be even, so that the grid of windows "
                          "sits where the arithmetic says it does")
        if not 0.0 <= self.overlap < 1.0:
            errors.append(f"OVERLAP={self.overlap} must be between 0 and 1")
        if not 1 <= self.passes <= 4:
            errors.append(f"PASSES={self.passes} must be between 1 and 4")
        if self.outlier_threshold <= 0:
            errors.append("OUTLIER_THRESHOLD must be positive")
        if self.smoothing < 0:
            errors.append(f"SMOOTH={self.smoothing} cannot be negative; 0 turns it off")
        if self.smoothing > 3:
            errors.append(f"SMOOTH={self.smoothing} would blur over {2 * self.smoothing:.0f} "
                          "grid points either side and flatten the real strain; 0.6 is the "
                          "default and 2 is already a lot")
        if self.mask_below is not None and self.mask_image:
            errors.append("MASK_BELOW and MASK_IMAGE do the same job; give only one")
        if self.region is not None:
            left, top, right, bottom = self.region
            if right - left < self.window or bottom - top < self.window:
                errors.append(f"REGION is smaller than one window of {self.window} px")
            if left < 0 or top < 0:
                errors.append("REGION cannot start before the edge of the photograph")
        if errors:
            raise PivSettingsError(f"{self.source}: " + "; ".join(errors))


def parse(text: str, source: str = "<memory>") -> PivSettings:
    """Interpret the contents of a ``.PIV`` file."""
    values = dict(DEFAULTS)
    unknown = []
    for number, line in enumerate(text.splitlines(), 1):
        clean = line.split("!", 1)[0].strip()
        if not clean:
            continue
        if "=" not in clean:
            raise PivSettingsError(f"{source}:{number}: expected 'KEY = value' and read "
                                   f"{line.strip()!r}")
        key, value = (part.strip() for part in clean.split("=", 1))
        key = key.upper()
        if key not in DEFAULTS:
            unknown.append(key)
        values[key] = value

    def number(key: str, default: float = 0.0) -> float:
        raw = values[key].strip().replace(",", ".")
        if not raw:
            return default
        try:
            return float(raw)
        except ValueError:
            raise PivSettingsError(f"{source}: {key} must be a number and it is "
                                   f"{values[key]!r}") from None

    try:
        raw_channel = values["CHANNEL"]
        channel = channel_number(int(raw_channel) if raw_channel.lstrip("-").isdigit()
                                 else raw_channel)
    except ValueError as error:
        raise PivSettingsError(f"{source}: {error}") from None

    region = None
    if values["REGION"].strip():
        parts = [p for p in values["REGION"].replace(",", " ").split() if p]
        if len(parts) != 4:
            raise PivSettingsError(f"{source}: REGION needs four numbers, "
                                   "left top right bottom, in pixels")
        try:
            region = tuple(int(round(float(p))) for p in parts)
        except ValueError:
            raise PivSettingsError(f"{source}: REGION must be four numbers and it is "
                                   f"{values['REGION']!r}") from None

    settings = PivSettings(
        image_pattern=values["IMAGES"],
        scale=number("SCALE"),
        channel=channel,
        window=int(number("WINDOW", 32)),
        overlap=number("OVERLAP", 0.5),
        passes=int(number("PASSES", 2)),
        first_image=int(number("FIRST_IMAGE", 1)),
        outlier_threshold=number("OUTLIER_THRESHOLD", 2.0),
        between_pixels=number("SUBPIXEL_OFFSET", 0.0) != 0.0,
        smoothing=number("SMOOTH", 0.6),
        mask_below=number("MASK_BELOW") if values["MASK_BELOW"].strip() else None,
        mask_image=values["MASK_IMAGE"],
        region=region,
        source=source,
        unknown_keys=tuple(unknown),
    )
    settings.validate()
    return settings


def read(path: Path) -> PivSettings:
    """Read the ``.PIV`` file of a case."""
    path = Path(path)
    if not path.exists():
        raise PivSettingsError(
            f"{path}: this case is set to measure its own displacements, but there is no "
            f"{path.name} saying how. It needs at least the name of the photographs and the "
            f"scale in metres per pixel."
        )
    return parse(path.read_text(encoding="latin-1"), str(path))
