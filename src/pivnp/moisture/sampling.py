"""From the PIV grid to the image pixel, and sampling of the gray value.

The grid nodes come in metres in the PIVlab files; the conversion factor sits in the header
of those same files. The MATLAB code rounded the coordinates to the nearest pixel, and so
does this, so that the values match.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Registration:
    """Transform from the PIV grid to the moisture image.

    With a single camera it is the identity. With two (say a visible and an infrared one) it
    is the homography that makes both images coincide, the same one the MATLAB code built by
    marking points by hand on the two of them::

        | column |   | scale_x  shear_xy       offset_x |   | x |
        | row    | ~ | shear_yx scale_y        offset_y | · | y |
        | 1      |   | perspective_x perspective_y   1  |   | 1 |

    With both cameras in the same place the scale and the offset are enough; the other four
    values are needed when they look from different angles.
    """

    scale_x: float = 1.0
    offset_x: float = 0.0
    scale_y: float = 1.0
    offset_y: float = 0.0
    shear_xy: float = 0.0
    shear_yx: float = 0.0
    perspective_x: float = 0.0
    perspective_y: float = 0.0

    @property
    def is_identity(self) -> bool:
        return (self.scale_x, self.scale_y) == (1.0, 1.0) and not any(
            (self.offset_x, self.offset_y, self.shear_xy, self.shear_yx,
             self.perspective_x, self.perspective_y))


NO_REGISTRATION = Registration()


def pixel_coordinates(x_m: np.ndarray, y_m: np.ndarray, metres_per_pixel: float,
                      registration: Registration = NO_REGISTRATION,
                      ) -> tuple[np.ndarray, np.ndarray]:
    """Column and row (0-based) of every node inside the image.

    Metres become pixels, the registration is applied and the result is rounded to the
    nearest pixel, as MATLAB did; 1 is subtracted because indices start at 1 there.
    """
    if metres_per_pixel <= 0:
        raise ValueError(f"the conversion factor must be positive, and it is {metres_per_pixel}")
    x = np.asarray(x_m, dtype=np.float64) / metres_per_pixel
    y = np.asarray(y_m, dtype=np.float64) / metres_per_pixel
    if not registration.is_identity:
        weight = registration.perspective_x * x + registration.perspective_y * y + 1.0
        weight = np.where(np.abs(weight) < 1e-12, np.nan, weight)
        x, y = ((registration.scale_x * x + registration.shear_xy * y
                 + registration.offset_x) / weight,
                (registration.shear_yx * x + registration.scale_y * y
                 + registration.offset_y) / weight)
    column = np.floor(np.nan_to_num(x, nan=-1e9) + 0.5)
    row = np.floor(np.nan_to_num(y, nan=-1e9) + 0.5)
    return column.astype(np.int64) - 1, row.astype(np.int64) - 1


def sample(image: np.ndarray, column: np.ndarray, row: np.ndarray,
           has_data: np.ndarray | None = None) -> np.ndarray:
    """Gray value of the image at every node; NaN where there is no data or it falls outside.

    ``has_data`` marks the nodes PIVlab did measure; the rest are dropped, the same way
    MATLAB did with the mask of the ``.mat`` file.
    """
    height, width = image.shape
    inside = (column >= 0) & (column < width) & (row >= 0) & (row < height)
    usable = inside if has_data is None else inside & np.asarray(has_data, dtype=bool)

    gray = np.full(column.shape, np.nan)
    gray[usable] = image[row[usable], column[usable]]
    return gray


def outside_the_image(image: np.ndarray, column: np.ndarray, row: np.ndarray) -> int:
    """How many nodes fall outside the image (a sign of a badly fitted registration)."""
    height, width = image.shape
    return int((~((column >= 0) & (column < width) & (row >= 0) & (row < height))).sum())
