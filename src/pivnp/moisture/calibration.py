"""Calibration curve: from normalized gray to degree of saturation and water content.

The curve comes from laboratory tests and changes with every soil, so it is one more input
file. It is a CSV with three columns::

    gray,saturation,moisture
    100,0.0079,0.0890
    99.5,0.0074,0.1084
    ...

* ``gray``: gray value normalized between the saturated reference (0) and the dry one (100).
* ``saturation``: degree of saturation, between 0 and 1.
* ``moisture``: water content, in %.

Between the table points it interpolates the same way the original MATLAB code did
(``pchip``), so results stay comparable with the earlier analyses. Outside the table range
it does **not** extrapolate: it clamps to the nearest end and counts how many values were
clamped, because extrapolating was the source of impossible water contents (around −1700 %).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


def _pchip_slopes(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Fritsch-Carlson slopes, the ones MATLAB's ``pchip`` uses."""
    h = np.diff(x)
    delta = np.diff(y) / h
    d = np.zeros_like(y)
    if len(y) == 2:  # with two points MATLAB interpolates linearly
        return np.full(2, delta[0])

    interior = slice(1, len(y) - 1)
    same_sign = delta[:-1] * delta[1:] > 0
    w1 = 2 * h[1:] + h[:-1]
    w2 = h[1:] + 2 * h[:-1]
    with np.errstate(divide="ignore", invalid="ignore"):
        weighted = (w1 + w2) / (w1 / delta[:-1] + w2 / delta[1:])
    d[interior] = np.where(same_sign, weighted, 0.0)

    for end, (i, j) in ((0, (0, 1)), (len(y) - 1, (-1, -2))):
        slope = ((2 * h[i] + h[j]) * delta[i] - h[i] * delta[j]) / (h[i] + h[j])
        if slope * delta[i] <= 0:
            slope = 0.0
        elif delta[i] * delta[j] < 0 and abs(slope) > abs(3 * delta[i]):
            slope = 3 * delta[i]
        d[end] = slope
    return d


def pchip_interpolate(x: np.ndarray, y: np.ndarray, query: np.ndarray) -> np.ndarray:
    """Monotone piecewise cubic interpolation, compatible with ``interp1(..., 'pchip')``.

    ``x`` must be sorted in increasing order. Query values outside the range are evaluated
    with the cubic of the nearest end segment (clamping, if wanted, belongs to the caller;
    :class:`Calibration` does it).
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    query = np.asarray(query, dtype=np.float64)
    if x.ndim != 1 or x.size < 2 or y.shape != x.shape:
        raise ValueError("x and y must be vectors of the same size, with 2 points or more")
    if np.any(np.diff(x) <= 0):
        raise ValueError("x must be sorted in increasing order, with no repeated values")

    d = _pchip_slopes(x, y)
    i = np.clip(np.searchsorted(x, query) - 1, 0, x.size - 2)
    h = x[i + 1] - x[i]
    t = (query - x[i]) / h
    t2, t3 = t * t, t * t * t
    return ((2 * t3 - 3 * t2 + 1) * y[i]
            + (t3 - 2 * t2 + t) * h * d[i]
            + (-2 * t3 + 3 * t2) * y[i + 1]
            + (t3 - t2) * h * d[i + 1])


@dataclass(frozen=True)
class Result:
    """Degree of saturation and water content of a set of gray values."""

    saturation: np.ndarray
    moisture: np.ndarray
    clamped: int  # values that fell outside the calibration table


@dataclass(frozen=True)
class Calibration:
    """Calibration table of a soil, sorted by increasing gray value."""

    gray: np.ndarray
    saturation: np.ndarray
    moisture: np.ndarray
    source: str = "<memory>"

    def __post_init__(self) -> None:
        if self.gray.size < 2:
            raise ValueError(f"{self.source}: the calibration needs at least 2 points")
        if not (self.gray.shape == self.saturation.shape == self.moisture.shape):
            raise ValueError(f"{self.source}: the three columns must have the same size")
        if np.any(np.diff(self.gray) <= 0):
            raise ValueError(f"{self.source}: the gray column must increase, with no repeats")
        if np.any(~np.isfinite(self.gray)) or np.any(~np.isfinite(self.saturation)) \
                or np.any(~np.isfinite(self.moisture)):
            raise ValueError(f"{self.source}: the calibration admits no empty or NaN values")

    @property
    def limits(self) -> tuple[float, float]:
        return float(self.gray[0]), float(self.gray[-1])

    @classmethod
    def from_csv(cls, path: Path) -> Calibration:
        """Read the calibration from a CSV with columns ``gray,saturation,moisture``.

        Blank lines and lines starting with ``#`` are ignored, and the table may come in any
        order.
        """
        path = Path(path)
        rows = []
        for number, line in enumerate(path.read_text(encoding="utf-8-sig").splitlines(), 1):
            clean = line.strip()
            if not clean or clean.startswith("#"):
                continue
            parts = [p.strip() for p in clean.replace(";", ",").split(",")]
            if len(parts) < 3:
                raise ValueError(f"{path}:{number}: expected 3 columns, found {len(parts)}")
            if number == 1 or not _is_number(parts[0]):
                continue  # header
            try:
                rows.append([float(p) for p in parts[:3]])
            except ValueError:
                raise ValueError(f"{path}:{number}: non-numeric values ({clean!r})") from None
        if not rows:
            raise ValueError(f"{path}: contains no data rows")

        table = np.array(rows, dtype=np.float64)
        order = np.argsort(table[:, 0])
        return cls(table[order, 0], table[order, 1], table[order, 2], source=str(path))

    def evaluate(self, gray: np.ndarray) -> Result:
        """Saturation and water content of every gray value, clamping outside the table.

        Non-finite values (nodes without data) propagate as NaN and do not count as clamped.
        """
        gray = np.asarray(gray, dtype=np.float64)
        valid = np.isfinite(gray)
        lowest, highest = self.limits
        outside = valid & ((gray < lowest) | (gray > highest))
        bounded = np.clip(gray, lowest, highest)

        saturation = np.full(gray.shape, np.nan)
        moisture = np.full(gray.shape, np.nan)
        if valid.any():
            saturation[valid] = pchip_interpolate(self.gray, self.saturation, bounded[valid])
            moisture[valid] = pchip_interpolate(self.gray, self.moisture, bounded[valid])
        return Result(saturation, moisture, int(outside.sum()))


def _is_number(text: str) -> bool:
    try:
        float(text)
    except ValueError:
        return False
    return True
