"""Reading the test images, channel selection and Gaussian filtering.

The filter reproduces MATLAB's (``imgaussfilt``), which is what makes the results comparable
with the earlier analyses: a kernel of size ``2·ceil(2σ)+1``, separable, with edge
replication padding.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
from numba import njit, prange

#: Channels accepted. The numbering is the one from the group's MATLAB code.
CHANNELS: dict[str, int] = {"red": 1, "green": 2, "blue": 3, "gray": 0}

#: Weights of MATLAB's ``rgb2gray`` (Rec. ITU-R BT.601).
GRAY_WEIGHTS = (0.298936021293776, 0.587043074451121, 0.114020904255103)


def channel_number(channel: int | str) -> int:
    """Translate ``1``/``"red"``/``"gray"``… into the internal number (0 = gray)."""
    if isinstance(channel, str):
        key = channel.strip().lower()
        if key not in CHANNELS:
            raise ValueError(f"unknown channel {channel!r}; use one of {sorted(CHANNELS)} "
                             "or the number 1, 2 or 3")
        return CHANNELS[key]
    if channel not in (0, 1, 2, 3):
        raise ValueError(f"unknown channel {channel!r}; use 1 (red), 2 (green), 3 (blue) "
                         "or 0 (gray)")
    return int(channel)


def matlab_rgb2gray(rgb: np.ndarray) -> np.ndarray:
    """Convert RGB to gray the way MATLAB's ``rgb2gray`` does, rounding to an integer."""
    weights = np.array(GRAY_WEIGHTS, dtype=np.float64)
    gray = rgb[:, :, :3].astype(np.float64) @ weights
    if rgb.dtype == np.uint8:
        return np.clip(np.floor(gray + 0.5), 0, 255).astype(np.uint8)
    return gray


def read_image(path: Path, channel: int | str = "gray") -> np.ndarray:
    """Read an image and return a single channel as a 2D array.

    ``channel`` may be 1 (red), 2 (green), 3 (blue) or 0/"gray" for the grayscale conversion.
    """
    try:
        from PIL import Image
    except ImportError:  # pragma: no cover - depends on the installation
        raise ImportError("reading images needs Pillow: pip install pillow") from None

    number = channel_number(channel)
    with Image.open(path) as image:
        array = np.asarray(image.convert("RGB"))
    return array[:, :, number - 1].copy() if number else matlab_rgb2gray(array)


def gaussian_kernel(sigma: float) -> np.ndarray:
    """MATLAB's 1D kernel: length ``2·ceil(2σ)+1``, normalized to sum 1."""
    if sigma <= 0:
        raise ValueError(f"sigma must be positive, and it is {sigma}")
    radius = int(math.ceil(2 * sigma))
    positions = np.arange(-radius, radius + 1, dtype=np.float64)
    kernel = np.exp(-(positions**2) / (2.0 * sigma * sigma))
    return kernel / kernel.sum()


@njit(parallel=True, cache=True)
def _convolve_horizontally(source, kernel, out):
    rows, columns = source.shape
    radius = kernel.size // 2
    for r in prange(rows):
        for c in range(columns):
            total = 0.0
            for k in range(-radius, radius + 1):
                origin = min(max(c + k, 0), columns - 1)  # edge replication padding
                total += source[r, origin] * kernel[k + radius]
            out[r, c] = total


@njit(parallel=True, cache=True)
def _convolve_vertically(source, kernel, out):
    rows, columns = source.shape
    radius = kernel.size // 2
    for r in prange(rows):
        for c in range(columns):
            total = 0.0
            for k in range(-radius, radius + 1):
                origin = min(max(r + k, 0), rows - 1)  # edge replication padding
                total += source[origin, c] * kernel[k + radius]
            out[r, c] = total


def gaussian_blur(image: np.ndarray, sigma: float, rounded: bool = False) -> np.ndarray:
    """Gaussian filter equivalent to MATLAB's ``imgaussfilt(image, sigma)``.

    The result comes back in floating point. MATLAB rounded it to the type of the image, that
    is to whole gray levels, and that is delicate here: the band between dry and saturated
    soil is narrow, so one level of difference turns into a large jump in saturation. With
    ``rounded=True`` the old behaviour is reproduced, to compare against older analyses.
    """
    kernel = gaussian_kernel(sigma)
    source = np.ascontiguousarray(image, dtype=np.float64)
    intermediate = np.empty_like(source)
    out = np.empty_like(source)
    _convolve_horizontally(source, kernel, intermediate)
    _convolve_vertically(intermediate, kernel, out)

    if rounded and np.issubdtype(image.dtype, np.integer):
        info = np.iinfo(image.dtype)
        return np.clip(np.floor(out + 0.5), info.min, info.max)
    return out
