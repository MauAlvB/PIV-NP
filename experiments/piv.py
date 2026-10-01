"""A minimal PIV, so that preprocessing can be judged by what it does to the velocities.

PIVlab cannot be driven from here, and measuring an image does not say what the correlation
will make of it. This is the smallest thing that closes the loop: window cross-correlation
by FFT, with a Gaussian three-point fit for the sub-pixel part, which is what every PIV
package does at bottom.

It is not meant to compete with PIVlab -- no multi-pass, no window deformation, no
validation. It is meant to be the *same* algorithm applied to two versions of an image, so
that the difference between them is the preprocessing and nothing else.

It is also a first sketch of the built-in source PIV-NP could offer one day.
"""
from __future__ import annotations

from collections.abc import Callable

import numpy as np

#: A preprocessing step: an image in, an image out, same shape.
Preprocess = Callable[[np.ndarray], np.ndarray]


def correlate(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Circular cross-correlation of two windows, peak centred."""
    a = a - a.mean()
    b = b - b.mean()
    spectrum = np.fft.rfft2(b) * np.conj(np.fft.rfft2(a))
    return np.fft.fftshift(np.fft.irfft2(spectrum, s=a.shape))


def _sub_pixel(line: np.ndarray, peak: int) -> float:
    """Gaussian three-point fit, the standard estimator; falls back to the integer peak."""
    if peak <= 0 or peak >= line.size - 1:
        return 0.0
    left, middle, right = line[peak - 1], line[peak], line[peak + 1]
    if left <= 0 or middle <= 0 or right <= 0:
        return 0.0
    low, mid, high = np.log(left), np.log(middle), np.log(right)
    bottom = 2.0 * (low - 2.0 * mid + high)
    return 0.0 if abs(bottom) < 1e-12 else float((low - high) / bottom)


def displacement(a: np.ndarray, b: np.ndarray) -> tuple[float, float, float]:
    """Displacement of one window, in pixels, and how far the peak stands above the rest."""
    power = correlate(a, b)
    flat = int(np.argmax(power))
    row, col = np.unravel_index(flat, power.shape)
    middle_r, middle_c = power.shape[0] // 2, power.shape[1] // 2

    dr = _sub_pixel(power[:, col], row)
    dc = _sub_pixel(power[row, :], col)

    ring = power.copy()
    r0, r1 = max(0, row - 2), min(power.shape[0], row + 3)
    c0, c1 = max(0, col - 2), min(power.shape[1], col + 3)
    ring[r0:r1, c0:c1] = -np.inf
    quality = float(power[row, col] / max(ring.max(), 1e-12))
    return float(col - middle_c + dc), float(row - middle_r + dr), quality


def analyse(first: np.ndarray, second: np.ndarray, window: int = 32,
            overlap: float = 0.5, floor: float = 25.0, reference: np.ndarray | None = None):
    """Velocity field of a pair of images, on a regular grid of interrogation windows.

    Windows whose material is mostly the masked-out background come back as NaN, the same
    way PIVlab leaves them, so that what follows sees the edge of the material where it is.

    ``reference`` is the image that decides which windows hold material, and it has to be
    the **unprocessed** one: after a high-pass the mean is zero everywhere, so judging by
    the processed image throws away the whole field.
    """
    if reference is None:
        reference = first
    step = max(1, int(window * (1.0 - overlap)))
    rows = range(0, first.shape[0] - window + 1, step)
    cols = range(0, first.shape[1] - window + 1, step)
    u = np.full((len(rows), len(cols)), np.nan)
    v = np.full((len(rows), len(cols)), np.nan)
    quality = np.full((len(rows), len(cols)), np.nan)

    for i, r in enumerate(rows):
        for j, c in enumerate(cols):
            judge = reference[r:r + window, c:c + window]
            if judge.mean() < floor or judge.std() < 1.0:
                continue
            a = first[r:r + window, c:c + window]
            b = second[r:r + window, c:c + window]
            du, dv, peak = displacement(a, b)
            u[i, j], v[i, j], quality[i, j] = du, dv, peak
    return u, v, quality


# --- preprocessing candidates ------------------------------------------------------------
def as_is(image: np.ndarray) -> np.ndarray:
    return image


def _blur(image: np.ndarray, sigma: float) -> np.ndarray:
    """Separable Gaussian, reflecting at the border."""
    half = max(1, int(np.ceil(3 * sigma)))
    offsets = np.arange(-half, half + 1)
    kernel = np.exp(-0.5 * (offsets / sigma) ** 2)
    kernel /= kernel.sum()
    padded = np.pad(image, half, mode="reflect")
    out = np.apply_along_axis(lambda line: np.convolve(line, kernel, "valid"), 1, padded)
    return np.apply_along_axis(lambda line: np.convolve(line, kernel, "valid"), 0, out)


def high_pass(sigma: float = 8.0) -> Preprocess:
    """Remove what varies slowly: lighting, shading, the background.

    None of it moves with the material, so all it does is weigh on the correlation.
    """
    def apply(image: np.ndarray) -> np.ndarray:
        return image - _blur(image, sigma)

    return apply


def local_contrast(window: int = 32) -> Preprocess:
    """Divide each neighbourhood by its own spread, so that dim areas count as much as bright.

    This is what CLAHE is for; done with a box here to keep it simple.
    """
    def apply(image: np.ndarray) -> np.ndarray:
        sigma = window / 4.0
        mean = _blur(image, sigma)
        spread = np.sqrt(np.maximum(_blur((image - mean) ** 2, sigma), 1e-6))
        return (image - mean) / np.maximum(spread, 1.0)

    return apply


def cap_intensity(percentile: float = 98.0) -> Preprocess:
    """Clip the brightest pixels, which otherwise dominate the correlation on their own."""
    def apply(image: np.ndarray) -> np.ndarray:
        return np.minimum(image, np.percentile(image, percentile))

    return apply


def anti_blocking() -> Preprocess:
    """Damp the JPEG 8x8 grid, which is locked to the pixels and does not move with the soil.

    Measured on this test's images: the gradients on the block edges are 37 % stronger than
    inside them. The correlation sees that grid in both images at the same place, which pulls
    the peak towards zero displacement.
    """
    def apply(image: np.ndarray) -> np.ndarray:
        out = image.astype(np.float64).copy()
        for axis in (0, 1):
            moved = np.moveaxis(out, axis, 0)
            edges = np.arange(7, moved.shape[0] - 1, 8)
            moved[edges] = 0.5 * (moved[edges - 1] + moved[edges + 1])
            out = np.moveaxis(moved, 0, axis)
        return out

    return apply


def chain(*steps: Preprocess) -> Preprocess:
    def apply(image: np.ndarray) -> np.ndarray:
        for step in steps:
            image = step(image)
        return image

    return apply


CANDIDATES: dict[str, Preprocess] = {
    "none": as_is,
    "high-pass s=8": high_pass(8.0),
    "high-pass s=4": high_pass(4.0),
    "local contrast": local_contrast(32),
    "cap at 98%": cap_intensity(98.0),
    "anti-blocking": anti_blocking(),
    "anti-block + high-pass": chain(anti_blocking(), high_pass(8.0)),
}
