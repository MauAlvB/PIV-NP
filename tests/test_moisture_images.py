"""Reading images, channel selection and Gaussian filtering."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pivnp.moisture.images import (
    channel_number,
    gaussian_blur,
    gaussian_kernel,
    matlab_rgb2gray,
    read_image,
)


# --- channels -----------------------------------------------------------------------------
@pytest.mark.parametrize(("given", "expected"), [
    (1, 1), (2, 2), (3, 3), (0, 0),
    ("red", 1), ("Green", 2), (" BLUE ", 3), ("gray", 0),
])
def test_channel_number(given, expected):
    assert channel_number(given) == expected


@pytest.mark.parametrize("given", [4, -1, "infrared"])
def test_unknown_channel(given):
    with pytest.raises(ValueError, match="unknown"):
        channel_number(given)


def test_gray_uses_the_matlab_weights():
    """``rgb2gray`` weighs 0.2989 R + 0.5870 G + 0.1140 B and rounds to an integer."""
    rgb = np.array([[[10, 200, 30]]], dtype=np.uint8)
    expected = round(0.298936021293776 * 10 + 0.587043074451121 * 200
                     + 0.114020904255103 * 30)
    assert matlab_rgb2gray(rgb)[0, 0] == expected == 124
    for tone in (0, 128, 255):  # pure grays are preserved
        assert matlab_rgb2gray(np.full((2, 2, 3), tone, dtype=np.uint8))[0, 0] == tone


def test_reads_every_channel(workdir: Path):
    pytest.importorskip("PIL")
    from PIL import Image
    data = np.zeros((4, 6, 3), dtype=np.uint8)
    data[..., 0], data[..., 1], data[..., 2] = 10, 200, 30
    path = workdir / "test.png"
    Image.fromarray(data).save(path)

    assert (read_image(path, 1) == 10).all()
    assert (read_image(path, "green") == 200).all()
    assert (read_image(path, 3) == 30).all()
    assert (read_image(path, "gray") == 124).all()
    assert read_image(path, 1).shape == (4, 6)


# --- kernel -------------------------------------------------------------------------------
@pytest.mark.parametrize(("sigma", "size"), [(0.5, 3), (1.0, 5), (15, 61), (40, 161)])
def test_kernel_size(sigma, size):
    """MATLAB uses 2·ceil(2σ)+1 points."""
    kernel = gaussian_kernel(sigma)
    assert kernel.size == size
    assert kernel.sum() == pytest.approx(1.0)
    np.testing.assert_allclose(kernel, kernel[::-1], atol=1e-15)
    assert np.argmax(kernel) == size // 2


def test_invalid_sigma():
    with pytest.raises(ValueError, match="positive"):
        gaussian_kernel(0)


# --- filter -------------------------------------------------------------------------------
def test_a_constant_image_does_not_change():
    """Checks both the kernel normalization and the edge replication padding."""
    image = np.full((30, 40), 173, dtype=np.uint8)
    np.testing.assert_allclose(gaussian_blur(image, 5.0), 173.0, atol=1e-9)


def test_matches_the_direct_convolution():
    """The separable version must give the same as convolving with the 2D kernel."""
    rng = np.random.default_rng(0)
    image = rng.uniform(0, 255, size=(37, 41))
    sigma = 3.0
    kernel = gaussian_kernel(sigma)
    radius = kernel.size // 2
    kernel2d = np.outer(kernel, kernel)

    rows, columns = image.shape
    expected = np.zeros_like(image)
    for r in range(rows):
        for c in range(columns):
            total = 0.0
            for i in range(-radius, radius + 1):
                for j in range(-radius, radius + 1):
                    ro = min(max(r + i, 0), rows - 1)
                    co = min(max(c + j, 0), columns - 1)
                    total += image[ro, co] * kernel2d[i + radius, j + radius]
            expected[r, c] = total
    np.testing.assert_allclose(gaussian_blur(image, sigma), expected, atol=1e-9)


def test_impulse_response():
    image = np.zeros((21, 21))
    image[10, 10] = 1.0
    filtered = gaussian_blur(image, 2.0)
    kernel = gaussian_kernel(2.0)
    np.testing.assert_allclose(filtered[10, 10], kernel[kernel.size // 2] ** 2, atol=1e-12)
    assert filtered.sum() == pytest.approx(1.0, abs=1e-9)
    np.testing.assert_allclose(filtered, filtered[::-1, :], atol=1e-12)


def test_returns_floating_point_without_rounding():
    """Rounding to whole levels amplified the noise, so it is no longer done."""
    rng = np.random.default_rng(1)
    image = rng.integers(0, 256, size=(20, 20), dtype=np.uint8)
    filtered = gaussian_blur(image, 2.0)
    assert filtered.dtype == np.float64
    assert np.any(filtered != np.round(filtered))


def test_the_legacy_mode_rounds():
    rng = np.random.default_rng(1)
    image = rng.integers(0, 256, size=(20, 20), dtype=np.uint8)
    exact = gaussian_blur(image, 2.0)
    legacy = gaussian_blur(image, 2.0, rounded=True)
    np.testing.assert_array_equal(legacy, np.clip(np.floor(exact + 0.5), 0, 255))
