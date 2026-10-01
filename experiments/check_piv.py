"""Validate the small PIV before trusting it: shift a real image by a known amount.

A tool has to be checked before it is used to judge something else. A real photograph is
shifted by a known displacement, whole pixels and fractions of one, and the PIV is asked to
find it back. If it cannot, nothing measured with it afterwards means anything.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from PIL import Image

sys.path.insert(0, str(Path(__file__).parent))
from piv import analyse  # noqa: E402
from test_images import masked  # noqa: E402

IMAGE = masked(1)


def shifted(image: np.ndarray, dx: float, dy: float) -> np.ndarray:
    """Shift by a real number of pixels, interpolating linearly between the four neighbours."""
    rows, cols = image.shape
    y, x = np.mgrid[0:rows, 0:cols].astype(np.float64)
    sx, sy = x - dx, y - dy
    x0 = np.clip(np.floor(sx).astype(int), 0, cols - 2)
    y0 = np.clip(np.floor(sy).astype(int), 0, rows - 2)
    fx = np.clip(sx - x0, 0.0, 1.0)
    fy = np.clip(sy - y0, 0.0, 1.0)
    return ((1 - fx) * (1 - fy) * image[y0, x0] + fx * (1 - fy) * image[y0, x0 + 1]
            + (1 - fx) * fy * image[y0 + 1, x0] + fx * fy * image[y0 + 1, x0 + 1])


original = np.asarray(Image.open(IMAGE).convert("L"), dtype=np.float64)
crop = original[300:800, 600:1400]       # a region with material, away from the mask
print(f"region {crop.shape[1]} x {crop.shape[0]} px, "
      f"mean {crop.mean():.1f}, contrast {crop.std():.1f}\n")

print(f"{'imposed dx':>11}{'imposed dy':>11}{'measured dx':>13}{'measured dy':>13}"
      f"{'error':>9}{'windows':>9}")
print("-" * 66)

worst = 0.0
for dx, dy in ((2.0, 0.0), (0.0, 3.0), (2.5, -1.5), (0.37, 0.62), (-4.2, 2.8)):
    moved = shifted(crop, dx, dy)
    u, v, _ = analyse(crop, moved, window=32, overlap=0.5)
    good = np.isfinite(u)
    mu, mv = float(np.median(u[good])), float(np.median(v[good]))
    error = float(np.hypot(mu - dx, mv - dy))
    worst = max(worst, error)
    print(f"{dx:>11.2f}{dy:>11.2f}{mu:>13.3f}{mv:>13.3f}{error:>9.3f}{good.sum():>9}")

print(f"\nworst error: {worst:.3f} px")
print("Sub-pixel PIV is usually quoted at about 0.1 px; anything near that is fit for")
print("comparing one preprocessing against another, which is all it is for.")
