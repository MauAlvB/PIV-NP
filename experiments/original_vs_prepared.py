"""Was the preprocessing already applied to these images worth it, and could it be better?

The ``Masked_*.jpg`` the analysis ran on carry 52 % more local contrast than the originals
they came from, so somebody enhanced them -- which is the right thing to do and is why
piling more preprocessing on top gained nothing.

The question left is the one that matters: starting from the untouched camera images, does
that enhancement improve the velocities, and is there a better way to do it? Same PIV, same
region, same pair of instants; only the image differs.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from PIL import Image

sys.path.insert(0, str(Path(__file__).parent))
from piv import analyse, anti_blocking, chain, high_pass, local_contrast  # noqa: E402
from preprocessing import neighbour_scatter, outlier_fraction  # noqa: E402

from test_images import masked, original as camera_image  # noqa: E402

REGION = (slice(200, 900), slice(400, 1500))


def load(path: Path) -> np.ndarray:
    return np.asarray(Image.open(path).convert("L"), dtype=np.float64)[REGION]


original = (load(camera_image("visual_375.jpg")), load(camera_image("visual_376.jpg")))
prepared = (load(masked(1)), load(masked(2)))

print(f"original : contrast {original[0].std():.1f}, mean {original[0].mean():.1f}")
print(f"prepared : contrast {prepared[0].std():.1f}, mean {prepared[0].mean():.1f}\n")

cases = [
    ("original, untouched", original, None),
    ("original + local contrast", original, local_contrast(32)),
    ("original + high-pass s=8", original, high_pass(8.0)),
    ("original + anti-blocking", original, anti_blocking()),
    ("original + anti-block + lc", original, chain(anti_blocking(), local_contrast(32))),
    ("as delivered (enhanced)", prepared, None),
    ("as delivered + anti-block", prepared, anti_blocking()),
]

header = f"{'image':<30}{'scatter px':>12}{'outliers':>10}{'coverage':>10}"
print(header)
print("-" * len(header))
for label, (first_raw, second_raw), preprocess in cases:
    first = first_raw if preprocess is None else preprocess(first_raw)
    second = second_raw if preprocess is None else preprocess(second_raw)
    u, v, _ = analyse(first, second, window=32, overlap=0.5, reference=first_raw)
    good = np.isfinite(u)
    print(f"{label:<30}{neighbour_scatter(u):>12.4f}{outlier_fraction(u, v):>9.1f}%"
          f"{100.0 * good.sum() / u.size:>9.1f}%")

print("\nThe two blocks are not strictly comparable on coverage: the delivered images have")
print("the background masked out, so fewer windows hold material. Scatter and outliers are")
print("measured only where there are vectors, so those do compare.")
