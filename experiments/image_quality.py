"""Are the images PIV ran on good material for it?

PIV correlates texture. What it needs is fine, high-contrast, locally distinctive detail;
what ruins it is anything that flattens that detail or adds structure of its own. The images
this test was analysed from, ``Masked_*.jpg``, are 177 KB for 1920x1080, while the originals
they came from are about 1 MB -- a five- to sixfold recompression, which falls exactly on
the detail the correlation lives off.

This measures both, with no PIV involved:

* **local contrast**: the standard deviation inside a window, averaged. The signal.
* **correlation sharpness**: how peaked a window's autocorrelation is. A flat autocorrelation
  means the window looks like its neighbours and the correlation peak will wander.
* **JPEG blocking**: how much the 8x8 grid of the compressor shows through. It is structure
  that does not move with the material, so it pulls the correlation towards zero displacement.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from PIL import Image

sys.path.insert(0, str(Path(__file__).parent))
from test_images import masked, original  # noqa: E402

PAIRS = (("as PIV saw it (recompressed)", masked(1)),
         ("the original", original("visual_375.jpg")))

WINDOW = 32          # a typical interrogation window
SAMPLES = 400


def load(path: Path) -> np.ndarray:
    image = Image.open(path).convert("L")
    return np.asarray(image, dtype=np.float64)


def windows(image: np.ndarray, rng: np.random.Generator) -> list[np.ndarray]:
    """Interrogation windows taken where there is actually material, not on the black mask."""
    out = []
    rows, cols = image.shape
    while len(out) < SAMPLES:
        r = int(rng.integers(0, rows - WINDOW))
        c = int(rng.integers(0, cols - WINDOW))
        patch = image[r:r + WINDOW, c:c + WINDOW]
        if patch.mean() > 25 and patch.std() > 1.0:     # skip the masked-out background
            out.append(patch)
    return out


def sharpness(patch: np.ndarray) -> float:
    """Height of the autocorrelation peak over its surroundings, as PIV would see it.

    The displacement is found by looking for a peak; the taller it stands over the rest, the
    less room there is for it to be found in the wrong place.
    """
    centred = patch - patch.mean()
    spectrum = np.fft.rfft2(centred)
    power = np.fft.irfft2(spectrum * np.conj(spectrum), s=patch.shape)
    power = np.fft.fftshift(power)
    peak = power.max()
    if peak <= 0:
        return 0.0
    middle = power.shape[0] // 2
    ring = power.copy()
    ring[middle - 2:middle + 3, middle - 2:middle + 3] = -np.inf
    return float(peak / max(ring.max(), 1e-12))


def blocking(image: np.ndarray) -> float:
    """How much the 8x8 JPEG grid shows: jumps on the block edges against jumps elsewhere."""
    diff = np.abs(np.diff(image, axis=1))
    columns = np.arange(diff.shape[1])
    on_edge = diff[:, (columns + 1) % 8 == 0].mean()
    inside = diff[:, (columns + 1) % 8 != 0].mean()
    return float(on_edge / max(inside, 1e-12))


rng = np.random.default_rng(0)
print(f"{'image':<32}{'local contrast':>16}{'peak sharpness':>16}{'JPEG blocking':>15}")
print("-" * 79)
results = {}
for label, path in PAIRS:
    image = load(path)
    patches = windows(image, np.random.default_rng(0))
    contrast = float(np.mean([p.std() for p in patches]))
    peak = float(np.median([sharpness(p) for p in patches]))
    block = blocking(image)
    results[label] = (contrast, peak, block)
    print(f"{label:<32}{contrast:>16.2f}{peak:>16.3f}{block:>15.3f}")

(c0, p0, b0), (c1, p1, b1) = results[PAIRS[0][0]], results[PAIRS[1][0]]
print(f"\nwhat the recompression cost: {100 * (1 - c0 / c1):.0f} % of the local contrast, "
      f"{100 * (1 - p0 / p1):.0f} % of the peak sharpness")
print(f"blocking: 1.00 would mean no JPEG grid at all; the originals are at {b1:.3f} "
      f"and the recompressed at {b0:.3f}")
