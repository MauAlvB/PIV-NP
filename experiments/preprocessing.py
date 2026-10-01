"""Does preprocessing the images give PIV a better velocity field?

This is the question phase 1 could not reach. There the velocity field was taken as given
and filtered afterwards, which turned out to remove more signal than noise. Here the
velocities are produced from the photographs, so a change to the images is judged by the
field that comes out of them.

On a real pair there is no known answer, so what is measured is coherence:

* **peak quality** -- how far the correlation peak stands above the rest of the plane. The
  taller it is, the less room the displacement has to be found in the wrong place.
* **scatter** -- how much each vector differs from its neighbours. Real soil moves
  smoothly over a few millimetres, so most of this is noise.
* **outliers** -- the fraction of vectors the normalised median test rejects, which is the
  standard way PIV counts the ones that are simply wrong.
* **coverage** -- how many windows produced a vector at all; a filter that wins by refusing
  to answer has not won.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from PIL import Image

sys.path.insert(0, str(Path(__file__).parent))
from piv import CANDIDATES, analyse  # noqa: E402
from test_images import masked  # noqa: E402

PAIR = (masked(10), masked(11))
REGION = (slice(200, 900), slice(400, 1500))


def load(path: Path) -> np.ndarray:
    return np.asarray(Image.open(path).convert("L"), dtype=np.float64)[REGION]


def neighbour_scatter(field: np.ndarray) -> float:
    """How much a vector differs from the mean of the neighbours it has."""
    known = np.isfinite(field)
    values = np.where(known, field, 0.0)
    weight = known.astype(float)
    total = np.zeros_like(values)
    count = np.zeros_like(weight)
    for di, dj in ((0, 1), (0, -1), (1, 0), (-1, 0)):
        total += np.roll(np.roll(values, di, 0), dj, 1)
        count += np.roll(np.roll(weight, di, 0), dj, 1)
    around = np.where(count > 0, total / np.maximum(count, 1.0), np.nan)
    gap = field - around
    return float(np.nanstd(gap))


def outlier_fraction(u: np.ndarray, v: np.ndarray, threshold: float = 2.0) -> float:
    """Normalised median test (Westerweel & Scarano), the usual way to count bad vectors."""
    bad = 0
    total = 0
    rows, cols = u.shape
    for i in range(1, rows - 1):
        for j in range(1, cols - 1):
            if not np.isfinite(u[i, j]):
                continue
            total += 1
            for field in (u, v):
                ring = np.delete(field[i - 1:i + 2, j - 1:j + 2].ravel(), 4)
                ring = ring[np.isfinite(ring)]
                if ring.size < 4:
                    continue
                middle = np.median(ring)
                spread = np.median(np.abs(ring - middle)) + 0.1
                if abs(field[i, j] - middle) / spread > threshold:
                    bad += 1
                    break
    return 100.0 * bad / max(total, 1)


first_raw, second_raw = load(PAIR[0]), load(PAIR[1])
print(f"pair {PAIR[0].name} -> {PAIR[1].name}, region {first_raw.shape[1]}x"
      f"{first_raw.shape[0]} px\n")

header = (f"{'preprocessing':<24}{'peak quality':>14}{'scatter px':>12}"
          f"{'outliers':>10}{'coverage':>10}")
print(header)
print("-" * len(header))

baseline = None
for label, preprocess in CANDIDATES.items():
    first, second = preprocess(first_raw), preprocess(second_raw)
    # the raw image decides where the material is; the processed one is only correlated
    u, v, quality = analyse(first, second, window=32, overlap=0.5, reference=first_raw)
    good = np.isfinite(u)
    row = (float(np.nanmedian(quality)), neighbour_scatter(u),
           outlier_fraction(u, v), 100.0 * good.sum() / u.size)
    if baseline is None:
        baseline = row
    print(f"{label:<24}{row[0]:>14.3f}{row[1]:>12.4f}{row[2]:>9.1f}%{row[3]:>9.1f}%")

print(f"\n{'(baseline)':<24}{baseline[0]:>14.3f}{baseline[1]:>12.4f}"
      f"{baseline[2]:>9.1f}%{baseline[3]:>9.1f}%")
print("\npeak quality: higher is better. scatter and outliers: lower is better.")
print("coverage has to stay put: a filter that answers less is not answering better.")
