"""Does the finite measure buy its lower scatter with signal, the way the filters did?

Every filter in phase 1 cut the noise by removing the feature along with it: the best trade
on offer kept 73 % of the shear band to remove 26 % of the quiet-zone noise. The finite
strain reaches a similar scatter on the real case. The question is whether it pays the same
price, and there is no reason it should -- it is not smoothing anything, it is computing a
different and more correct quantity.
"""
from __future__ import annotations

import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))

from filters import CANDIDATES, FilteredSource  # noqa: E402
from finite_strain import (  # noqa: E402
    deformation_gradient,
    equivalent_shear,
    green_lagrange,
    lattice,
)
from test_images import DAM_BREAK, case_folder  # noqa: E402

from pivnp.config import load_case  # noqa: E402
from pivnp.simulation import Simulation  # noqa: E402
from pivnp.solver import output_mask  # noqa: E402
from pivnp.sources import build_source  # noqa: E402

CASE = case_folder(DAM_BREAK, "the dam-break case and its photographs",
                   Path(__file__).resolve().parents[1] / "examples" / "dam-break-swir")


def run(velocity_filter=None):
    work = Path(tempfile.mkdtemp())
    shutil.copytree(CASE, work, dirs_exist_ok=True)
    name, cfg = load_case(work)
    sim = Simulation.from_directory(work)
    if velocity_filter is not None:
        inner = build_source("pivlab", work, cfg, prefetch=0)
        source = FilteredSource(inner, cfg.n_cols, cfg.n_rows, velocity_filter)
        source.images = sim.frames.images
        source.moisture = False
        sim.frames = source
    sim.run()
    return sim


def scatter(values: np.ndarray, row: np.ndarray, col: np.ndarray) -> float:
    grid = np.full((row.max() + 1, col.max() + 1), np.nan)
    grid[row, col] = values
    return float(np.nanstd(np.concatenate([np.diff(grid, axis=a).ravel() for a in (0, 1)])))


base = run()
p = base.particles
usable = output_mask(p, base.grid, base.config.total_steps) & (p.nan_initial == 0)
row, col, _, _ = lattice(p.initial_position)

reference = p.eq_strain.copy()
strong = usable & (reference >= np.nanpercentile(reference[usable], 95))
quiet = usable & (reference <= np.nanpercentile(reference[usable], 25))
peak = float(np.nanmax(reference[usable]))
band = float(np.mean(reference[strong]))
noise = float(np.mean(reference[quiet]))
print(f"baseline: peak {peak:.4f}, strongest 5 % {band:.4f}, quiet quarter {noise:.5f}\n")

header = (f"{'method':<28}{'peak kept':>11}{'band kept':>11}{'quiet cut':>11}"
          f"{'scatter':>10}")
print(header)
print("-" * len(header))


def report(label: str, values: np.ndarray) -> None:
    print(f"{label:<28}{100 * np.nanmax(values[usable]) / peak:>10.0f}%"
          f"{100 * np.mean(values[strong]) / band:>10.0f}%"
          f"{100 * (1 - np.mean(values[quiet]) / noise):>10.0f}%"
          f"{scatter(np.where(usable, values, np.nan), row, col):>10.5f}")


report("incremental (as shipped)", reference)

exx, eyy, gxy = green_lagrange(*deformation_gradient(p, usable))
report("finite, same velocities", equivalent_shear(exx, eyy, gxy))

for label in ("edge: median 5x5", "edge: gaussian s=1.5", "edge: average 3x3"):
    filtered = run(CANDIDATES[label])
    report(f"incremental + {label}", filtered.particles.eq_strain)

print("\npeak and band kept: how much of the real feature survives, higher is better.")
print("quiet cut: how much of the noise where nothing moved goes, higher is better.")
print("A method that cuts the noise without touching the feature is not trading anything.")
