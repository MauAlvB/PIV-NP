"""The test phase 1 still owed: a field with curvature.

Both cases with a known answer so far have a *linear* velocity field, and a symmetric kernel
reproduces a linear field exactly. That is part of why the edge-aware filters came out
looking perfect, so it proves less than it seems.

``horizontal_1P`` is the case that answers it: the imposed velocity goes 0, 0.1, 0.3, 0.9,
2.4 cm/s from one column of cells to the next, which is strongly curved, and it leaves a
different strain in every column (0.03, 0.06, 0.18...). Flattening peaks is exactly what a
smoothing kernel does to a field like that, so this is where the real cost shows up.

The strain of this case is also computed independently in ``tests/test_synthetic.py`` by
finite differences of the velocities, so the unfiltered run is not just a self-comparison.
"""
from __future__ import annotations

import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from filters import CANDIDATES, FilteredSource  # noqa: E402

from pivnp.config import load_case  # noqa: E402
from pivnp.simulation import Simulation  # noqa: E402
from pivnp.solver import output_mask  # noqa: E402
from pivnp.sources import build_source  # noqa: E402

CASE = Path(__file__).resolve().parents[1] / "tests" / "data" / "synthetic" / "horizontal_1P"


def run(velocity_filter) -> Simulation:
    work = Path(tempfile.mkdtemp())
    shutil.copytree(CASE, work, dirs_exist_ok=True, ignore=shutil.ignore_patterns("expected"))
    name, cfg = load_case(work)
    inner = build_source("pivlab", work, cfg, prefetch=0)
    sim = Simulation.from_directory(work)
    sim.frames = FilteredSource(inner, cfg.n_cols, cfg.n_rows, velocity_filter)
    sim.run()
    return sim


def profile(sim: Simulation) -> tuple[np.ndarray, np.ndarray]:
    """The strain of each column of cells, which is what the case imposes."""
    p = sim.particles
    shown = output_mask(p, sim.grid, sim.config.total_steps) & (p.nan_initial == 0)
    x = p.position[shown, 0]
    exx = p.strain[shown, 0]
    columns = np.round(x / sim.grid.dx).astype(int)
    values = np.array([exx[columns == c].mean() for c in np.unique(columns)])
    return np.unique(columns), values


truth_columns, truth = profile(run(CANDIDATES["none"]))
print("imposed strain per column of cells (the unfiltered run, which test_synthetic")
print("checks independently by finite differences):")
print("  " + "  ".join(f"{v:.4f}" for v in truth))
print(f"\n  peak {truth.max():.4f}   contrast between the two strongest columns "
      f"{abs(truth[-1] - truth[-2]):.4f}\n")

header = (f"{'filter':<22}{'peak':>10}{'peak lost':>11}{'contrast':>10}"
          f"{'contrast lost':>15}{'rms error':>11}")
print(header)
print("-" * len(header))

for label, velocity_filter in CANDIDATES.items():
    _, values = profile(run(velocity_filter))
    peak = values.max()
    contrast = abs(values[-1] - values[-2])
    rms = float(np.sqrt(np.mean((values - truth) ** 2)))
    print(f"{label:<22}{peak:>10.4f}{100 * (1 - peak / truth.max()):>10.1f}%"
          f"{contrast:>10.4f}"
          f"{100 * (1 - contrast / abs(truth[-1] - truth[-2])):>14.1f}%"
          f"{rms:>11.5f}")

print("\npeak lost / contrast lost: how much of the real feature the filter flattened.")
print("A smoothing kernel cannot both remove noise and keep a sharp gradient; this is")
print("where the trade stops being free.")
