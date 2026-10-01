"""The synthetic cases said the median keeps features; the real case says it halves the peak.

The difference is how wide the feature is. A median removes anything narrower than about
half its window, and the columns of the synthetic case are many cells wide while a real
shear band can be one or two. This measures what each filter does to the peak of a real
case, which is the number the synthetic tests cannot give.

There is no ground truth here: a high strain on real data may be a real band or may be
noise. So this is read as a cost, not as an error -- how much of whatever is there survives.
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

CASE = Path(__file__).resolve().parents[1] / "examples" / "dam-break-swir"


def run(velocity_filter) -> tuple[np.ndarray, np.ndarray]:
    work = Path(tempfile.mkdtemp())
    shutil.copytree(CASE, work, dirs_exist_ok=True)
    name, cfg = load_case(work)
    inner = build_source("pivlab", work, cfg, prefetch=0)
    sim = Simulation.from_directory(work)
    source = FilteredSource(inner, cfg.n_cols, cfg.n_rows, velocity_filter)
    source.images = sim.frames.images
    source.moisture = False
    sim.frames = source
    sim.run()
    p = sim.particles
    shown = output_mask(p, sim.grid, sim.config.total_steps) & (p.nan_initial == 0)
    return p.eq_strain[shown], p.position[shown]


baseline, position = run(CANDIDATES["none"])
peak = float(baseline.max())
strong = baseline >= np.percentile(baseline, 95)
print(f"real case: {len(baseline)} particles, peak strain {peak:.4f}")
print(f"the strongest 5 % average {baseline[strong].mean():.4f}\n")

header = (f"{'filter':<22}{'peak':>9}{'kept':>8}{'top 5%':>9}{'kept':>8}"
          f"{'quiet zone':>12}{'cut':>8}")
print(header)
print("-" * len(header))

quiet_mask = baseline <= np.percentile(baseline, 25)
quiet_base = float(baseline[quiet_mask].mean())

for label, velocity_filter in CANDIDATES.items():
    values, _ = run(velocity_filter)
    top = float(values[strong].mean())
    quiet = float(values[quiet_mask].mean())
    print(f"{label:<22}{values.max():>9.4f}{100 * values.max() / peak:>7.0f}%"
          f"{top:>9.4f}{100 * top / baseline[strong].mean():>7.0f}%"
          f"{quiet:>12.5f}{100 * (1 - quiet / quiet_base):>7.0f}%")

print("\nkept: how much of the feature survives. cut: how much of the quiet-zone noise goes.")
print("A filter that cuts the noise and keeps the feature is doing its job; one that cuts")
print("both equally is just scaling the map down.")
