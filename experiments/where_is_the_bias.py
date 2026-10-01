"""Every filter broke the known answer. Where exactly does the damage happen?

Suspicion: not in the interior. A symmetric kernel reproduces a linear field exactly there,
whatever its width. At the edge of the material it has no neighbours on one side, so the
renormalised average leans inwards and flattens the very gradient we are measuring -- and
the edge is where PIV data is already worst.

If that is right, the bias should vanish as soon as the particles near the boundary are
left out of the average.
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

CASE = Path(__file__).resolve().parents[1] / "examples" / "shear-block"
TRUTH = 0.10
CELL = 0.01


def run(velocity_filter) -> Simulation:
    work = Path(tempfile.mkdtemp())
    shutil.copytree(CASE, work, dirs_exist_ok=True)
    name, cfg = load_case(work)
    inner = build_source("pivlab", work, cfg, prefetch=0)
    sim = Simulation.from_directory(work)
    sim.frames = FilteredSource(inner, cfg.n_cols, cfg.n_rows, velocity_filter)
    sim.run()
    return sim


header = f"{'filter':<20}{'all':>12}{'1 cell in':>12}{'2 cells in':>12}{'3 cells in':>12}"
print("mean gamma_xy, dropping the particles within N cells of the material edge")
print(f"(the truth is {TRUTH})\n")
print(header)
print("-" * len(header))

for label, velocity_filter in CANDIDATES.items():
    sim = run(velocity_filter)
    p = sim.particles
    shown = output_mask(p, sim.grid, sim.config.total_steps)
    x, y = p.position[:, 0], p.position[:, 1]
    # the block spans 0.01 to 0.13 in x and 0.01 to 0.09 in y
    row = [f"{label:<20}"]
    for margin in (0, 1, 2, 3):
        inside = shown & (x > 0.01 + margin * CELL) & (x < 0.13 - margin * CELL) \
            & (y > 0.01 + margin * CELL) & (y < 0.09 - margin * CELL)
        row.append(f"{np.mean(p.strain[inside, 2]):>12.6f}" if inside.sum()
                   else f"{'-':>12}")
    print("".join(row))
