"""Phase 0: the yardstick. For each filter, the bias and the noise, side by side.

    python experiments/measure.py

Bias is how much a known answer moves; noise is how much scatter is left in a real case.
A filter is only worth having if the second drops a lot and the first barely moves.
"""
from __future__ import annotations

import math
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from filters import CANDIDATES, FilteredSource  # noqa: E402

from pivnp.config import load_case  # noqa: E402
from pivnp.contour import make_contour_correction  # noqa: E402
from pivnp.simulation import Simulation  # noqa: E402
from pivnp.solver import output_mask  # noqa: E402
from pivnp.sources import build_source  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
SHEAR = ROOT / "examples" / "shear-block"
ROTATION = ROOT / "tests" / "data" / "synthetic" / "rotation_1P"
REAL = ROOT / "examples" / "dam-break-swir"

TRUE_SHEAR_STRAIN = 0.10
TRUE_EQ_STRAIN = 0.10 / math.sqrt(3.0)
TRUE_ROTATION = 50.0


def run(case: Path, velocity_filter, contour: int | None = None) -> Simulation:
    """Run a case with a filter between the files and the solver."""
    work = Path(tempfile.mkdtemp())
    shutil.copytree(case, work, dirs_exist_ok=True, ignore=shutil.ignore_patterns("expected"))
    name, cfg = load_case(work)
    inner = build_source("pivlab", work, cfg, prefetch=0)
    source = FilteredSource(inner, cfg.n_cols, cfg.n_rows, velocity_filter)
    sim = Simulation.from_directory(work)
    if cfg.moisture_from_images:
        source.images = sim.frames.images
        source.moisture = False
    sim.frames = source
    if contour is not None:
        sim.contour = make_contour_correction(contour)
    sim.run()
    return sim


def alive(sim: Simulation) -> np.ndarray:
    p = sim.particles
    return output_mask(p, sim.grid, sim.config.total_steps) & (p.nan_initial == 0)


def roughness(values: np.ndarray, position: np.ndarray, grid_step: float) -> float:
    """Scatter between each particle and the ones next to it: the short-wavelength part."""
    order = np.lexsort((position[:, 0], np.round(position[:, 1] / grid_step)))
    return float(np.std(np.diff(values[order])))


def main() -> None:
    header = (f"{'filter':<22}{'gamma_xy':>11}{'eq strain':>11}{'rotation':>10}"
              f"{'strain noise':>14}{'quiet zone':>12}")
    print(header)
    print("-" * len(header))

    # The quiet zone has to be the SAME particles for every filter, chosen on the unfiltered
    # run. Picking each filter's own lowest quarter instead compares a filter against itself,
    # and makes "it scaled everything down" look like "it removed noise".
    reference = run(REAL, CANDIDATES["none"])
    ref_mask = alive(reference)
    ref_strain = reference.particles.eq_strain[ref_mask]
    quiet = ref_strain <= np.percentile(ref_strain, 25)

    for label, velocity_filter in CANDIDATES.items():
        # --- bias: the cases whose answer is known
        shear = run(SHEAR, velocity_filter)
        mask = alive(shear)
        gxy = float(np.mean(shear.particles.strain[mask, 2]))
        eq = float(np.mean(shear.particles.eq_strain[mask]))

        turn = run(ROTATION, velocity_filter, contour=3)
        rotation = float(np.mean(np.abs(turn.particles.rotation[alive(turn)])))

        # --- noise: the real case
        real = run(REAL, velocity_filter)
        p = real.particles
        mask = alive(real)
        eq_real = p.eq_strain[mask]
        noise = roughness(eq_real, p.position[mask], real.grid.dy)
        quiet_level = float(np.mean(eq_real[quiet]))

        print(f"{label:<22}{gxy:>11.6f}{eq:>11.6f}{rotation:>10.3f}"
              f"{noise:>14.5f}{quiet_level:>12.5f}")

    print(f"\n{'truth':<22}{TRUE_SHEAR_STRAIN:>11.6f}{TRUE_EQ_STRAIN:>11.6f}"
          f"{TRUE_ROTATION:>10.3f}{'lower=better':>14}{'lower=better':>12}")
    print("\nstrain noise: scatter between neighbouring particles in the real case")
    print("quiet zone  : mean strain of the particles that deformed least WITHOUT a filter,")
    print("              the same ones for every row, which should be near zero")


if __name__ == "__main__":
    main()
