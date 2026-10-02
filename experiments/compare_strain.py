"""Incremental strain against finite strain, on cases whose answer is known.

Three questions, in order of how much they matter:

1. Does the finite measure stop inventing strain for a rigid rotation? The answer there is
   exactly zero, so there is nothing to argue about.
2. Does it still get a real deformation right? Simple shear has an exact finite answer too,
   and it is *not* the same as the incremental one: for engineering shear gamma the
   Green-Lagrange strain is E_xx = 0, E_yy = gamma^2/2, 2E_xy = gamma. The gamma^2/2 is a
   real term the linear theory drops, not an error.
3. Does it leave less noise on a real case?
"""
from __future__ import annotations

import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))

from finite_strain import (  # noqa: E402
    deformation_gradient,
    equivalent_shear,
    green_lagrange,
    rotation_degrees,
    volume_change,
)
from test_images import DAM_BREAK, case_folder  # noqa: E402

from pivnp.contour import make_contour_correction  # noqa: E402
from pivnp.simulation import Simulation  # noqa: E402
from pivnp.solver import output_mask  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
CASES = {
    "rotation_1P": (ROOT / "tests" / "data" / "synthetic" / "rotation_1P", 3),
    "shear-block": (ROOT / "examples" / "shear-block", None),
    "dam-break": (case_folder(DAM_BREAK, "the dam-break case and its photographs",
                              ROOT / "examples" / "dam-break-swir"), None),
}


def run(case: Path, contour: int | None):
    work = Path(tempfile.mkdtemp())
    shutil.copytree(case, work, dirs_exist_ok=True, ignore=shutil.ignore_patterns("expected"))
    sim = Simulation.from_directory(work)
    if contour is not None:
        sim.contour = make_contour_correction(contour)
    sim.run()
    return sim


def roughness(values: np.ndarray, row: np.ndarray, col: np.ndarray) -> float:
    """Scatter between neighbours on the initial lattice, which is a regular grid."""
    shape = (row.max() + 1, col.max() + 1)
    grid = np.full(shape, np.nan)
    grid[row, col] = values
    gaps = []
    for axis in (0, 1):
        gaps.append(np.diff(grid, axis=axis).ravel())
    return float(np.nanstd(np.concatenate(gaps)))


for name, (path, contour) in CASES.items():
    sim = run(path, contour)
    p = sim.particles
    usable = output_mask(p, sim.grid, sim.config.total_steps) & (p.nan_initial == 0)

    f = deformation_gradient(p, usable)
    exx, eyy, gxy = green_lagrange(*f)
    turn = rotation_degrees(*f)
    volume = volume_change(*f)
    eq_finite = equivalent_shear(exx, eyy, gxy)

    good = usable & np.isfinite(eq_finite)
    eq_incremental = p.eq_strain

    print(f"\n=== {name}  ({good.sum()} particles) ===")
    print(f"{'':<26}{'incremental':>14}{'finite':>14}")
    print(f"{'equivalent shear':<26}{np.mean(eq_incremental[good]):>14.6f}"
          f"{np.mean(eq_finite[good]):>14.6f}")
    print(f"{'volumetric':<26}{np.mean(p.vol_strain[good]):>14.6f}"
          f"{np.mean(volume[good]):>14.6f}")
    print(f"{'rotation (deg)':<26}{np.mean(p.rotation[good]):>14.4f}"
          f"{np.mean(turn[good]):>14.4f}")
    print(f"{'engineering shear':<26}{np.mean(p.strain[good, 2]):>14.6f}"
          f"{np.mean(gxy[good]):>14.6f}")

    from finite_strain import lattice
    row, col, _, _ = lattice(p.initial_position)
    print(f"{'scatter of eq. shear':<26}"
          f"{roughness(np.where(good, eq_incremental, np.nan), row, col):>14.6f}"
          f"{roughness(np.where(good, eq_finite, np.nan), row, col):>14.6f}")
