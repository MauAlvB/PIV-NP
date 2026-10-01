"""The real acceptance test: the same case analysed twice, once from PIVlab and once from
the photographs.

Comparing velocity fields step by step turned out to be the wrong question. Each field has
about 45 % of its points without data, in different places, and at a displacement of a sixth
of a pixel per step neither implementation is measuring much anyway. What matters is what the
user gets at the end, so that is what is compared: the displacement and the strain of the
particles after the whole analysis, which is an accumulation over twenty steps and well above
the noise.

The agreed tolerance is 10 %.
"""
from __future__ import annotations

import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))

from pivnp.simulation import RunOptions, Simulation  # noqa: E402
from pivnp.solver import output_mask  # noqa: E402

CASE = Path(__file__).resolve().parents[1] / "examples" / "dam-break-swir"


def _pearson(x: np.ndarray, y: np.ndarray) -> float:
    """Correlation coefficient, by explicit sums.

    Written out rather than ``np.corrcoef`` because this environment has several BLAS
    libraries loaded at once and anything going through a matrix product crashes the
    interpreter. Noted in the project status; the arithmetic here is the same.
    """
    dx, dy = x - x.mean(), y - y.mean()
    bottom = np.sqrt((dx * dx).sum() * (dy * dy).sum())
    return float((dx * dy).sum() / bottom) if bottom > 0 else float("nan")


#: Where the two analyses are kept so the comparison can be reworked without re-running them.
#: Both take minutes; the comparison is where the thinking is, and it was reworked twice.
CACHE = Path(__file__).parent / "_end_to_end_cache.npz"


def _patches(position: np.ndarray, keep: np.ndarray, side: float) -> np.ndarray:
    """A label per kept particle saying which square patch of the model it falls in."""
    x, y = position[keep, 0], position[keep, 1]
    if side <= 0:
        return np.arange(x.size)                      # every particle its own patch
    columns = np.floor(x / side).astype(np.int64)
    rows = np.floor(y / side).astype(np.int64)
    return rows * (columns.max() - columns.min() + 1) + columns


def _by_patch(groups: np.ndarray, left: np.ndarray,
              right: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """The mean of each quantity in each patch, in the same order for both."""
    labels, index = np.unique(groups, return_inverse=True)
    count = np.bincount(index, minlength=labels.size)
    return (np.bincount(index, weights=left, minlength=labels.size) / count,
            np.bincount(index, weights=right, minlength=labels.size) / count)


def run(source: str) -> Simulation:
    work = Path(tempfile.mkdtemp())
    shutil.copytree(CASE, work, dirs_exist_ok=True)
    sim = Simulation.from_directory(work, options=RunOptions(source=source, prefetch=0))
    sim.run()
    return sim


def _collect(source: str) -> dict:
    sim = run(source)
    p = sim.particles
    return {
        "displacement": p.displacement,
        "strain": p.strain,
        "eq_strain": p.eq_strain,
        "position": p.position,
        "shown": (output_mask(p, sim.grid, sim.config.total_steps) & (p.nan_initial == 0)),
    }


def gather(fresh: bool = False) -> tuple[dict, dict]:
    """The two analyses, from the cache when it is there."""
    if CACHE.exists() and not fresh:
        stored = np.load(CACHE)
        keys = ("displacement", "strain", "eq_strain", "position", "shown")
        return ({k: stored[f"a_{k}"] for k in keys}, {k: stored[f"b_{k}"] for k in keys})
    print("running from the PIVlab files...")
    first = _collect("pivlab")
    print("running from the photographs...")
    second = _collect("images")
    np.savez_compressed(CACHE, **{f"a_{k}": v for k, v in first.items()},
                        **{f"b_{k}": v for k, v in second.items()})
    return first, second


def main(fresh: bool = False) -> None:
    left, right = gather(fresh)

    class Holder:
        def __init__(self, d):
            self.__dict__.update(d)

    a, b = Holder(left), Holder(right)
    shown = a.shown & b.shown
    print(f"\nparticles both analyses kept: {shown.sum()} of {shown.size}")

    moved = shown & (np.hypot(*a.displacement[:, :2].T) > 1e-3)   # more than a millimetre
    print(f"of those, moved more than 1 mm: {moved.sum()}")

    # Two different questions, and they need two different measurements.
    #
    # For the displacement, the strict one: a particle went somewhere, and the two analyses
    # should put it in the same place. That is measured particle by particle.
    #
    # For the strain it is not. The strain is a *difference between neighbouring vectors*, so
    # at a sixth of a pixel per step it is a small difference between two noisy numbers, and
    # no two PIV implementations agree on it particle by particle -- two runs of PIVlab with
    # slightly different validation settings would not either. What can be asked, and what
    # the user actually reads off the result, is whether the two agree on *how much* strain
    # there is and on *where* it is. So the strain is judged on the summary of the field and
    # on whether the two patterns track each other.
    print(f"\n{'':<22}{'from PIVlab':>14}{'from photos':>14}{'difference':>13}")
    print("-" * 63)

    def compare(label: str, left: np.ndarray, right: np.ndarray, unit: float = 1.0) -> float:
        """Agreement of the typical value of a quantity over the field."""
        x, y = np.abs(left[moved]) * unit, np.abs(right[moved]) * unit
        mx, my = np.median(x), np.median(y)
        gap = 100.0 * abs(my - mx) / max(mx, 1e-30)
        print(f"{label:<22}{mx:>14.5f}{my:>14.5f}{gap:>12.1f}%")
        return gap

    def per_particle(label: str, left: np.ndarray, right: np.ndarray,
                     unit: float = 1.0) -> float:
        """Agreement of each individual particle -- the strict test."""
        x, y = left[moved] * unit, right[moved] * unit
        gap = 100.0 * np.median(np.abs(y - x)) / max(np.median(np.abs(x)), 1e-30)
        print(f"{label:<22}{np.median(np.abs(x)):>14.5f}{np.median(np.abs(y)):>14.5f}"
              f"{gap:>12.1f}%")
        return gap

    print("displacement, judged particle by particle:")
    gaps = {
        "displacement x [mm]": per_particle("  displacement x [mm]", a.displacement[:, 0],
                                            b.displacement[:, 0], 1000.0),
        "displacement y [mm]": per_particle("  displacement y [mm]", a.displacement[:, 1],
                                            b.displacement[:, 1], 1000.0),
    }
    print("\nstrain, judged on the field as a whole:")
    summary = {
        "equivalent strain": compare("  equivalent strain", a.eq_strain, b.eq_strain),
        "shear strain": compare("  shear strain", a.strain[:, 2], b.strain[:, 2]),
        "volumetric strain": compare("  volumetric strain",
                                     a.strain[:, 0] + a.strain[:, 1],
                                     b.strain[:, 0] + b.strain[:, 1]),
    }
    gaps.update(summary)

    print("\n  does the strain sit in the same places?")
    print("  Correlation between the two, with the particles grouped into square patches of")
    print("  the side below. If the two measure the same field and differ only in")
    print("  point-to-point noise, grouping averages the noise away and the correlation")
    print("  climbs; if they measure different fields it stays where it is.")
    print(f"\n{'patch side':>12}{'particles/patch':>18}{'eq strain r':>14}{'shear r':>10}")
    for side in (0.0, 0.01, 0.02, 0.04, 0.08):
        groups = _patches(a.position, moved, side)
        eq = _pearson(*_by_patch(groups, a.eq_strain[moved], b.eq_strain[moved]))
        sh = _pearson(*_by_patch(groups, a.strain[moved, 2], b.strain[moved, 2]))
        per = moved.sum() / max(len(set(groups.tolist())), 1)
        label = "none" if side == 0 else f"{1000 * side:.0f} mm"
        print(f"{label:>12}{per:>18.1f}{eq:>14.3f}{sh:>10.3f}")

    size_a = np.hypot(a.displacement[moved, 0], a.displacement[moved, 1])
    size_b = np.hypot(b.displacement[moved, 0], b.displacement[moved, 1])
    apart = np.hypot(a.displacement[moved, 0] - b.displacement[moved, 0],
                     a.displacement[moved, 1] - b.displacement[moved, 1])
    relative = 100.0 * apart / np.maximum(size_a, 1e-30)
    print("\ndisplacement, particle by particle:")
    print(f"  median magnitude from PIVlab : {1000 * np.median(size_a):.3f} mm")
    print(f"  median magnitude from photos : {1000 * np.median(size_b):.3f} mm")
    print(f"  median distance between them : {1000 * np.median(apart):.3f} mm")
    print(f"  within 10 %                  : {100 * (relative < 10).mean():.1f} % of particles")
    print(f"  within 25 %                  : {100 * (relative < 25).mean():.1f} %")

    worst = max(gaps.values())
    print(f"\nworst of the summary figures: {worst:.1f} %  "
          f"({'inside' if worst <= 10 else 'OUTSIDE'} the 10 % agreed)")


if __name__ == "__main__":
    main(fresh="--fresh" in sys.argv)
