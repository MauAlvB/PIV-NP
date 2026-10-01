"""What is the ceiling? How well does a PIV measurement of this test agree with itself?

The built-in PIV and PIVlab correlate at about 0.5 on the strain of each particle, and that
number does not improve when particles are averaged into patches of up to 48. I read that as
a *structured* disagreement, on the grounds that point-to-point noise would average away.

That conclusion has an untested assumption in it: that 0.5 is low. Low compared with what?
Nothing here says what two independent measurements of this field *could* agree to. If
PIVlab disagrees with itself by about as much, then 0.5 is the reproducibility floor of the
measurement, there is nothing structured to explain, and the open question in STATUS.md is
not a question.

So this measures the floor, from PIVlab's own data and with no model of the noise. The twenty
steps are split into the odd ones and the even ones, each run as a complete ten-step
analysis. The two sample the same physical event, interleaved, so they should give the same
strain *pattern* at about half the magnitude. How well they agree is how well this
measurement reproduces -- and nothing we build can beat it.

The same split is then applied to the built-in PIV, which answers a second question: is our
field *less* reproducible than PIVlab's, or merely different from it?

    python experiments/piv_reproducibility.py
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
STEPS = 20
#: Steps in each half. Ten of the twenty, taken every other one.
HALF = STEPS // 2


def pearson(x: np.ndarray, y: np.ndarray) -> float:
    """Correlation by explicit sums; this environment crashes on matrix products."""
    keep = np.isfinite(x) & np.isfinite(y)
    if keep.sum() < 3:
        return float("nan")
    dx, dy = x[keep] - x[keep].mean(), y[keep] - y[keep].mean()
    bottom = np.sqrt((dx * dx).sum() * (dy * dy).sum())
    return float((dx * dy).sum() / bottom) if bottom > 0 else float("nan")


def half_case(which: str, source: str) -> Path:
    """A ten-step case built from every other step of the twenty.

    The ``datos`` files are renumbered so the analysis sees a normal sequence, and ``dt`` is
    doubled because each step now spans two frames. For the built-in PIV there are no files
    to renumber: ``FIRST_IMAGE`` and a doubled step would be needed, which the settings
    cannot express, so the photographs are renamed instead.
    """
    work = Path(tempfile.mkdtemp())
    shutil.copytree(CASE, work, dirs_exist_ok=True)
    first = 1 if which == "odd" else 2

    par = work / "dambreak.PAR"
    lines = par.read_text(encoding="latin-1").splitlines()
    for number, line in enumerate(lines):
        if line.strip().startswith("1  20"):
            lines[number] = line.replace("1  20", "2  10", 1)
            break
        if line.strip().startswith("1 20"):
            lines[number] = line.replace("1 20", "2 10", 1)
            break
    par.write_text("\n".join(lines) + "\n", encoding="latin-1")

    if source == "pivlab":
        folder = work / "pivlab"
        wanted = [first + 2 * k for k in range(HALF)]          # ten measurements
        kept = {}
        for new, old in enumerate(wanted, 1):
            kept[new] = (folder / f"datos ({old}).txt").read_bytes()
        for path in folder.glob("datos (*).txt"):
            path.unlink()
        for new, data in kept.items():
            (folder / f"datos ({new}).txt").write_bytes(data)
    else:
        # the built-in PIV reads image n and n+1, so step k of this half must read the
        # photographs of the original steps; renaming into a scratch folder does that
        settings = work / "dambreak.PIV"
        text = settings.read_text(encoding="latin-1")
        source_folder = None
        for line in text.splitlines():
            if line.strip().upper().startswith("IMAGES"):
                source_folder = line.split("=", 1)[1].split("!")[0].strip()
        images = work / "half_images"
        images.mkdir()
        pattern = Path(source_folder)
        # ten steps need eleven photographs, every other one: 1,3..21 or 2,4..22
        for new, old in enumerate([first + 2 * k for k in range(HALF + 1)], 1):
            original = Path(str(pattern).replace("{n:03d}", f"{old:03d}"))
            if not original.exists():
                raise FileNotFoundError(f"{original} is needed for the {which} half")
            shutil.copy(original, images / f"img_{new:03d}.jpg")
        new_pattern = str((images / "img_{n:03d}.jpg").as_posix())
        out = []
        for line in text.splitlines():
            if line.strip().upper().startswith("IMAGES"):
                out.append(f"IMAGES = {new_pattern}")
            else:
                out.append(line)
        settings.write_text("\n".join(out) + "\n", encoding="latin-1")
    return work


def analyse(directory: Path, source: str) -> dict:
    sim = Simulation.from_directory(directory,
                                    options=RunOptions(source=source, prefetch=0))
    sim.run()
    p = sim.particles
    return {
        "eq": np.array(p.eq_strain),
        "shear": np.array(p.strain[:, 2]),
        "ux": np.array(p.displacement[:, 0]),
        "uy": np.array(p.displacement[:, 1]),
        "position": np.array(p.position),
        "shown": output_mask(p, sim.grid, sim.config.total_steps) & (p.nan_initial == 0),
    }


def patches(position: np.ndarray, keep: np.ndarray, side: float) -> np.ndarray:
    x, y = position[keep, 0], position[keep, 1]
    if side <= 0:
        return np.arange(x.size)
    columns = np.floor(x / side).astype(np.int64)
    rows = np.floor(y / side).astype(np.int64)
    return rows * (columns.max() - columns.min() + 1) + columns


def by_patch(groups: np.ndarray, left: np.ndarray, right: np.ndarray):
    labels, index = np.unique(groups, return_inverse=True)
    count = np.bincount(index, minlength=labels.size)
    return (np.bincount(index, weights=left, minlength=labels.size) / count,
            np.bincount(index, weights=right, minlength=labels.size) / count)


def report(name: str, odd: dict, even: dict) -> None:
    moved = (odd["shown"] & even["shown"]
             & (np.hypot(odd["ux"], odd["uy"]) > 5e-4))     # half the steps, half the move
    print(f"\n{name}: {moved.sum()} particles in both halves and moving")
    print(f"{'patch side':>12}{'per patch':>11}{'eq strain r':>14}{'shear r':>10}"
          f"{'displacement r':>16}")
    for side in (0.0, 0.02, 0.04, 0.08):
        groups = patches(odd["position"], moved, side)
        eq = pearson(*by_patch(groups, odd["eq"][moved], even["eq"][moved]))
        sh = pearson(*by_patch(groups, odd["shear"][moved], even["shear"][moved]))
        ux = pearson(*by_patch(groups, odd["ux"][moved], even["ux"][moved]))
        per = moved.sum() / max(len(set(groups.tolist())), 1)
        label = "none" if side == 0 else f"{1000 * side:.0f} mm"
        print(f"{label:>12}{per:>11.1f}{eq:>14.3f}{sh:>10.3f}{ux:>16.3f}")


def main() -> None:
    for source in ("pivlab", "images"):
        print(f"\n{'=' * 72}\nsplitting the twenty steps, measured with --source {source}")
        halves = {}
        for which in ("odd", "even"):
            halves[which] = analyse(half_case(which, source), source)
        report(f"{source} odd half vs its own even half", halves["odd"], halves["even"])

    print("\nThe displacement column is the control: it is well above the noise, so it has")
    print("to correlate near 1. If the strain columns are far below it, the strain of a")
    print("single particle is simply not reproducible on this test, by either method.")


if __name__ == "__main__":
    main()
