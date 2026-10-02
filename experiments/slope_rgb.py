"""The Slope_RGB test analysed from its own photographs, and against PIVlab.

A different regime from the dam-break, and that is why it is worth running: this slope moves
up to 19.5 px between photographs where the dam break moves 0.16. Everything measured about
the built-in PIV so far was measured where the displacement is a fraction of a pixel, and
nothing said whether that carries over to a test that actually moves.

    python experiments/slope_rgb.py [--fresh] [--steps N]

The case is named by ``SLOPE_RGB_CASE`` and is **never written to**: it holds results from
earlier work. Everything runs on a copy.

Setting it up took three things the case does not carry:

* a ``.PIV``. The scale comes from PIVlab's own export header, 0.0042423 m/px, and the
  window and region are chosen to land on the grid the ``.PAR`` describes: 60 x 35 points,
  50 px apart.
* ``moisture = 0``. The case is set to read ``Moist_<n>.TXT``, which the photograph source
  does not do, and the question here is about displacement and strain.
* bands in the correlation, since a large window over a 4-megapixel photograph asks for
  gigabytes of correlation planes at once.

**The window is what makes or breaks this case, and the first attempt got it wrong.** A
64 px window also steps by 50, at an overlap of 0.219, and it looks reasonable until the
answer is compared: it gives a field nearly three times rougher than PIVlab's, and the strain
is a difference between neighbouring vectors, so that came out at twice PIVlab's, with the
rotation doubled too. Measured over four steps, against PIVlab on the same grid:

=======  ========  ============  ===========
window   measured  our roughness  PIVlab's
=======  ========  ============  ===========
64       64 %      0.257          0.096
100      74 %      0.162          0.096
128      80 %      0.141          0.096
160      86 %      0.115          0.096
**200**  **89 %**  **0.091**      0.096
=======  ========  ============  ===========

PIVlab reaching a 50 px step with its usual half overlap would have used a 100 px window;
200 at 0.75 overlap is what actually matches what it produced. At that window one pass is
enough -- a second changes the roughness from 0.096 to 0.095 and costs 5.7 times the time --
because 19.5 px of movement is a tenth of the window rather than a third of it.
"""
from __future__ import annotations

import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))

from test_images import SLOPE_RGB, case_folder  # noqa: E402

from pivnp.simulation import RunOptions, Simulation  # noqa: E402
from pivnp.solver import output_mask  # noqa: E402

HERE = Path(__file__).parent
NAME = "Slope_Vis"

SCALE = 0.0042423          # m/px, from the PIVlab export header
WINDOW = 200               # see the note below: this is what makes or breaks the case
OVERLAP = 1 - 50 / 200     # so the grid steps by 50 px, which is what PIVlab used
PASSES = 1                 # at this window a second pass changes nothing and costs 5.7x
REGION = (200 - (WINDOW - 1) // 2, 112 - (WINDOW - 1) // 2,
          200 - (WINDOW - 1) // 2 + 59 * 50 + WINDOW + 4,
          112 - (WINDOW - 1) // 2 + 34 * 50 + WINDOW + 4)

PIV_SETTINGS = f"""! Built-in PIV for the Slope_RGB test, written by experiments/slope_rgb.py
IMAGES  = vis_{{n}}.jpg
SCALE   = {SCALE}
CHANNEL = gray
WINDOW  = {WINDOW}
OVERLAP = {OVERLAP:.6f}
PASSES  = {PASSES}
REGION  = {REGION[0]}, {REGION[1]}, {REGION[2]}, {REGION[3]}
"""

RUNS = {
    "PIVlab": dict(source="pivlab", subpixel=False),
    "PIV-NP from the photographs": dict(source="images", subpixel=False),
}
KEYS = ("x", "y", "ux", "uy", "eq", "shear", "rotation", "vorticity", "shown")


def prepare(steps: int, subpixel: bool) -> Path:
    """A working copy of the case, set up for the run. The original is never touched."""
    case = case_folder(SLOPE_RGB, "the Slope_RGB case and its photographs")
    work = Path(tempfile.mkdtemp(prefix="slope_"))
    shutil.copytree(case, work, dirs_exist_ok=True)
    for stale in list(work.glob("*.POST.*")) + list(work.glob("*.REC")):
        stale.unlink()
    (work / "PIV-NP.TXT").write_text(f"{NAME}\n", encoding="ascii")

    settings = PIV_SETTINGS
    if subpixel:
        settings += "SUBPIXEL_OFFSET = 1\n"
    (work / f"{NAME}.PIV").write_text(settings, encoding="ascii")

    # moisture off, and however many steps were asked for
    par = work / f"{NAME}.PAR"
    lines = par.read_text(encoding="latin-1").splitlines()
    for number, line in enumerate(lines):
        if line.strip().startswith("1  149") or line.strip().startswith("1 149"):
            parts = line.split()
            parts[1] = str(steps)
            parts[3] = "0"                        # moisture
            lines[number] = "         " + " ".join(parts)
            break
    par.write_text("\n".join(lines) + "\n", encoding="latin-1")
    return work


def collect(source: str, subpixel: bool, steps: int) -> dict:
    work = prepare(steps, subpixel)
    sim = Simulation.from_directory(work, options=RunOptions(source=source, prefetch=0))
    sim.run()
    p = sim.particles
    out = {
        "x": p.position[:, 0].copy(),
        "y": p.position[:, 1].copy(),
        "ux": p.displacement[:, 0].copy(),
        "uy": p.displacement[:, 1].copy(),
        "eq": np.asarray(p.eq_strain).copy(),
        "shear": p.strain[:, 2].copy(),
        "shown": output_mask(p, sim.grid, sim.config.total_steps) & (p.nan_initial == 0),
    }
    for field in ("rotation", "vorticity"):
        value = getattr(p, field, None)
        out[field] = (np.asarray(value).copy() if value is not None
                      else np.zeros(out["ux"].shape))
    shutil.rmtree(work, ignore_errors=True)
    return out


def run_one(number: int, steps: int) -> None:
    name = list(RUNS)[number]
    data = collect(steps=steps, **RUNS[name])
    np.savez_compressed(HERE / f"_slope_part{number}.npz", **{k: data[k] for k in KEYS})


def gather(fresh: bool, steps: int) -> dict[str, dict]:
    # ask for the case here, before any subprocess is started: otherwise the message about
    # where the data should be comes out of a child and is followed by the parent reporting
    # that the child failed, which buries the only line worth reading
    case_folder(SLOPE_RGB, "the Slope_RGB case and its photographs")
    out = {}
    for number, name in enumerate(RUNS):
        part = HERE / f"_slope_part{number}.npz"
        if fresh or not part.exists():
            print(f"running {name} over {steps} steps...", flush=True)
            done = subprocess.run([sys.executable, str(Path(__file__).resolve()),
                                   "--one", str(number), "--steps", str(steps)],
                                  check=False)
            if done.returncode != 0 or not part.exists():
                raise SystemExit(f"{name} did not finish (exit {done.returncode})")
        stored = np.load(part)
        out[name] = {k: stored[k] for k in KEYS}
    return out


def main(fresh: bool, steps: int) -> None:
    results = gather(fresh, steps)

    shown = np.ones_like(results["PIVlab"]["shown"])
    for data in results.values():
        shown = shown & data["shown"]
    print(f"\nparticles kept by all three: {shown.sum()} of {shown.size}")

    print(f"\n{'':<24}{'|u| [mm]':>11}{'equivalent':>12}{'|shear|':>10}{'rotation':>11}")
    for name, data in results.items():
        magnitude = np.hypot(data["ux"][shown], data["uy"][shown]) * 1000.0
        print(f"{name:<24}{np.nanmedian(magnitude):>11.2f}"
              f"{np.nanmedian(data['eq'][shown]):>12.5f}"
              f"{np.nanmedian(np.abs(data['shear'][shown])):>10.5f}"
              f"{np.nanmedian(np.abs(data['rotation'][shown])):>11.4f}")

    import compare_fields as draw

    for data in results.values():
        data["magnitude"] = np.hypot(data["ux"], data["uy"])

    draw.HERE = HERE
    common = dict(title="Slope_RGB", steps=steps)
    draw.figure(results, "magnitude", 1000.0, "displacement [mm]",
                "slope_displacement.png", **common)
    draw.figure(results, "eq", 1.0, "equivalent strain [-]", "slope_strain.png", **common)
    # signed: a rotation that turns one way over the scarp and the other at the toe is the
    # thing to see, and the size alone would hide it
    draw.figure(results, "rotation", 1.0, "rotation [deg]", "slope_rotation.png",
                signed=True, **common)


if __name__ == "__main__":
    steps = 149
    if "--steps" in sys.argv:
        steps = int(sys.argv[sys.argv.index("--steps") + 1])
    if "--one" in sys.argv:
        run_one(int(sys.argv[sys.argv.index("--one") + 1]), steps)
    else:
        main("--fresh" in sys.argv, steps)
