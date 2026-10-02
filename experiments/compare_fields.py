"""Pictures of the same test analysed three ways, because a table hides *where* things differ.

The dam-break case is run from PIVlab's files, from the photographs as PIV-NP ships, and
from the photographs reading each window between pixels. What comes out is the accumulated
displacement and the accumulated strain of every particle, which is what a user looks at.

    python experiments/compare_fields.py [--fresh]

Writes PNG files next to itself. Needs ``PIVNP_TEST_IMAGES``; the photographs are not in the
repository.

Two things about how this is built are forced by the machine rather than chosen, and both are
noted at the end of ``results.md``:

* each analysis runs in a **process of its own**. Three in one interpreter brings it down
  with ``0xc06d007f``, and freeing the arrays between them is not enough.
* the figures are drawn with **numpy and PIL**, not matplotlib. Matplotlib imports and even
  plots here, then crashes the interpreter inside ``savefig`` -- PNG and SVG alike, so it is
  the font and raster layer, not the format. PIL writes PNGs perfectly well, so the few
  things needed (a colour ramp, a colour bar, a label) are done by hand below.
"""
from __future__ import annotations

import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont

sys.path.insert(0, str(Path(__file__).parent))

from test_images import DAM_BREAK, case_folder  # noqa: E402

from pivnp.simulation import RunOptions, Simulation  # noqa: E402
from pivnp.solver import output_mask  # noqa: E402

CASE = case_folder(DAM_BREAK, "the dam-break case and its photographs",
                   Path(__file__).resolve().parents[1] / "examples" / "dam-break-swir")
HERE = Path(__file__).parent

RUNS = {
    "PIVlab": dict(source="pivlab", extra=None),
    "PIV-NP as shipped": dict(source="images", extra=None),
    "PIV-NP between pixels": dict(source="images", extra="\nSUBPIXEL_OFFSET = 1\n"),
}
KEYS = ("x", "y", "ux", "uy", "eq", "shear", "shown")

#: Cells across a panel, and how much each is blown up to be visible. The grid is 177 x 78
#: points, so a raster near that count gives about one cell per measurement: finer leaves
#: holes between particles and reads as noise that is not there, coarser throws measurements
#: away. The blow-up is nearest-neighbour on purpose -- nothing here is interpolated that was
#: not measured.
RASTER_WIDTH = 180
MAGNIFY = 3


# --- running the analyses ---------------------------------------------------------------
def collect(source: str, extra: str | None) -> dict:
    work = Path(tempfile.mkdtemp())
    shutil.copytree(CASE, work, dirs_exist_ok=True)
    if extra:
        settings = work / "dambreak.PIV"
        settings.write_text(settings.read_text(encoding="latin-1") + extra,
                            encoding="latin-1")
    sim = Simulation.from_directory(work, options=RunOptions(source=source, prefetch=0))
    sim.run()
    p = sim.particles
    return {
        "x": p.position[:, 0].copy(),
        "y": p.position[:, 1].copy(),
        "ux": p.displacement[:, 0].copy(),
        "uy": p.displacement[:, 1].copy(),
        "eq": np.asarray(p.eq_strain).copy(),
        "shear": p.strain[:, 2].copy(),
        "shown": output_mask(p, sim.grid, sim.config.total_steps) & (p.nan_initial == 0),
    }


def run_one(number: int) -> None:
    """Run a single analysis and leave it on disk. This is what each subprocess does."""
    name = list(RUNS)[number]
    data = collect(**RUNS[name])
    np.savez_compressed(HERE / f"_compare_part{number}.npz",
                        **{k: data[k] for k in KEYS})


def gather(fresh: bool) -> dict[str, dict]:
    out = {}
    for number, name in enumerate(RUNS):
        part = HERE / f"_compare_part{number}.npz"
        if fresh or not part.exists():
            print(f"running {name}...")
            done = subprocess.run([sys.executable, str(Path(__file__).resolve()),
                                   "--one", str(number)], check=False)
            if done.returncode != 0 or not part.exists():
                raise SystemExit(f"the analysis of {name} did not finish "
                                 f"(exit {done.returncode})")
        stored = np.load(part)
        out[name] = {k: stored[k] for k in KEYS}
    return out


# --- drawing ----------------------------------------------------------------------------
def ramp(anchors: list[tuple[float, float, float]]) -> np.ndarray:
    """A 256-colour lookup table interpolated between a few anchor colours."""
    anchor = np.array(anchors, dtype=np.float64)
    position = np.linspace(0.0, 1.0, anchor.shape[0])
    wanted = np.linspace(0.0, 1.0, 256)
    table = np.empty((256, 3))
    for channel in range(3):
        table[:, channel] = np.interp(wanted, position, anchor[:, channel])
    return np.clip(table, 0, 255).astype(np.uint8)


#: A perceptual dark-to-bright ramp, for a quantity that goes from none to a lot.
SEQUENTIAL = ramp([(68, 1, 84), (59, 82, 139), (33, 145, 140), (94, 201, 98),
                   (253, 231, 37)])
#: Blue-white-red, for a difference that has a sign and should read as nothing in the middle.
DIVERGING = ramp([(5, 48, 97), (67, 147, 195), (247, 247, 247), (214, 96, 77),
                  (103, 0, 31)])
#: Pale to dark, for a scatter: one particle in a cell should be faint and a pile of them
#: solid. A bright ramp cannot do that -- with six thousand particles over a grid most cells
#: hold exactly one, every cell saturates, and the cloud comes out a uniform slab.
DENSITY = ramp([(206, 220, 236), (107, 152, 201), (33, 76, 138), (10, 20, 50)])
EMPTY = (248, 248, 248)        # where no particle landed


def raster(x: np.ndarray, y: np.ndarray, values: np.ndarray,
           extent: tuple[float, float, float, float], width: int) -> np.ndarray:
    """Average the particles into a picture, ``NaN`` where none landed in a pixel."""
    left, right, bottom, top = extent
    height = max(1, int(round(width * (top - bottom) / (right - left))))
    column = np.clip(((x - left) / (right - left) * (width - 1)).astype(int), 0, width - 1)
    # y upwards in the model, downwards in a picture
    row = np.clip(((top - y) / (top - bottom) * (height - 1)).astype(int), 0, height - 1)
    flat = row * width + column

    good = np.isfinite(values)
    total = np.bincount(flat[good], weights=values[good], minlength=width * height)
    count = np.bincount(flat[good], minlength=width * height)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = np.where(count > 0, total / np.maximum(count, 1), np.nan)
    return mean.reshape(height, width)


def colour(picture: np.ndarray, low: float, high: float, table: np.ndarray) -> np.ndarray:
    """Turn a picture of numbers into one of colours, grey where there is no number."""
    scaled = (picture - low) / (high - low if high > low else 1.0)
    index = np.clip(np.nan_to_num(scaled, nan=0.0) * 255.0, 0, 255).astype(np.uint8)
    out = table[index]
    return np.where(np.isfinite(picture)[..., None], out, np.array(EMPTY, dtype=np.uint8))


def colour_bar(low: float, high: float, table: np.ndarray, width: int,
               height: int = 14) -> np.ndarray:
    strip = np.linspace(0.0, 1.0, width)
    index = (strip * 255).astype(np.uint8)
    return np.repeat(table[index][None, :, :], height, axis=0)


def panel(picture: np.ndarray, low: float, high: float, table: np.ndarray, title: str,
          note: str, unit: str) -> Image.Image:
    """One map, titled, with its scale under it."""
    small = Image.fromarray(colour(picture, low, high, table), mode="RGB")
    body = small.resize((small.width * MAGNIFY, small.height * MAGNIFY), Image.NEAREST)
    width = body.width
    out = Image.new("RGB", (width + 16, body.height + 76), (255, 255, 255))
    out.paste(body, (8, 34))
    out.paste(Image.fromarray(colour_bar(low, high, table, width), mode="RGB"),
              (8, body.height + 42))

    draw = ImageDraw.Draw(out)
    font = ImageFont.load_default()
    draw.text((8, 4), title, fill=(0, 0, 0), font=font)
    draw.text((8, 18), note, fill=(90, 90, 90), font=font)
    draw.text((8, body.height + 60), f"{low:.3g}", fill=(0, 0, 0), font=font)
    draw.text((width // 2 - 20, body.height + 60), unit, fill=(0, 0, 0), font=font)
    right = f"{high:.3g}"
    draw.text((width - 6 * len(right), body.height + 60), right, fill=(0, 0, 0), font=font)
    return out


def stack(panels: list[Image.Image], heading: str) -> Image.Image:
    """Put the panels one under another with a heading."""
    width = max(p.width for p in panels)
    height = sum(p.height for p in panels)
    out = Image.new("RGB", (width, height + 26), (255, 255, 255))
    ImageDraw.Draw(out).text((8, 8), heading, fill=(0, 0, 0),
                             font=ImageFont.load_default())
    offset = 26
    for piece in panels:
        out.paste(piece, (0, offset))
        offset += piece.height
    return out


def figure(results: dict[str, dict], quantity: str, unit: float, label: str,
           filename: str, signed: bool = False, title: str = "Dam break",
           steps: int = 20) -> Path:
    """The three analyses of one quantity, and the two differences from PIVlab.

    ``signed`` is for a quantity that turns both ways, such as a rotation: it gets the
    blue-white-red ramp centred on zero, because showing the size of a rotation and throwing
    away which way it went loses the thing worth seeing.
    """
    shown = np.ones_like(results["PIVlab"]["shown"])
    for data in results.values():
        shown = shown & data["shown"]

    names = list(results)
    x, y = results["PIVlab"]["x"][shown], results["PIVlab"]["y"][shown]
    extent = (x.min(), x.max(), y.min(), y.max())
    values = {n: results[n][quantity][shown] * unit for n in names}

    # The top of the scale is a percentile, not the largest value: a handful of particles on
    # the face of the slope are several times the rest, and letting them set the scale paints
    # everything else black and hides the comparison the figure is for.
    if signed:
        high = float(np.nanpercentile(np.abs(values["PIVlab"]), 97))
        low, table, capped = -high, DIVERGING, "97th percentile of the size"
    else:
        high = float(np.nanpercentile(values["PIVlab"], 95))
        low, table, capped = 0.0, SEQUENTIAL, "95th percentile"
    panels = []
    for name in names:
        picture = raster(x, y, values[name], extent, RASTER_WIDTH)
        panels.append(panel(picture, low, high, table, name,
                            f"median {np.nanmedian(values[name]):.4g}, "
                            f"scale capped at the {capped}", label))

    reference = values[names[0]]
    span = max(float(np.nanpercentile(np.abs(values[n] - reference), 90))
               for n in names[1:]) or 1.0
    scale_of_it = float(np.nanmedian(np.abs(reference)))
    for name in names[1:]:
        gap = values[name] - reference
        picture = raster(x, y, gap, extent, RASTER_WIDTH)
        typical = np.nanmedian(np.abs(gap))
        # as a fraction of the typical *size*, not of the median, which for a signed
        # quantity sits near zero and would turn any difference into a huge percentage
        share = (f" ({100 * typical / scale_of_it:.0f} % of its typical size)"
                 if scale_of_it > 0 else "")
        panels.append(panel(picture, -span, span, DIVERGING,
                            f"{name} minus {names[0]}",
                            f"median difference {typical:.4g}{share}, "
                            f"blue = lower than PIVlab", label))

    out = stack(panels, f"{title}: {label} accumulated over {steps} steps")
    path = HERE / filename
    out.save(path)
    print(f"wrote {path}  ({out.width}x{out.height})")
    return path


def agreement(results: dict[str, dict]) -> Path:
    """Particle against particle: whether the two agree, or only look alike."""
    shown = np.ones_like(results["PIVlab"]["shown"])
    for data in results.values():
        shown = shown & data["shown"]
    reference = results["PIVlab"]
    moving = shown & (np.hypot(reference["ux"], reference["uy"]) > 1e-3)

    # Coarse on purpose. Six thousand particles spread over a fine grid leave almost every
    # cell holding one of them, every cell then has the same density, and the picture says
    # nothing about where the cloud is dense. Cells a few times larger is what makes a cloud
    # look like a cloud.
    size = 150
    panels = []
    for quantity, unit, label in (("ux", 1000.0, "displacement x [mm]"),
                                  ("eq", 1.0, "equivalent strain [-]")):
        base = reference[quantity][moving] * unit
        limit = float(np.nanpercentile(np.abs(base), 99))
        for name in list(results)[1:]:
            mine = results[name][quantity][moving] * unit
            low = 0.0 if quantity == "eq" else -limit
            cloud = raster(np.clip(base, low, limit), np.clip(mine, low, limit),
                           np.ones_like(base), (low, limit, low, limit), size)
            # how many particles fell in each pixel, on a log scale so a thin tail shows
            count = np.nan_to_num(cloud, nan=0.0)
            picture = np.where(count > 0, count, np.nan)
            # one particle is the palest colour and four or more the darkest, which is what
            # makes a thin tail visible without the dense middle swamping it
            image = Image.fromarray(colour(picture, 1.0, 4.0, DENSITY), mode="RGB")
            image = image.resize((size * 2, size * 2), Image.NEAREST)
            draw = ImageDraw.Draw(image)
            draw.line([(0, image.height - 1), (image.width - 1, 0)], fill=(200, 60, 60),
                      width=2)
            frame = Image.new("RGB", (image.width + 16, image.height + 60),
                              (255, 255, 255))
            frame.paste(image, (8, 34))
            text = ImageDraw.Draw(frame)
            font = ImageFont.load_default()
            text.text((8, 4), f"{name}: {label}", fill=(0, 0, 0), font=font)
            agree = 100.0 * np.mean(np.abs(mine - base)
                                    <= 0.1 * np.maximum(np.abs(base), 1e-30))
            text.text((8, 18), f"within 10 % of PIVlab: {agree:.0f} % of particles",
                      fill=(90, 90, 90), font=font)
            text.text((8, image.height + 40), f"PIVlab {low:.3g} to {limit:.3g}  "
                                      "(red line = perfect agreement)",
                      fill=(0, 0, 0), font=font)
            panels.append(frame)

    row_width = 2 * panels[0].width
    out = Image.new("RGB", (row_width, 26 + 2 * panels[0].height), (255, 255, 255))
    ImageDraw.Draw(out).text((8, 8), "Every particle against what PIVlab said about it",
                             fill=(0, 0, 0), font=ImageFont.load_default())
    for number, piece in enumerate(panels):
        out.paste(piece, ((number % 2) * piece.width,
                          26 + (number // 2) * piece.height))
    path = HERE / "compare_agreement.png"
    out.save(path)
    print(f"wrote {path}  ({out.width}x{out.height})")
    return path


def main(fresh: bool = False) -> None:
    results = gather(fresh)

    shown = np.ones_like(results["PIVlab"]["shown"])
    for data in results.values():
        shown = shown & data["shown"]
    print(f"\nparticles kept by all three: {shown.sum()} of {shown.size}")

    print(f"\n{'':<24}{'|u| [mm]':>10}{'equivalent':>12}{'|shear|':>10}")
    for name, data in results.items():
        magnitude = np.hypot(data["ux"][shown], data["uy"][shown]) * 1000.0
        print(f"{name:<24}{np.nanmedian(magnitude):>10.3f}"
              f"{np.nanmedian(data['eq'][shown]):>12.5f}"
              f"{np.nanmedian(np.abs(data['shear'][shown])):>10.5f}")

    for data in results.values():
        data["magnitude"] = np.hypot(data["ux"], data["uy"])

    figure(results, "magnitude", 1000.0, "displacement [mm]", "compare_displacement.png")
    figure(results, "eq", 1.0, "equivalent strain [-]", "compare_strain.png")
    agreement(results)


if __name__ == "__main__":
    if "--one" in sys.argv:
        run_one(int(sys.argv[sys.argv.index("--one") + 1]))
    else:
        main(fresh="--fresh" in sys.argv)
