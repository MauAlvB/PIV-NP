"""Build the files of the ``shear-block`` example.

The example is already in the repository, so you do not need to run this: it is here so that
you can see exactly what the data is, and change it to try something else. It writes

* ``datos (1).txt`` ... ``datos (20).txt`` -- the velocity field, in the format PIVlab
  exports, as if a PIV analysis had measured it;
* ``images/wet_1.png`` ... ``wet_20.png`` -- photographs of the same test, in which a
  wetting front rises through the block.

The field is a **simple shear**: the horizontal velocity grows linearly with the height
above the base, and the vertical velocity is zero. That has an exact solution, which is what
makes it useful as a first example (see README.md).

Run it with::

    python examples/shear-block/make_example.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from PIL import Image

HERE = Path(__file__).parent

# --- the grid ---------------------------------------------------------------------------
CELL = 0.01          # m, the side of a grid cell
N_COLS, N_ROWS = 12, 8   # cells; the points are one more in each direction
ORIGIN_X, ORIGIN_Y = 0.01, 0.01   # m, where the first point sits on the image

# --- the test ---------------------------------------------------------------------------
STEPS = 20
DT = 0.1             # s between images
SHEAR_RATE = 0.05    # 1/s; du/dh, with h the height above the base

#: Points where PIVlab measured nothing, as happens past the edge of the material. The
#: example ships with none, so that its solution is exact everywhere. Set them to, say,
#: ``range(10, 13)`` and ``range(0, 2)``, regenerate, and you get a case with a hole in the
#: top-right corner: a good way to see what the ``contour`` of the ``.PAR`` is for, because
#: the particles around the hole stop giving the exact answer. Rows count from the top.
NO_DATA_COLUMNS: range = range(0)
NO_DATA_ROWS: range = range(0)

# --- the images -------------------------------------------------------------------------
METRES_PER_PIXEL = 0.0005        # 0.5 mm per pixel, so a cell is 20 px
WIDTH, HEIGHT = 280, 200         # px
DRY_GRAY, SATURATED_GRAY = 150, 100   # the band the .HUM declares
FRONT_BLUR = 10                  # px over which the front goes from dry to saturated


def grid_points() -> tuple[np.ndarray, np.ndarray]:
    """x and y of every point, in metres, in the order PIVlab writes them.

    PIVlab goes by columns: x stays put while y grows **downwards**, which is the image
    axis. PIV-NP flips it round when it reads the file.
    """
    column, row_from_top = np.divmod(np.arange((N_COLS + 1) * (N_ROWS + 1)), N_ROWS + 1)
    return ORIGIN_X + column * CELL, ORIGIN_Y + row_from_top * CELL


def velocity_field() -> tuple[np.ndarray, np.ndarray]:
    """Horizontal and vertical velocity of every point, in m/s, with the gaps as NaN."""
    x, y = grid_points()
    base_y = ORIGIN_Y + N_ROWS * CELL          # the base of the block, lowest in the image
    u = SHEAR_RATE * (base_y - y)              # grows with the height above the base
    v = np.zeros_like(u)

    column, row_from_top = np.divmod(np.arange(x.size), N_ROWS + 1)
    missing = (np.isin(column, list(NO_DATA_COLUMNS))
               & np.isin(row_from_top, list(NO_DATA_ROWS)))
    u[missing] = np.nan
    v[missing] = np.nan
    return u, v


def write_velocity_files() -> None:
    """One file per step, with the three header lines PIVlab writes."""
    x, y = grid_points()
    u, v = velocity_field()
    # PIVlab states both factors; their ratio is the interval between images, which PIV-NP
    # compares with the DT of the .PAR and warns about if they disagree.
    header = (
        "PIVlab by W.Th. & E.J.S., ASCII chart output - synthetic example\n"
        "FRAME: {step}, filenames: A: wet_{step}.png & B: wet_{next}.png, "
        f"conversion factor xy (px -> m): {METRES_PER_PIXEL}, "
        f"conversion factor uv (px/frame -> m/s): {METRES_PER_PIXEL / DT}\n"
        "x [m],y [m],u [m/s],v [m/s]\n"
    )
    for step in range(1, STEPS + 1):
        rows = [
            f"{xi:.6f},{yi:.6f},"
            f"{'NaN' if np.isnan(ui) else format(ui, '.9g')},"
            f"{'NaN' if np.isnan(vi) else format(vi, '.9g')}"
            for xi, yi, ui, vi in zip(x, y, u, v, strict=True)
        ]
        text = header.format(step=step, next=step + 1) + "\n".join(rows) + "\n"
        (HERE / f"datos ({step}).txt").write_text(text, encoding="ascii")


def wet_image(step: int) -> np.ndarray:
    """The block photographed at ``step``, with the wetting front part way up.

    Wet soil looks darker. The front climbs from the base to the top of the block over the
    twenty steps, with a soft edge so that it does not look like a drawn line.
    """
    base_row = (ORIGIN_Y + N_ROWS * CELL) / METRES_PER_PIXEL   # the base, in pixels
    height_risen = N_ROWS * CELL * step / STEPS                # m climbed by the front
    front_row = base_row - height_risen / METRES_PER_PIXEL

    row = np.arange(HEIGHT)[:, None]
    # 0 well above the front (dry) and 1 well below it (saturated)
    wetness = np.clip((row - front_row) / (2 * FRONT_BLUR) + 0.5, 0.0, 1.0)
    gray = DRY_GRAY + wetness * (SATURATED_GRAY - DRY_GRAY)
    return np.broadcast_to(gray, (HEIGHT, WIDTH))


def write_images() -> None:
    folder = HERE / "images"
    folder.mkdir(exist_ok=True)
    for step in range(1, STEPS + 1):
        array = np.rint(wet_image(step)).astype(np.uint8)
        Image.fromarray(array, mode="L").save(folder / f"wet_{step}.png", optimize=True)


def main() -> None:
    write_velocity_files()
    write_images()
    points = (N_COLS + 1) * (N_ROWS + 1)
    print(f"{STEPS} velocity files of {points} points and {STEPS} images "
          f"of {WIDTH}x{HEIGHT} px written in {HERE}")


if __name__ == "__main__":
    main()
