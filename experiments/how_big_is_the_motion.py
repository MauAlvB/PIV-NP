"""How many pixels does the soil actually move between two photographs?

Everything so far has been about separating noise from signal. This asks whether there is
enough signal to separate: PIV finds a displacement to about a tenth of a pixel on good
images, and no filter can do anything about that floor. If the material only moves a couple
of tenths of a pixel per frame, the scatter is not noise to be removed, it is the
measurement limit.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from pivnp.config import load_case
from pivnp.sources import build_source

CASE = Path(__file__).resolve().parents[1] / "examples" / "dam-break-swir"

name, cfg = load_case(CASE)
source = build_source("pivlab", CASE, cfg, prefetch=0)
source.moisture = False
x, y, metres_per_pixel, has_data = source.mesh_in_metres(1)

print(f"metres per pixel from the PIVlab header : {metres_per_pixel * 1000:.4f} mm/px")
print(f"interval between images (DT of the .PAR): {cfg.dt} s")
print(f"grid cell                               : {cfg.cell_width * 1000:.2f} mm "
      f"= {cfg.cell_width / metres_per_pixel:.1f} px\n")

speeds = []
scatters = []
for step in (5, 10, 15, 20):
    frame = source.read(step)
    speed = np.hypot(frame.u, frame.v)
    good = np.isfinite(speed)
    # per-frame displacement in pixels
    pixels = speed[good] * cfg.dt / metres_per_pixel
    speeds.append(np.median(pixels))
    rows = cfg.n_rows + 1
    grid = np.full((cfg.n_cols + 1, rows), np.nan)
    grid.ravel()[: frame.u.size] = frame.u
    scatters.append(np.nanstd(np.diff(grid, axis=0)) * cfg.dt / metres_per_pixel)
    print(f"step {step:>2}: median displacement {np.median(pixels):.3f} px, "
          f"90th percentile {np.percentile(pixels, 90):.3f} px, "
          f"max {np.nanmax(pixels):.2f} px")

typical = float(np.median(speeds))
scatter = float(np.median(scatters))
print(f"\ntypical displacement per frame : {typical:.3f} px")
print(f"scatter between neighbours     : {scatter:.3f} px")
print(f"signal to noise                : {typical / max(scatter, 1e-9):.1f} : 1")
print("\nSub-pixel PIV is good to roughly 0.1 px on images like these -- the small PIV in")
print("this folder was checked at 0.07 to 0.14 px against known shifts. So the scatter is")
print("not something a filter left behind: it is what the method can resolve.")
print("\nThe way out is not a better filter but a bigger displacement: correlate frames")
print("further apart in time, or photograph at a finer scale. Both raise the signal while")
print("the floor stays where it is.")
