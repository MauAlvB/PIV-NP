# Visualising the results with ParaView

[ParaView](https://www.paraview.org/download/) is free, open source (BSD) and runs on
Windows, Linux and Mac. It is the standard tool for scientific post-processing and covers
what used to be done with GiD: animation over time, colouring by any result, displacement
arrows, slices, plots along lines and export of images and video.

ParaView does not read GiD's `.POST.RES`, so PIV-NP converts them to VTK (its native
format): one `.vtu` per instant and a `.pvd` that groups them as an animation. It takes
about 6 times less space (768 MB → 130 MB on the centrifuge case) and the data is stored in
single precision, the same as the 6 digits of the GiD text.

![Equivalent shear strain](img/centrifuge_equi_strain.png)

## 1. Generating the VTK files

While computing:

```bash
pivnp path/to/case --vtk
```

From results that already exist:

```bash
pivnp-vtk path/to/case/zapatak.POST.RES
```

For results of the original Fortran executable (or computed with `--legacy-compat`) you
have to add `--legacy-msh`, because its mesh carries the positions of the first step
instead of the initial ones:

```bash
pivnp-vtk path/to/case/zapatak.POST.RES --legacy-msh
```

A folder `zapatak_vtk/` is created with `zapatak.pvd`. With `--every 5` one instant out of
every 5 is exported, and with `-o folder` you choose where they are stored.

## 2. Opening it in ParaView

1. **File → Open** → `zapatak.pvd` → **Apply**.
2. Press the **2D** button of the view (or **Reset Camera** and look from −Z).
3. Under **Representation**, pick **Points** (or **Point Gaussian** for round points) and
   set **Point Size** to 3–5.
4. In the colour drop-down, pick the result: `Equi_strain`, `Displacement` (magnitude or
   the X/Y components), `Vol_strain`…
5. **Play** (▶) on the top bar to watch the animation. The time is GiD's (`Isochrones`).

Tips:

* **Fixed colour scale**: *Edit color map → Rescale to custom range* (for example, from 0
  to 0.5 for `Equi_strain`); otherwise the scale changes at every instant.
* **Displacement arrows**: *Filters → Glyph*, Orientation = `Displacement`,
  Scale Array = `No scale array`, and *Masking → Every Nth Point* so it does not saturate.
* **Profile along a line** (for instance, the displacement along a vertical):
  *Filters → Plot Over Line*.
* **Evolution of one particle over time**: select it with *Select Points On* and apply
  *Filters → Plot Selection Over Time*. The `id` array is GiD's particle number.
* **Filtering particles**: *Filters → Threshold* on any result.
* **Is it shearing or just turning?** Colour by `Vorticity` beside `Equi_strain`. A rigid
  rotation shows up in the first and, misleadingly, in the second too: see the note in the
  [README](../README.md) on telling rotation apart from shear. Use a diverging colour map
  (*Cool to Warm*) centred on zero, since the sign is the direction of the turn.
* **Which regime is each zone in?** Colour by `Rot_angle` with the range fixed to 0–90:
  blue is material deforming without turning, 45 is a shear band, 90 is a block rotating
  rigidly. **Threshold by `Equi_strain` first** (*Filters → Threshold*, keeping the upper
  part): where the material barely moved the field is a ratio of two near-zero numbers and
  shows only noise. Without that filter the map looks like static — that is the single most
  common way to misread it.
* **Saving the setup** to reuse it with other tests: *File → Save State* (`.pvsm`) and,
  when loading it, *Search files under specified directory* with the new folder.
* **Video**: *File → Save Animation* (`.avi`/`.ogv`) or a series of `.png`.

## 3. What each file holds

The coordinates are the **current** positions of the particles at each instant (no need for
*Warp By Vector*). Point arrays:

| Array | Components |
|---|---|
| `id`, `material` | GiD particle number and material (1 active, 2 without data) |
| `Displacement`, `Inst_displacement`, `Velocity`, `Acceleration` | vector (x, y, 0) |
| `Total_strain`, `Inc_strain` | xx, yy, xy (engineering γxy) |
| `Equi_strain`, `In_E_strain`, `Vol_strain`, `Ins_vol_strain` | scalar |
| `Vorticity` | scalar (s⁻¹): the curl, which tells rotation apart from shear |
| `Rotation` | scalar (degrees): rotation accumulated by the particle |
| `Vorticity_num`, `Rot_angle` | scalar: how the deformation splits between shear and rotation |
| `Finite_strain` | xx, yy, xy from the deformation gradient: zero for a rigid rotation |
| `Fin_equi_strain`, `Finite_rotation`, `Finite_area` | scalar: the finite counterparts |
| `E_potential`, `E_kinetic`, `E_total` | scalar |
| `NaNs` | data missing around the particle (0 to 4 with IVERSION=1; 0 or 1 with IVERSION=2) |
| `Moisture`, `Saturation` | scalar (only with `MOISTER=1`) |

## 4. From Python (optional)

With the `viz` extra (`pip install -e ".[viz]"`) the same files can be opened with
[PyVista](https://pyvista.org), which uses the same engine as ParaView:

```python
import pyvista as pv

reader = pv.get_reader("zapatak_vtk/zapatak.pvd")
reader.set_active_time_point(len(reader.time_values) - 1)
grid = reader.read()[0]
grid.plot(scalars="Equi_strain", clim=(0, 0.5), cmap="turbo", point_size=4,
          render_points_as_spheres=True, cpos="xy")
```

![Accumulated displacement](img/centrifuge_displacement.png)

## Other free alternatives

* **VisIt** (LLNL, BSD): similar to ParaView, reads the same `.vtu`/`.pvd`.
* **GiD**: it has a free evaluation version with a limit on the mesh size (depending on the
  version and the licence); the `.POST.RES` files are still compatible.
