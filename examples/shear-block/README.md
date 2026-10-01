# Example: a block in simple shear

A small, complete case that runs straight from the repository, with no data to download.
Start here: [`docs/GUIDE.md`](../../docs/GUIDE.md) walks through it step by step.

```bash
pivnp examples/shear-block
```

## What it is

A block of soil 12 × 8 cm, measured on a grid of 12 × 8 cells of 1 cm, over 20 steps of
0.1 s. The velocity field is **simple shear**: the horizontal velocity grows linearly with
the height above the base, from 0 at the base to 4 mm/s at the top, and nothing moves
vertically.

It is not a real test. It was made up precisely because simple shear has an exact answer,
so after running it you can tell whether your installation is working — not just that it
did not crash. [`make_example.py`](make_example.py) is the script that builds the data, so
you can see exactly what it is and change it.

## What has to come out

After 2 s of shearing at 0.05 s⁻¹:

| Result | Expected value |
|---|---|
| `Total_strain` xy (γxy) | **0.100** on every particle |
| `Total_strain` xx and yy | 0 |
| `Vol_strain` | 0 (shear does not change the volume) |
| `Equi_strain` | **0.0577** = 0.100/√3 on every particle |
| `Displacement` x | from **0.25 mm** at the base to **7.75 mm** at the top |
| `Displacement` y | 0 |
| `Vorticity` | **−0.0500** s⁻¹ on every particle |
| `Rotation` | **−2.8648°** on every particle |
| `Vorticity_num` | **1.0000** — the textbook value for simple shear |
| `Rot_angle` | **45.00°** — the same thing on the bounded scale |
| particles published | **372** of 384 |

The vorticity is not zero, which surprises people: simple shear is pure shear **plus** a
rotation, so the block turns while it deforms. Here ω = ∂v/∂x − ∂u/∂y = 0 − 0.05 = −0.05
s⁻¹, clockwise, and over the 2 s that is ω/2 × 2 s = −0.05 rad = −2.8648°.

The 12 missing particles are the ones on the right that leave the grid as the block shears;
once a particle leaves it stops being computed, which is what happens in a real test too.

## The files

| File | What it is |
|---|---|
| `PIV-NP.TXT` | the case name, `shearblock` |
| `shearblock.PAR` | what to analyse: the grid, the 20 steps, and the options |
| `shearblock.HUM` | moisture settings, only read when `moisture = 2` |
| `calibration.csv` | the curve of this soil: gray → saturation and water content |
| `pivlab/datos (1..20).txt` | the velocity field of each step, in the format PIVlab exports |
| `images/wet_1..20.png` | photographs of the test, with a wetting front climbing the block |
| `make_example.py` | regenerates the data files and the images |

The inputs are in subfolders so that the case root holds only the configuration and, after
running it, the results. PIV-NP reads them from the root as well, which is how older cases
are laid out.

## Turning the moisture on

Edit `shearblock.PAR` and change `moisture` in block 3 from `0` to `2`:

```
BLOCK 3: dt  total_steps print_every moisture mesh_version pivlab contour restart track
         0.1 20          1           2        1            1      0       0       0
```

Run it again and the results carry two more fields, `Moisture` and `Saturation`, read from
the images through `shearblock.HUM`. Everything else — displacements, strains, velocities —
comes out **exactly the same**: the moisture is measured alongside, it does not feed back
into the motion.

At the last step the wetting front has climbed the whole block, so the saturation goes from
**1.00** at the base to about **0.55** at the top, and the water content from 25 % to
around 13 %. The log ends with

```
Moisture: of 2340 values with data: 2340 (100 %) measured
```

which is the point of the calibration in `calibration.csv` carrying one row past each end
of the range: a node whose gray lands exactly on the dry or the saturated reference still
counts as a measurement instead of a bound. See
[`docs/GUIDE.md`](../../docs/GUIDE.md#read-that-quality-line).
