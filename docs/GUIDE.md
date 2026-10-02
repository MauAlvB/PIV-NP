**English** · [`README.md`](../README.md)

# Getting started

This guide takes you from a fresh clone to analysing your own test. No prior knowledge of
PIV-NP is assumed. It is in four parts:

1. [Install it](#1-install-it)
2. [Run the example](#2-run-the-example) — and check that the answer is right
3. [Turn the moisture on](#3-turn-the-moisture-on)
4. [Make your own case](#4-make-your-own-case) — every input file, with and without moisture

If something goes wrong, [What usually goes wrong](#what-usually-goes-wrong) at the end
lists the messages PIV-NP actually prints and what to do about each.

## What PIV-NP does, in one paragraph

You photograph a test — a slope in a centrifuge, a dam breaking, a soil sample being
sheared. A PIV program (today [PIVlab](https://pivlab.blogspot.com/)) compares consecutive
photographs and reports **velocities on a fixed grid of points**. That is useful, but it
does not tell you how far a given lump of soil has travelled or how much it has deformed,
because the grid stays put while the soil moves. PIV-NP fills the grid with **numerical
particles**, moves each one with the velocity around it, step by step, and accumulates its
displacement and strain. It can also read the colour of the soil in the same photographs
and turn it into **water content**. The output opens in GiD or, free of charge, in ParaView.

---

## 1. Install it

You need Python 3.10 or newer. With [conda](https://docs.conda.io/):

```bash
conda env create -f environment.yml
```

```bash
conda activate pivnp
```

```bash
pip install -e .
```

With pip alone, inside a virtual environment:

```bash
pip install -e .
```

Check that it answers:

```bash
pivnp --help
```

The first analysis takes a few seconds longer while Numba compiles the computational
kernels. That is cached, so every later run starts immediately.

---

## 2. Run the example

The repository ships a complete case, data included:

```bash
pivnp examples/shear-block
```

It takes well under a second and prints something like:

```
READING DATA... case shearblock: 384 particles, 20 steps
ANALYZING datos (1).txt
ANALYZING datos (10).txt
ANALYZING datos (20).txt
ANALYSIS FINISHED in 0.4 s
```

Three files appear in `examples/shear-block/`:

| File | What it holds |
|---|---|
| `shearblock.POST.MSH` | the particles at their starting positions |
| `shearblock.POST.RES` | every result, at every printed step |
| `shearblock.REC` | the final state, so a second stage can continue from it |

### Check that the answer is right

This is the part that matters, and the reason the example is a made-up test rather than a
real one. The case is a block of soil in **simple shear**: the horizontal velocity grows
linearly with height, from 0 at the base to 4 mm/s at the top, for 2 seconds. That has an
exact answer, so you can tell the difference between "it did not crash" and "it works".

Paste this into a Python prompt:

```python
import numpy as np
from pivnp import Simulation
from pivnp.solver import output_mask

sim = Simulation.from_directory("examples/shear-block")
sim.run()

p = sim.particles
published = output_mask(p, sim.grid, sim.config.total_steps)
print("shear strain :", np.unique(np.round(p.strain[published, 2], 9)))
print("eq. strain   :", np.unique(np.round(p.eq_strain[published], 6)))
print("displacement : %.2f to %.2f mm" % (p.displacement[published, 0].min() * 1000,
                                          p.displacement[published, 0].max() * 1000))
print("particles    :", published.sum(), "of", p.lost.size)
```

You should get **exactly** this:

```
shear strain : [0.1]
eq. strain   : [0.057735]
displacement : 0.25 to 7.75 mm
particles    : 372 of 384
```

Every particle shows the same shear strain, 0.100, and the same equivalent strain,
0.100/√3 = 0.057735. The normal strains, the volumetric strain and the vertical
displacement are all zero. If you get those numbers, your installation is correct.

The 12 particles that are missing from the 384 are the ones on the right-hand edge: as the
block shears they travel out of the measured grid, and a particle that leaves stops being
computed. That happens in real tests too.

### Look at the results

`.POST.RES` is a GiD text file. To see it in [ParaView](https://www.paraview.org/), which
is free, convert it:

```bash
pivnp examples/shear-block --vtk
```

That writes `shearblock_vtk/shearblock.pvd`; open it with **File → Open**, press **Apply**,
switch the view to **2D**, choose **Points** as the representation and colour by
`Equi_strain`. [`VISUALIZATION.md`](VISUALIZATION.md) covers the rest, including animation
and displacement arrows.

---

## 3. Turn the moisture on

PIV-NP can measure how wet the soil is from the same photographs, because wet soil looks
darker. The switch is a **single number** in block 3 of the `.PAR`.

Open `examples/shear-block/shearblock.PAR` and change `moisture` from `0` to `2`:

```
BLOCK 3: dt  total_steps print_every moisture mesh_version pivlab contour restart track
         0.1 20          1           2        1            1      0       0       0
```

Run it again:

```bash
pivnp examples/shear-block
```

The log now also says where the moisture came from and how much of it is a real
measurement:

```
Moisture from the images: images/wet_{n}.png, channel 0, sigma 3, calibration calibration.csv
Moisture: of 2340 values with data: 2340 (100 %) measured
```

The results carry two extra fields, `Moisture` (water content, %) and `Saturation` (0 to 1).
Everything else — displacements, strains, velocities — comes out **exactly the same**: the
moisture is measured alongside the motion, it does not feed back into it.

In this example a wetting front climbs the block over the 20 steps, so at the end the
saturation runs from 1.00 at the base to about 0.55 at the top.

### The three values of `moisture`

| Value | Where the moisture comes from | Extra files needed |
|---|---|---|
| `0` | nowhere; it is not computed | none |
| `1` | read from `Moist_<n>.TXT`, one per step | the `Moist_<n>.TXT` files |
| `2` | measured from the test photographs | `<case>.HUM`, the calibration CSV, the images |

Use `1` when another program (for example an older MATLAB workflow) already produced the
moisture fields. Use `2` to let PIV-NP measure them.

To look at the moisture fields on their own, without running the analysis:

```bash
pivnp examples/shear-block --write-moisture
```

That writes `Moist_1.TXT` … `Moist_20.TXT` and stops.

### Read that quality line

That last line is worth understanding, because it is easy to over-read a moisture map. Every
node gets one of three marks:

* **measured** — the gray fell inside the calibration table. The value is a measurement.
* **at the wet limit (an 'at least')** — the gray fell at or past the saturated end of the
  table. All you know is that the node is *at least* that wet.
* **at the dry limit (an 'at most')** — the same at the other end.

In the example every value is a measurement, which is how it should be. If a large share of
yours land at a limit, it is usually one of two things: the calibration table does not reach
far enough past the range your test produces, or the band in the `.HUM` is too narrow or
badly placed. The calibration of the example deliberately carries a row beyond each end for
exactly that reason — have a look at the comments in `calibration.csv`.

This matters more than it sounds. On the published test used as a reference, 82 % of the
values with data sat at the wet limit: they were bounds, not measurements.
[`VALIDATION.md`](VALIDATION.md) measures how much each of these decisions moves the answer.

Remember to set `moisture` back to `0` if you want the plain example again.

---

## 4. Make your own case

A case is **one directory**. PIV-NP never looks outside it. File names are
case-insensitive.

### Without moisture: three things

```
my-test/
├── PIV-NP.TXT          one line: the case name, say  mytest
├── mytest.PAR          what to analyse
└── pivlab/             the velocity field of each step
    ├── datos (1).txt
    ├── datos (2).txt
    └── ...
```

The velocity files go in a `pivlab/` subfolder. A test of 149 steps has 149 of them, and
leaving them in the case root buries the two configuration files among them — and the
results land there too once you run it. PIV-NP also reads them straight from the root,
which is how every case was laid out before, so an old case needs no reorganising. What it
will not do is take them from both places at once: that is an error, not a silent choice.

#### `PIV-NP.TXT`

A single line with the case name, no extension. Every other file is named after it. If you
would rather not have this file, pass `pivnp my-test --case mytest` instead.

#### `pivlab/datos (n).txt` — the velocity field

These are what PIVlab exports. One file per step, numbered from 1 with no gaps, with
**three header lines** and then one line per grid point:

```
PIVlab by W.Th. & E.J.S., ASCII chart output - 14-Mar-2017
FRAME: 1, filenames: A: img001.jpg & B: img002.jpg, conversion factor xy (px -> m): 0.0042423, conversion factor uv (px/frame -> m/s): 0.0042423
x [m],y [m],u [m/s],v [m/s]
5.939275568,4.717481737,NaN,NaN
5.939275568,4.929598722,0.000274531,-1.313988693e-05
...
```

Four things have to be right:

* **The point order.** PIVlab writes them by columns: `x` stays fixed while `y` grows
  *downwards* (the image axis), then `x` moves to the next column. PIV-NP flips the vertical
  axis when it reads them. If you produce these files yourself, keep that order.
* **Points with no measurement carry `NaN`** in `u` and `v`. Do not replace them with zeros:
  PIV-NP needs to know the difference between "it is not moving" and "nobody measured it".
* **The second header line must carry the two conversion factors.** `px -> m` places the
  grid on the image, and is required when you use `moisture = 2`. The ratio of the two is
  the time between photographs, which PIV-NP compares with your `DT` and warns about if they
  disagree — a mistake whose only symptom is results at the wrong scale.
* **The number of points must match the `.PAR`.**

Export them from PIVlab as ASCII with the calibration applied, so the header carries those
factors; PIVlab's own documentation covers the menu. If your export has **five** columns
instead of four (newer PIVlab adds the vector type), set `pivlab = 2` in the `.PAR`.

#### `<case>.PAR` — what to analyse

Four blocks. The comment lines carry the field names above the values, and a legend after
block 4 explains every option; neither is read, both are there to be read by you.

```
My test: what it is
BLOCK 2: n_cells n_nodes n_part_cell n_rows width height
         96      117     2           8      0.01  0.01
BLOCK 3: dt  total_steps print_every moisture mesh_version pivlab contour restart track
         0.1 20          1           0        1            1      0       0       0
BLOCK 4: s_density porosity
         2650.0    0.4
```

**Block 2, the grid.** Take these from your PIV settings:

| Field | What it is |
|---|---|
| `n_cells` | cells of the grid (points minus one, in each direction, multiplied) |
| `n_nodes` | grid points; it **must** equal `(n_cells/n_rows + 1) × (n_rows + 1)` |
| `n_part_cell` | particles per cell side, 1 to 6. `2` gives 4 particles per cell, `3` gives 9 |
| `n_rows` | rows of cells |
| `width`, `height` | size of one cell, in metres |

If your PIVlab grid is 13 × 9 points, then it is 12 × 8 cells: `n_rows = 8`,
`n_cells = 96`, `n_nodes = 117`. Getting this wrong is the most common mistake, and PIV-NP
refuses to run rather than produce nonsense.

**Block 3, the analysis.**

| Field | What it is |
|---|---|
| `dt` | seconds between photographs. Must match what you exported from PIVlab with |
| `total_steps` | how many `datos (n).txt` to process |
| `print_every` | results are written at step 1 and then every `print_every` steps |
| `moisture` | `0` none, `1` from `Moist_<n>.TXT`, `2` from the photographs |
| `mesh_version` | `1` velocities sit at the grid points; `2` at the centre of each cell |
| `pivlab` | `1` four-column export, `2` five-column |
| `contour` | boundary correction: `0` none, `1` neighbour average, `2` particle average, `3` extrapolation |
| `restart` | `0` new analysis, `1` continue from `<case>.REC` |
| `track` | must be `0` |

For `contour`, use **`1`**. At the edge of the material PIVlab measures nothing, so those
points would otherwise hold zero velocity and drag the particles near them. Rebuilding them
from the neighbours that do have data was measured on five real case/step combinations and
beat leaving the zero every time; `2` was often worse than no correction at all. The
numbers are in [`VALIDATION.md`](VALIDATION.md), under the boundary correction.

**Block 4, the soil.** `s_density` in kg/m³ and the initial `porosity`. They are read and
validated but not yet used in the computation, so any sensible value will do for now.

The quickest way to start is to copy the `.PAR` of the example and change the numbers.

### With moisture: three more files

```
my-test/
├── PIV-NP.TXT
├── mytest.PAR          with  moisture = 2
├── mytest.HUM          how to measure the moisture
├── calibration.csv     the curve of your soil
├── pivlab/             datos (1..n).txt
└── images/
    └── wet_1.png ...   one photograph per step
```

With `moisture = 1` instead, where the fields come from another program rather than from
the photographs, the `Moist_<n>.TXT` go in a `moisture/` subfolder (or in the root) and
neither the `.HUM`, the calibration nor the images are needed.

#### The photographs

One per step, numbered the same way as the `datos` files, and **the same photographs the
PIV analysis used** (or ones registered to them — see below). PNG or JPEG. Prefer PNG or
high-quality JPEG: compression changes the gray values, and the whole measurement is based
on them. One gray level can shift the saturation of a node appreciably.

#### `calibration.csv` — the curve of your soil

This comes from the laboratory, and it is specific to each soil and each lighting setup.
Wet samples of your soil to known degrees of saturation, photograph them under the light of
the test, and record the gray value of each one:

```
gray,saturation,moisture
0,1.00,25.0
25,0.75,18.0
50,0.45,11.0
75,0.20,5.0
100,0.00,0.0
```

* `gray` is **normalized**: 0 means "as dark as saturated soil", 100 means "as light as dry
  soil", on the band you declare in the `.HUM`. It must increase down the column.
* `saturation` is the degree of saturation, 0 to 1.
* `moisture` is the water content in %.

Between the points PIV-NP interpolates the way MATLAB's `pchip` does, so results stay
comparable with earlier analyses. Outside the table it does **not** extrapolate: it clamps
to the nearest end and tells you how many values were clamped.

#### `<case>.HUM` — how to measure it

Named values, any order, `!` starts a comment. The minimum:

```
IMAGES         = images/wet_{n}.png   ! {n} is replaced by the step number
CHANNEL        = gray                 ! 1 red, 2 green, 3 blue, 0 or "gray"
SIGMA          = 3                    ! averaging radius, in pixels
CALIBRATION    = calibration.csv
DRY_BAND       = 150                  ! gray of dry soil in your photographs
SATURATED_BAND = 100                  ! gray of saturated soil
```

The two decisions that matter most:

**The band — `DRY_BAND` and `SATURATED_BAND`.** These are the two gray values that anchor
the scale. They are the single most sensitive numbers in the whole method: halving the width
of the band shifts the mean saturation by about 0.06, and moving it by two gray levels
changes the fraction of nodes that end up saturated by ten points. Measure them, do not
guess them. Keep the band **wider than 20 gray levels** — PIV-NP warns below that, because
one level then weighs more than 5 % of the scale.

The alternative is to leave both out and give `DRY_REFERENCE = ref.jpg`, a photograph of the
soil while **dry**. Then the dry value is taken pixel by pixel from that photograph and the
band is set by `DRY_OFFSET` and `SATURATED_OFFSET` around it. That cancels the fixed texture
of the soil, which is its advantage, but it only works if the reference really is a dry
state — and it was measured to fall apart when it is not.

**`SIGMA`,** the radius of the blur applied to each photograph before sampling. It is not
cosmetic: changing it between 15 and 80 px moved the saturation as much as getting the band
width wrong by 30-50 %. Choose it from the grain size and the image scale.

Two more worth knowing:

* `INCREMENTAL = 1` makes a node never get drier than it has been. It is a hypothesis about
  the test, not a measurement, so it is **off** by default; with it on you cannot measure a
  drying soil at all.
* `SCALE_X`, `OFFSET_X`, `SCALE_Y`, `OFFSET_Y` and the four projective terms place the PIV
  grid on the moisture photograph. Leave them alone when it is the same camera. They are
  there for a second camera — an infrared one, say — looking from a different angle.

The full list of keys, with every default, is in the [README](../README.md), under
"`<case>.HUM`".

### Running a test in two stages

If a test was filmed at one frame rate and then another, analyse it in two runs. Do the
first normally, then renumber the second batch of `datos` files from 1, set
`restart = 1`, and run again: PIV-NP picks up the particles where the `.REC` left them and
appends to the same results file.

### Where the displacements come from

Nothing ties PIV-NP to PIVlab. The analysis only ever sees one velocity field per step, so
the data can come from a PIV analysis built in, another PIV package, or a simulation.
`pivnp --source <name>` chooses, and `pivlab` is the default. Writing another one is one
class and one registration — the contract, with a worked example, is in
[`src/pivnp/sources.py`](../src/pivnp/sources.py) and in
[`ARCHITECTURE.md`](ARCHITECTURE.md).

---

## 5. A case with no PIVlab at all

Everything above assumes you have PIVlab's `datos (n).txt` files. If you have only the
photographs, PIV-NP can measure the displacements itself. There is a complete example of
this in the repository, so try it before you build your own:

```bash
pivnp examples/piv-from-images --source images
```

It runs in about a second. Its photographs were **generated** with a shear this much and no
more — γxy = 0.15625 — so you can check the answer rather than trust it:
[`examples/piv-from-images/README.md`](../examples/piv-from-images/README.md) lists every
value that has to come out, and `make_example.py` in that folder is the script that put them
there. If the shear comes back near 0.1545 your installation is measuring correctly.

Read the next paragraph before you use this on real data, because it decides whether the
results are any use to you. The built-in PIV
is a **way in**, not a replacement for PIVlab. Checked against PIVlab on a real dam-break
test, it gets the displacements within 8 % and the *amount* of strain within 7 %, but it
agrees only broadly on *where* the strain is. So it will show you where a slope failed and
roughly how hard; if you are going to publish a strain pattern, check it against PIVlab
first. The numbers and the reasoning are in
[`VALIDATION.md`](VALIDATION.md#8-the-built-in-piv).

### The one file you need

A `<case>.PIV` next to the `.PAR` — the only file that is new, and the example's
[`shearphotos.PIV`](../examples/piv-from-images/shearphotos.PIV) is a working one to copy.
Two values have no default and you must supply both:

```
IMAGES = images/shot_{n:03d}.jpg   ! {n} is the image number; {n:03d} makes it 001, 002...
SCALE  = 0.00054081                ! metres per pixel
```

`SCALE` is the one nobody can guess for you. Measure something whose length you know in one
photograph — the width of the box, a ruler in shot — and divide metres by pixels. Everything
the analysis reports is in metres *because of this number*, so if it is wrong every result is
wrong by the same factor and nothing will look amiss.

Then three that decide what the analysis can see:

```
WINDOW  = 16     ! the patch of soil each vector describes, in pixels
OVERLAP = 0.5    ! so there is a vector every 8 px
PASSES  = 2      ! 32 px first to find the movement, then 16 to pin it down
```

**`WINDOW` is the setting that decides whether your result is worth anything**, so it is
worth a paragraph. The grid step is `WINDOW × (1 − OVERLAP)`, which means any step you want
can be reached with a small window and little overlap, or a large window and a lot of it.
Those two are *not* equivalent and nothing will warn you: the first gives a noisy field and
the second a measurement. Pick the window for how far the soil moves between photographs —
the movement should stay under about a quarter of it — and let the overlap follow.

Getting it wrong is quiet. On a slope test analysed here, a 64 px window ran without
complaint and gave a field nearly three times rougher than PIVlab's; since the strain is a
difference between neighbouring vectors, that came out at *twice* the real strain, with the
rotation doubled too. A 200 px window on the same case matched PIVlab. If your strains look
implausibly large, suspect the window before you suspect the soil.

**If your test changes pace, analyse it in two stages.** The same slope moved 7 px per
photograph while it was failing and well under 1 px once it settled, and no single window
serves both: the 200 px window that measures the failure correctly reads only a third of the
quiet movement afterwards, because a pixel is half a per cent of its width. Run the fast part
with a large window, then set `restart = 1` and run the slow part with a small one —
[§4](#running-a-test-in-two-stages) describes the mechanics. One window for a test that
changes regime will be wrong for one half of it.

A window must be even; a power of two is fastest. And two settings that keep the grid off
the background:

```
REGION     = 232, 213, 1656, 845   ! left, top, right, bottom: where the material is
MASK_BELOW = 25                    ! anything darker than this is not soil
```

Without `REGION` the grid covers the whole photograph, including sky and apparatus that will
never move — more points, no more information, and a slower run. `MASK_BELOW` suits
photographs whose background is already blacked out; if instead you have drawn a mask image,
use `MASK_IMAGE = mask.png` and not both.

The full list of keys, with defaults, is in the
[README](../README.md#input-files).

### Making the grid and the `.PAR` agree

The grid of interrogation windows *is* the PIV grid, so `BLOCK 2` of the `.PAR` has to
describe it. You do not have to work it out: run it, and if the two disagree PIV-NP stops and
prints the `BLOCK 2` line your photographs and settings need.

```
the photographs and the .PIV settings give a grid of 177 x 78 points, and the .PAR
describes 100 x 50. Block 2 of the .PAR should read:
    n_cells 13552  n_nodes 13806  n_rows 77  width 0.00432648  height 0.00432648
```

Paste that in and run again.

### If the results look noisy

The field is smoothed before it is used, because the strain is a difference between
neighbouring vectors and sub-pixel scatter dominates it. `SMOOTH = 0.6` is the default and
was chosen by measurement, not taste. Raising it to `1.0` quietens the strain further at the
cost of flattening genuine detail; `SMOOTH = 0` shows you the raw correlation, which is
worth looking at once to see what the smoothing is for.

If the displacements themselves look wrong rather than noisy, the usual causes are a wrong
`SCALE`, a `WINDOW` too small for how far the soil moves between photographs (keep the
movement under about a quarter of the window), or a `DT` in the `.PAR` that does not match
the time between the photographs — nothing in the images says what that interval was, so
PIV-NP cannot check it for you.

---

## What usually goes wrong

These are the messages PIV-NP actually prints.

**`NN=... does not match (NC/NFIL+1)*(NFIL+1)=...`**
Block 2 does not describe a rectangular grid. Count the *points* of your PIV grid in each
direction; cells are one fewer. For 13 × 9 points: `n_rows = 8`, `n_cells = 96`,
`n_nodes = 117`.

**`NC=... is not a multiple of NFIL=...`**
Same thing: `n_cells` divided by `n_rows` has to give the number of columns exactly.

**`block 3 carries N values and the single format has 9`**, or
**`the comment line of block 3 names the fields in another order`**
An old `.PAR`, from a previous version of the program. Convert it once:

```bash
pivnp my-test --convert-par
```

That rewrites every `.PAR` below the directory in the current format and keeps each
original as `<case>.PAR.orig`. It does not change a single value.

**`block 4 is missing, with the soil density and the porosity`**
Same cause, same fix.

**`The .PAR uses DT=0.8 s, but the ... interval between images of 1 s: displacements will
come out multiplied by 0.8`**
Your `dt` disagrees with what the PIVlab export says. One of the two is wrong, and the
results scale with the ratio. Fix the `.PAR` or re-export. PIV-NP does not touch your `dt`,
because changing it silently would change your results.

**`the .PAR says IPIVLAB=1 (four columns), but the file carries 5 per line`**
Set `pivlab = 2` in block 3.

**`FileNotFoundError: ... datos (7).txt`**
A gap in the numbering, or `total_steps` larger than the number of files you have. The files
must run 1, 2, 3 … with nothing missing. The path in the message tells you where it looked.

**`datos (1).txt is both in ... and in its pivlab/ subfolder`**
You have two copies of the same input. Delete the one you do not want: PIV-NP refuses to
guess, because editing a file that turns out not to be the one being read costs hours.

**`MOISTER=... must be 0 ..., 1 ... or 2 ...`**
Only those three values. Note that in some pre-2024 cases `2` meant "read the files", which
is `1` today; `--convert-par` corrects that for you and says so.

**`the moisture settings file does not exist`**
`moisture = 2` needs `<case>.HUM` beside the `.PAR`, named after the case.

**`CALIBRATION is missing, the file with the curve of the soil`**
Add `CALIBRATION = ...` to the `.HUM`. There is no default: the curve belongs to your soil.

**`the gray column must increase, with no repeats`**
Sort your calibration CSV by the `gray` column and remove duplicates.

**`N grid nodes fall outside the moisture image: check the registration`**
The grid does not land on the photograph. Check that the `px -> m` factor in the `datos`
header is the one for *these* images, and the `SCALE`/`OFFSET` keys of the `.HUM`.

**`the band between dry and saturated soil is 11 gray levels, so one level is 9 % of the
saturation scale`**
Fewer than 20 gray levels between `DRY_BAND` and `SATURATED_BAND`. It still runs, but every
gray level then moves the result a lot, and image noise alone can shift a node's saturation
by a tenth. Improve the lighting or the contrast of the test if you can, and treat the
numbers with care.

**`ITR: PTV tracking particles are no longer supported`**
Set `track = 0` and delete block 5 of the `.PAR`.

---

## Where to go next

* [`README.md`](../README.md) — the reference: every file, every result, every option.
* [`VISUALIZATION.md`](VISUALIZATION.md) — ParaView, step by step.
* [`VALIDATION.md`](VALIDATION.md) — what was checked against earlier versions and against
  a published test, how much each moisture decision moves the answer, and the limits of the
  method that are still open. Worth reading before you trust a moisture map.
* [`ARCHITECTURE.md`](ARCHITECTURE.md) — how the code is laid out, and how to extend it.
