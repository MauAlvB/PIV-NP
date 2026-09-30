**English** · [Español](README.es.md)

# PIV-NP: Particle Image Velocimetry with numerical particles

PIV-NP computes **displacements, strains, velocities, accelerations, energies and water
content** of numerical particles that move with the velocity field measured by
[PIVlab](https://pivlab.blogspot.com/) on a sequence of images of a test (for example, a
slope in a geotechnical centrifuge). Unlike classic PIV, which reports velocities on a
fixed (Eulerian) grid, PIV-NP follows each material point (Lagrangian), which makes it
suitable for **large displacements** and accumulated strains.

> Pinyol, N.M. & Alvarado, M. (2017). *Novel PIV-based analysis for large displacement*.
> Canadian Geotechnical Journal 54(7): 933-944.

Version 2.0 is written in Python and is about **19 times faster** than the original Fortran
version (`legacy/`): 38 s → 2 s for the centrifuge case.

![Equivalent shear strain in the centrifuge test](docs/img/centrifuga_equi_strain.png)

---

## Contents

1. [Installation](#installation)
2. [Usage](#usage)
3. [Input files](#input-files)
4. [Results](#results)
5. [How it works](#how-it-works)
6. [Visualisation](#visualisation)
7. [Boundary correction](#boundary-correction)
8. [Performance](#performance)
9. [Tests](#tests)
10. [Project layout](#project-layout)

## Installation

With [conda](https://docs.conda.io/) (recommended):

```bash
conda env create -f environment.yml
```

```bash
conda activate pivnp
```

```bash
pip install -e .
```

With pip only (Python ≥ 3.10):

```bash
pip install -e .
```

The first run takes a few extra seconds while Numba compiles the computational kernels.
The compilation is cached, so later runs start immediately.

## Usage

### Command line

The case directory must contain `PIV-NP.TXT` (with the case name), `<case>.PAR` and the
files exported by PIVlab:

```bash
pivnp path/to/case
```

Options:

| Option | Description |
|---|---|
| `--case NAME` | case name (otherwise read from `PIV-NP.TXT`) |
| `--threads N` | computation threads (default: all cores) |
| `--prefetch N` | PIVlab files read ahead (default 4) |
| `--contour-min-neighbors N` | with `ICONTOUR=1`: neighbours with data required (default 3) |
| `--contour-layers N` | with `ICONTOUR=1` and `3`: layers of points to rebuild (default 1) |
| `--contour-min-particles N` | with `ICONTOUR=2`: particles required around the point (default 1) |
| `--legacy-compat` | reproduce the behaviour of the Fortran version exactly |
| `--vtk` | also export to VTK for ParaView |
| `-q` | errors only |

`python -m pivnp path/to/case` works as well.

### From Python

```python
from pivnp import Simulation

sim = Simulation.from_directory("path/to/case")
summary = sim.run()

# Final state in memory, without reading files:
sim.particles.displacement   # (n, 2) accumulated displacement
sim.particles.eq_strain      # (n,)   equivalent shear strain
```

### Comparing results

```bash
python -m pivnp.compare reference.POST.RES new.POST.RES
```

Compares line by line and tolerates differences of one unit in the sixth digit.

## Input files

### `PIV-NP.TXT`

A single line with the case name (no extension), for example `zapatak`.

### `<case>.PAR`

Fortran list-directed format: values are separated by spaces, tabs or commas. The comment
lines are required and carry the field names above their values, so the file explains
itself.

```
Analysis title
BLOQUE 2: N_cel N_nod N_part_celda N_fil Ancho    Alto
          2006  2100  3            34    0.212115 0.212115
BLOQUE 3: del_t total_steps impresion moister version pivlab contour rec track
          0.8   149         1         0       1       1      0       0   0
BLOQUE 4: s_density porosity
          2650.0    0.4
```

| Parameter | Meaning |
|---|---|
| `NC`, `NN` | cells and nodes (points) of the PIVlab grid; NN = (NC/NFIL + 1)(NFIL + 1) must hold |
| `NPC` | rows and columns of particles per cell (NPC² per cell), from 1 to 6 |
| `NFIL` | rows of cells |
| `AXC`, `AYC` | cell width and height [m] |
| `DT` | time between images [s] |
| `TOTAL_STEPS` | number of PIVlab files to process |
| `IMPPAS` | results are written at step 1 and every `IMPPAS` steps |
| `MOISTER` | 0 = no moisture; 1 = read it from `Moist_<n>.TXT`; 2 = compute it from the test images (`<case>.HUM`) |
| `IVERSION` | 1 = PIVlab velocities sit at the grid nodes; 2 = they sit at the centre of each element (grid shifted half a cell) |
| `IPIVLAB` | 1 = 4-column files (older PIVlab); any other value = 5 columns |
| `ICONTOUR` | boundary correction: 0 = none, 1 = neighbour average, 2 = particle average, 3 = extrapolation |
| `IREC` | 0 = new analysis; 1 = continue from `<case>.REC` |
| `ITR` | must be 0 |
| `S_DENSITY`, `POROSITY` | read but not yet used in the computation (particle mass is 1) |

Input is validated on reading: if something does not add up (for instance, NC not a
multiple of NFIL) you get an explicit message instead of meaningless results.

#### Cases from earlier versions

Block 3 changed order several times over the years, and there were even files with the same
values arranged differently. There is a single input format, so an old case is moved to it
once:

```bash
pivnp <directory> --convert-par
```

This rewrites every `.PAR` below the directory and keeps each original next to it as
`<case>.PAR.orig`. The conversion reorders the values as they are, without reformatting
them, so not one decimal changes; the two places where it does correct something —`IPIVLAB`,
from the number of columns the case files actually have, and the `MOISTER` of the versions
where 2 meant something else— are reported on screen.

### PIVlab data

* `datos (1).txt`, `datos (2).txt`, …: PIVlab ASCII export with 3 header lines and columns
  `x, y, u, v` (or 5 columns when `IPIVLAB ≠ 1`). Points without a measurement hold `NaN`.
* `Moist_1.TXT`, … (only with `MOISTER=1`): 1 header line and columns
  `x, y, water content, saturation`.

File names are case-insensitive.

## Results

| File | Contents |
|---|---|
| `<case>.POST.MSH` | point mesh for GiD (one per particle), at their initial positions. Material 1 = active, 2 = no data at step 1 |
| `<case>.POST.RES` | per-particle results at every printed instant (GiD) |
| `<case>.REC` | final state, to continue the analysis with `IREC=1` |

Results in `.POST.RES`:

| Name | Type | Description |
|---|---|---|
| `Displacement` | vector | accumulated displacement |
| `Inst_displacement` | vector | displacement of the last step |
| `NaNs` | scalar | missing data around the particle: element nodes without a measurement (IVERSION=1) or central point without one (IVERSION=2) |
| `Velocity`, `Acceleration` | vector | particle velocity and acceleration |
| `Total_strain`, `Inc_strain` | 3 comp. | εxx, εyy, γxy accumulated / per step |
| `Equi_strain`, `In_E_strain` | scalar | equivalent shear strain, accumulated / per step |
| `Vol_strain`, `Ins_vol_strain` | scalar | volumetric strain, accumulated / per step |
| `E_potential`, `E_kinetic`, `E_total` | scalar | energies per unit mass (`E_total` = potential + kinetic) |
| `Moisture`, `Saturation` | scalar | only with `MOISTER=1` |

Since the mesh holds the initial positions, `mesh + Displacement` is always the current
position of each particle.

### Continuing an analysis

With `IREC=1` the analysis carries on from the previous `.REC`: time and step numbering
continue, and the new results are **appended** to the existing `.POST.RES`. The PIVlab
files of the continuation are numbered from 1 again. Splitting an analysis in two gives
exactly the same result as running it in one go.

## How it works

```
 images ──PIVlab──► velocities on a fixed grid ──PIV-NP──► trajectories and strains
                    (datos (n).txt)                        of material particles
```

**1. Grid and particles.** The PIVlab points form a rectangular grid of
(NC/NFIL + 1) × (NFIL + 1) nodes. NPC × NPC numerical particles are created in each cell,
either at the Gauss points (NPC from 4 to 6) or evenly spread (NPC 2 and 3). With
`IVERSION=2` the grid is shifted half a cell, so that each PIVlab point sits at the centre
of an element and shares its velocity among the 4 nodes of that element.

**2. Nodal velocities.** At each step *n*, `datos (n).txt` is read, the points are
reordered (PIVlab numbers them by columns and top-down; PIV-NP by rows and bottom-up) and
the sign of *v* is flipped (the image *y* axis points downwards).

**3. Interpolation to the particles.** Each particle is located in its cell and its local
coordinates (ξ, η) ∈ [−1, 1] are computed. With the bilinear shape functions

  Nᵢ(ξ, η) = ¼ (1 + ξ ξᵢ)(1 + η ηᵢ)

the code obtains the velocity **v**ₚ = Σ Nᵢ **v**ᵢ, the acceleration
**a**ₚ = Σ Nᵢ (**v**ᵢⁿ − **v**ᵢⁿ⁻¹) / Δt and the step displacement Δ**u**ₚ = **v**ₚ Δt.

**4. Strains.** Using the shape function derivatives at the element centre (strain is
uniform within each cell):

  Δεxx = Σ ∂Nᵢ/∂x · vx,ᵢ Δt  Δεyy = Σ ∂Nᵢ/∂y · vy,ᵢ Δt  Δγxy = Σ (∂Nᵢ/∂y · vx,ᵢ + ∂Nᵢ/∂x · vy,ᵢ) Δt

They are accumulated, and the volumetric strain εv = εxx + εyy and the equivalent shear
strain εq = ⅔ q are computed, with q = √(3 J₂) of the tensor
(εxx, εyy, εzz = 0, εxy = γxy/2).

**5. Update.** The particle moves to **x**ₚ + Δ**u**ₚ. If it leaves the grid it is no
longer computed. Particles in cells whose 4 nodes have no data at the first step are taken
to be outside the material (air, background) and are excluded from the results.

**6. Results.** All results are written at step 1 and every `IMPPAS` steps; the final state
is saved to `<case>.REC`.

## Visualisation

Recommended: [ParaView](https://www.paraview.org/download/), free and open source. PIV-NP
converts the results to VTK, its native format, including those of the Fortran version:

```bash
pivnp-vtk path/to/case/zapatak.POST.RES
```

Then, in ParaView: **File → Open → `zapatak_vtk/zapatak.pvd` → Apply**, colour by
`Equi_strain` and press **Play**. Full guide (in Spanish) in
[`docs/VISUALIZACION.md`](docs/VISUALIZACION.md).

## Boundary correction

At the edge of the material PIVlab reports no velocity (NaN) for the points whose
interrogation window falls partly outside. Those points enter the interpolation without a
measured value, so the particles at the boundary move less than they should and some end up
detaching from the material.

`ICONTOUR` selects how the velocity of those points is rebuilt:

| ICONTOUR | Method | Parameters |
|---|---|---|
| 0 | none | — |
| 1 | average of the neighbouring nodes with data (8 neighbours) | `--contour-min-neighbors`, `--contour-layers` |
| 2 | average velocity of the particles in the surrounding elements | `--contour-min-particles` |
| 3 | linear extrapolation from the interior outwards | `--contour-layers` |

With any method other than 0, the shifted grid also normalises how each point is shared
among the nodes of its cell, so boundary nodes get the average of the contributing points
instead of a fraction of it.

Rebuilt points are flagged separately and do **not** change the "no data" marks, so the
correction never activates particles in the air: it only improves the interpolated velocity
of the boundary particles.

Effect on the centrifuge case (149 steps, 9378 active particles):

| Configuration | Isolated particles at the end | Maximum displacement |
|---|---|---|
| No correction | 14 | 3.61 m |
| ICONTOUR=1, 3 neighbours, 1 layer | 0 | 3.67 m |
| ICONTOUR=2, particle average | 0 | 3.10 m |
| ICONTOUR=3, extrapolation | 0 | 3.39 m |

"Isolated particles" are those left without any neighbour within one cell, that is, those
that detached from the material.

![Comparison of the boundary correction methods](docs/img/contorno_comparativa.png)

Which one is best is a physical question rather than a programming one: it is worth
checking against PTV markers or photographs of the test. Adding another method only takes a
class with an `apply(ctx)` method in `pivnp/contour.py`.

## Performance

Centrifuge case (18 054 particles, 149 steps, 768 MB of results), on a 16-core processor:

| Version | Time | Speed-up |
|---|---|---|
| Fortran (gfortran -O2) | 38.0 s | 1× |
| Python, 1 thread | 3.1 s | 12× |
| Python, 4 threads | 2.3 s | 17× |
| Python, 16 threads | 2.0 s | 19× |

Excluding the first Numba compilation (about 6 s, once). To measure it on your machine:

```bash
python benchmarks/benchmark.py path/to/case --threads 1 4 0 --legacy legacy/pivnp_legacy.exe
```

Where the speed-up comes from is explained (in Spanish) in
[`docs/ARQUITECTURA.md`](docs/ARQUITECTURA.md). The comparison uses the Fortran code built
with gfortran.

## Tests

```bash
pytest
```

* **Synthetic cases** (`tests/data/sinteticos/`): velocity fields with a known solution
  (displacement without strain, uniform and varying horizontal strain, shear and rigid-body
  rotation). The computed strains are checked against the analytical ones, and the analyses
  the group ran at the time are reproduced.
* **Regression**: 7 scenarios also run with the Fortran version (both grids, NPC = 2, 3 and
  4, water content, 5-column format and restart). In `--legacy-compat` mode, `.POST.RES`,
  `.POST.MSH` and `.REC` are required to be **byte-for-byte identical**.
* **Behaviour**: properties checked against known solutions (constant acceleration, linear
  strain field, an analysis split in two being equivalent to a single run).
* **Unit**: each module is compared with a literal translation of the Fortran loops
  (`tests/legacy_reference.py`) or with cases whose solution is known.

Coverage (Numba kernels can only be measured with compilation disabled):

```bash
NUMBA_DISABLE_JIT=1 pytest --cov=pivnp
```

In PowerShell: `$env:NUMBA_DISABLE_JIT=1; pytest --cov=pivnp`. Current coverage: 96 %.

To regenerate the reference results with the Fortran version (requires gfortran):

```bash
python tools/make_reference.py --source "path/to/caso paper centrifuga"
```

## Project layout

```
piv-np/
├── src/pivnp/
│   ├── config.py          reading and validation of PIV-NP.TXT and .PAR
│   ├── mesh.py            grids and particle location
│   ├── particles.py       particle generation
│   ├── state.py           particle and nodal arrays
│   ├── pivlab_io.py       reading of PIVlab and water content data
│   ├── nodal.py           nodal fields
│   ├── contour.py         boundary correction
│   ├── solver.py          motion and strains
│   ├── gid_writer.py      GiD output
│   ├── fortran_format.py  fast, exact I/E14.6 formatting
│   ├── restart.py         .REC file
│   ├── simulation.py      main loop
│   ├── compare.py         result comparison
│   ├── vtk_export.py      VTK export for ParaView
│   └── cli.py             command line
├── tests/                 unit, behaviour and regression tests
├── legacy/                original Fortran code, unmodified
├── tools/                 generation of reference results
├── benchmarks/            performance measurement
├── examples/              centrifuge case
└── docs/
    ├── ARQUITECTURA.md    design and improvement proposals (Spanish)
    └── VISUALIZACION.md   ParaView guide (Spanish)
```

Source code comments and the documents under `docs/` are written in Spanish.

## License and citation

BSD 4-clause (see [`LICENSE`](LICENSE)). All advertising materials mentioning features or
use of this software must cite:

> Pinyol, N.M. & Alvarado, M. (2017). Novel PIV-based analysis for large displacement.
> Canadian Geotechnical Journal 54(7): 933-944.
