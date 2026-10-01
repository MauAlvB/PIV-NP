# Architecture: current state and proposals

## 1. What has been done

The Fortran program was a single file with 7 subroutines that talked to each other through
global `COMMON` blocks (more than 60 fixed-size arrays). Version 2.0 splits the
responsibilities into small modules that are tested separately:

```
             ┌────────────┐   ┌──────────────────────┐
 .PAR ─────► │ config.py  │   │ sources.py           │ ◄── datos (n).txt, Moist_n.TXT
             └─────┬──────┘   │  pivlab_io.py (PIVlab)│    (read ahead on worker threads)
                   │          │  ...another source    │
                   │ CaseConfig└──────┬───────────────┘
                   │                  │ Frame
                   ▼                  ▼
             ┌──────────────────────────────────┐
             │ simulation.py   (loop of steps)  │
             └──┬──────────┬──────────┬──────┬──┘
                │          │          │      │
      mesh.py   │  nodal.py│ contour.py  solver.py        restart.py
      particles │  (nodes) │ (CONTOUR)   (particles,      (binary .REC,
      state.py  │          │             parallel Numba)   compatible)
                ▼
          gid_writer.py + fortran_format.py  ──► .POST.MSH / .POST.RES
          (results table, parallel formatting, written on another thread)
```

| Fortran | Python |
|---|---|
| `COMMON /PARTICULAS/`, `/NODOS/` | `state.Particles`, `state.Nodes` (named arrays) |
| `COMMON /GEOMETRIA/` + `UCELDA` | `mesh.Grid`, `mesh.locate_point` (O(1) instead of a linear search) |
| `PIVLAB_DATA` | `config.parse_par`, `particles.create_particles`, `simulation._load_restart` |
| `VELOCIDADES` | `sources.DisplacementSource` (`pivlab_io.PivlabSource`) + `nodal.*` |
| `CONTOUR` | `contour.ContourCorrection` (interface) |
| `SOLMOV`, `INVAR2` | `solver.advance_particles`, `solver.update_strains`, `solver.deviatoric_q` |
| `IMPRES_GiD` (16 copied loops) | `gid_writer.RESULTS` (a table) + `GidWriter` |
| `RECOM` | `restart.write_restart` |

Where the speed-up comes from (38 s → 2 s on the centrifuge case):

1. **Writing the results** (what cost the most): the original formats number by number
   with `WRITE`. Now every block is formatted in parallel straight into a byte buffer and
   the disk is written on a separate thread while the next step is being computed.
2. **Locating the particles**: `UCELDA` walked rows and columns (O(rows+columns)) and was
   called about 20 times per particle and step (once per results block). It is now O(1) and
   done once per phase.
3. **Per-particle loops** in parallel with `numba.prange` (each particle is independent).
4. **Reading ahead** the PIVlab files on worker threads.

The loops that accumulate over nodes (nodal mass, distribution on the staggered mesh) are
deliberately left sequential: summing in a different order would change the rounding and
would stop matching the original byte for byte.

## 2. Proposals for going further

Ordered by benefit/effort ratio.

### P1. Compatibility mode with the original code — **done**

`RunOptions(legacy_compat=True)` (`pivnp --legacy-compat`) reproduces the behaviour of the
Fortran exactly, so the regression tests keep comparing against it byte for byte and the
old analyses can be repeated with this version.

### P2. Self-describing configuration file — **partly done**

The `.PAR` is positional: one value out of place changes the meaning of everything that
follows, and there used to be several dialects (the one in the code, the one in the
tutorial and the one of the example case). Done so far: every case was rewritten in a
**single format** with the field names above their values and a legend after block 4
(`pivnp --convert-par`, which keeps the original as `.PAR.orig`).

Still pending: a file with named keys, which would make the order irrelevant. Something
like a `case.toml`:

```toml
[mesh]
cells = 2006
nodes = 2100
rows = 34
cell_width = 0.212115    # m
cell_height = 0.212115   # m
particles_per_side = 3
version = 1

[analysis]
dt = 0.8                 # s
steps = 149
print_every = 1

[contour]
method = "neighbour_average"   # "none" | "neighbour_average"
min_neighbours = 3
layers = 1
```

### P3. Interchangeable strategies for the physical decisions (medium effort)

In the original, the alternatives were handled by commenting lines out (for example, the
strain at the centre of the element versus the strain at the position of the particle, or a
mass of `1.0` versus `S_DENSITY*(1-POROSITY)*VVP`). Proposal: turn every decision into an
interface with named implementations, chosen from the configuration, just as has already
been done with `ContourCorrection`:

| Decision | Options | State |
|---|---|---|
| Contour correction | none · neighbour average · particle average · extrapolation | **done** (`ICONTOUR`) |
| Particle mass | unit (current) · dry density × volume (`S_DENSITY`, `POROSITY`) | pending |
| Strain computation | centre of the element (current) · at the position of the particle | pending |
| Initial particle layout | Gauss (NPC 4-6) · uniform (NPC 2-3) | fixed by NPC |

That way every variant is documented, tested and selectable without touching the code, as
already happens with `ContourCorrection`.

### P4. Modern output formats (medium effort, large gain)

The text `.POST.RES` of the example case takes **768 MB**. Proposal: a `ResultWriter`
interface with several implementations:

* `GidAsciiWriter` (the current one, for compatibility).
* `Hdf5Writer` / `NetCdfWriter`: compressed binary, typically 10 to 20 times smaller; opens
  directly from Python, MATLAB or ParaView.
* `VtkWriter` (`.vtu` per step + `.pvd`): visualisation in ParaView without a GiD licence.
* `TrajectoryCsvWriter` with the trajectory of specific points, to compare against
  laboratory PTV tracking markers.

### P5. Data input decoupled from PIVlab — **done**

Where the displacements come from is now an explicit, interchangeable contract. The solver
never sees a file: it only ever sees one `Frame` per step.

[`sources.py`](../src/pivnp/sources.py) holds `DisplacementSource`, the five members a
source has to offer, and the registry that `--source` selects from. PIVlab is one
implementation of it (`PivlabSource` in `pivlab_io.py`), registered under `pivlab`, which is
the default. Reading the interval between images from the PIVlab header — the last thing
that tied the analysis to the format — now sits behind `frame_interval()`.

A new source is one class and one registration:

```python
import numpy as np
from pivnp.sources import Frame, register_source


class MySource:
    """Velocities straight from a simulation: no files involved."""

    def __init__(self, n_points: int) -> None:
        self.n_points = n_points
        self.moisture = False   # are there Moist_<n>.TXT to read?
        self.images = None      # moisture from images; the analysis sets this

    def frames(self, steps):
        for step in steps:                      # in order: the moisture model has memory
            u = np.full(self.n_points, 0.2)     # m/s at every grid point, NaN if missing
            zeros = np.zeros(self.n_points)
            yield Frame(step, f"my data {step}", u, zeros, zeros, zeros)

    def mesh_in_metres(self, step=1):           # only needed with MOISTER=2
        x = y = np.zeros(self.n_points)
        return x, y, 0.001, np.ones(self.n_points, dtype=bool)

    def frame_interval(self):
        return None                             # unknown: the DT check is skipped


@register_source("mine")
def _build(case_dir, config, prefetch=4):
    return MySource(config.n_nodes)
```

Then `pivnp <case> --source mine`. The points must come in the order of the PIVlab export
(by columns, y growing downwards), which is the order `pivlab_to_node` maps to PIV-NP nodes.
`tests/test_sources.py` runs a whole analysis through a source like this one, so the seam
stays checked.

Still worth adding, now that it is cheap: a reader for PIVlab's `.mat` files (which would
avoid exporting thousands of `.txt`), OpenPIV, DaVis, and a PIV analysis built into PIV-NP
so that no external package is needed at all.

### P6. API for notebooks (low effort)

`Simulation` already exposes `particles` and `nodes` in memory. Adding a method that
returns the results as an `xarray.Dataset` (dimensions: time × particle) would allow
analysing and plotting in Jupyter without going through files.

### P7. Continuous integration and publishing (low effort)

* GitHub/GitLab with CI running `pytest` and `ruff` on every change (Windows and Linux).
* `pre-commit` with `ruff format`.
* Semantic versioning, a `CHANGELOG.md` and publishing on PyPI or conda-forge.
* A Zenodo DOI for each version, citable alongside the 2017 paper.

### P8. Documentation (low effort)

The modules already carry docstrings. With MkDocs + mkdocstrings a website can be generated
with the API reference, the laboratory tutorial (the current `.docx`) and an example
notebook with the centrifuge case.

## 3. Suggested roadmap

1. **P2** (named keys) + **P3**: configuration with names and selectable computation options.
2. **P4** (HDF5/VTK) and **P7** (CI) to share the code with other groups.
3. A PIV analysis built into PIV-NP, as a source (the seam of P5 is already there), so that
   a case can be run without an external PIV package.
4. **P6** and **P8** as the group needs them.
