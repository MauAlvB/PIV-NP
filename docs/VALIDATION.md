# Validation

How this version was checked against the earlier ones, what came out of it, and what is
still open. It is here so that anyone can see on what grounds the results are trusted, and
so that the known limits of the method are not rediscovered from scratch.

Every case was run on a hard-linked copy in a temporary directory, never on the original
folder: an analysis writes `.POST.MSH`, `.POST.RES` and `.REC`, and the test data is never
touched.

## 1. The port is identical to the original Fortran

**The Python version is byte-for-byte identical to the compiled 2024 Fortran.** Checked by
running both programs on three different cases — PIV-NP mesh and staggered mesh, 20 and 149
steps, from 5432 to 13806 nodes — and comparing the fifteen results and the mesh: not one
difference.

That is what makes the rest of this document meaningful. When a stored result does not
match, the question is no longer "is our code right?" but "which version produced that
analysis?".

The regression tests (`tests/test_regression.py`) keep that property under guard: seven
scenarios, generated with the Fortran executable by `tools/make_reference.py`, are required
to come out identical byte for byte in `--legacy-compat` mode.

## 2. Reproducing the old analyses

| Case | Result |
|---|---|
| Dam failure, stage 1 | **Exact** |
| Dam failure, stage 1, staggered mesh | **Exact** with `--legacy-2023-average` |
| Dam failure, stage 2 (restart) | **Exact** |
| Dam failure, stage 2 (restart, staggered mesh) | Exact in every value; 4 particles of difference in which ones get printed |
| Centrifuge | First step exact with `--legacy-2023-average`; after that 10 to 15 % of the trajectories diverge |

### Finding: the `.PAR` had no fixed order

Block 3 was laid out differently by each version of the program, and the *same number of
values* means different things depending on the year:

| Dialect | Block 3 | Block 4 |
|---|---|---|
| 2024 | `DT STEPS IMPPAS MOISTER IVERSION IPIVLAB ICONTOUR IREC ITR` | yes |
| Dam failure, stages | `DT STEPS IMPPAS MOISTER S_DENSITY POROSITY IVERSION PTV IREC` | no |
| Slope RGB | `DT STEPS IMPPAS MOISTER S_DENSITY POROSITY IVERSION PTV` | no |
| Dam failure, tests | `DT STEPS IMPPAS MOISTER IVERSION IPIVLAB IREC PTV` | yes |
| Centrifuge 2022 | `DT STEPS IMPPAS IVERSION IPIVLAB MOISTER IREC ITR` | yes |
| Oldest | `DT STEPS IMPPAS [MOISTER]` | no |

The last two carry the same eight values and the same block 4 and mean different things: by
position alone they cannot be told apart.

**The fix was a single input format**, not a cleverer reader. Files are converted once with
`pivnp <directory> --convert-par`, which keeps each original as `<case>.PAR.orig`, and from
then on there is one format, with the field names above their values and a legend after
block 4. The 2024 Fortran still reads it, which is what allows the comparison against it to
continue.

To read the old files, the converter does look at the comment line, which is where each
dialect documents its own order — the only way to separate the two eight-value dialects.
That logic now lives only in the converter; the reader accepts the single format and, when
it sees anything else, says so and points at the converter instead of guessing.

**The conversion changes no result**, and that is checked with a simple guarantee: reading
the old file and reading the converted one must give the same configuration, field by
field. It was verified on all 46 `.PAR` files available. Two deliberate corrections fall
outside that and are always reported:

* **`MOISTER=2` of the 2022 centrifuge**, where it meant "read the moisture files" and in
  the current format means "compute it from the images". It is written as 1.
* **`IPIVLAB` of the dam-failure stages**, whose `.PAR` does not carry the field, so 1
  (4-column files) was assumed while theirs have 5 columns. The value is measured by
  counting the columns of the case's own `datos (1).txt`. These are the only two files of
  the 46 whose configuration changes, and the change is exactly the one that was needed:
  with it, that case gives exactly the same as the 2024 Fortran.

A `DT` that does not match the PIVlab header is reported but **not touched**: changing it
would change the results, and that is a separate decision.

### Finding: the distribution on the staggered mesh changed between versions

Results stored with `IVERSION=2` did not match. The ratio between their values and ours
came out as 4/3 and 3/2, which is what you get by dividing by the number of points
contributing to each node. With that average **all fifteen results match exactly**, with
one detail: that version did not apply the average to the momentum increment, so its
acceleration was not consistent with its velocity. Reproducing that too, the agreement is
complete.

The change is visible in the 2024 Fortran, with the original line commented out beside it:

```fortran
F1=0.25d0*1.0D0 !AMASSINI(I) ... !AM(JJ)      !Cuidado con la masa
AMASSNODE(JJ)=1.0D0 !AMASSNODE(JJ)+0.25d0*AMASSINI(I)
```

The accumulation of the nodal weight was replaced by the constant 1, and the average went
with it. That is the origin of the under-valued velocity at the boundary nodes: not a
long-standing oversight but a regression introduced between versions — and it confirms that
averaging is what used to be done.

`--legacy-2023-average` exists to repeat those analyses.

### Finding: the PIVlab file format cannot be taken from the `.PAR`

`.PAR` files from before 2024 carry no `IPIVLAB`, so it defaulted to 1, which means reading
4·NN values in a row without looking at the line breaks. Files with five columns shift
every value out of place from the first one on, and the analysis comes out meaningless. The
format is now deduced from the file itself. With four-column files both paths agree, so no
analysis that already worked changes.

### Note on the oldest stored analyses

The original dam-failure analyses were made with a version well before the others.
Comparing the particles both print, the displacements, velocities, energies, moisture and
saturation **match exactly**. They differ in that:

* that version printed all 13552 particles, including the ones outside the material; this
  one prints only the 5916 that are located;
* the strains agree at step 1 and diverge afterwards, which points to a different criterion
  for marking a particle as NaN rather than to a different formula.

## 3. Moisture from the images (`MOISTER=2`)

The reference is the test of the published paper, which has SWIR images and the moisture
analysis made from them. **99.8 % of the nodes come out with the exact moisture**, and the
saturation with a mean difference of 0.0024. Getting there needed four things that were not
written down anywhere:

1. **Which calibration table.** The (moisture, saturation) pairs of the case did not fall on
   the table used for the RGB flow. Since both fields come from the same gray level through
   the same table, the pairs identify it: the moisture matches to the last decimal.
2. **How wide the band is.** Across the 40 instants only **19 distinct pairs** appear, and
   solving for the normalised gray of each gives 0, 2.5, 5, 7.5 … 45: a perfect staircase of
   2.5. As the image gray is an integer, that fixes the band width at **exactly 40 levels**
   (100/2.5). It also confirms that the old code did round the filter, which is what
   quantises the values.
3. **The registration between the two cameras.** MATLAB solved it with a four-point
   homography marked by hand, which was not saved. It was recovered from the data itself:
   the gray measured on the SWIR image has to be `saturated + gv_n/2.5`, so one looks for
   the homography leaving the least scatter in that difference. The residual is **0.294 gray
   levels**, below the quantisation, and three different starting points converge on it.
4. **The band in numbers.** Saturated **92**, dry **132**.

That case's moisture field is correctly computed at every node and carries **no incremental
policy**: no node has saturation exactly 1, and the minimum falls and rises again through
the test. It is the method of the published paper, which is why it is the reference.

### A bug in the MATLAB code the old data carries

In the RGB case the moisture was computed **in a single column of nodes**: in `gv2sat_JC2.m`
the moisture computation sits *outside* the inner loop, so it runs once per row with the
last value of the column. Verified in the case files: only 34 of 2100 nodes (1.6 %) have a
non-zero moisture, and they are always the last ones of the file. The saturation is computed
correctly at every node, because that computation is inside both loops.

So **the moisture field of those analyses is essentially empty**. This implementation
computes it at every node, which is what the code meant to do.

## 4. How much each decision moves the result

Measured on the paper case, the only one with a reference computed without a ratchet and
without a threshold. The reference variant is the one that reproduces its files: global band
92–132 (40 levels), σ = 40, rounded filter and no incremental policy. Each variant changes
one single thing, compared over the 40 instants and every node with data.

| Variant | mean saturation diff. | p95 | max | mean moisture diff. | nodes moving > 0.05 | final mean saturation | nodes with sat > 0.99 |
|---|---|---|---|---|---|---|---|
| *(reference)* | — | — | — | — | — | 0.937 | 80.4 % |
| Ratchet, threshold 0.8 | 0.022 | 0.131 | 0.53 | **0.000** | 11.3 % | 0.972 | 92.7 % |
| Ratchet, threshold 0.95 | 0.015 | 0.056 | 0.43 | **0.000** | 6.0 % | 0.968 | 88.8 % |
| Ratchet, no threshold | 0.008 | 0.056 | 0.43 | **0.000** | 5.2 % | 0.961 | 85.3 % |
| Band of 11 levels | **0.105** | 0.613 | 0.68 | 2.52 | 24.0 % | 0.856 | 80.4 % |
| Band of 20 levels | 0.062 | 0.383 | 0.44 | 1.47 | 21.6 % | 0.888 | 80.4 % |
| Band of 60 levels | 0.028 | 0.234 | 0.26 | 0.59 | 15.5 % | 0.961 | 80.4 % |
| Band of 80 levels | 0.041 | 0.328 | 0.41 | 0.82 | 17.2 % | 0.972 | 80.4 % |
| Band shifted −2 levels | 0.021 | 0.095 | 0.12 | 0.52 | 21.6 % | 0.920 | 74.8 % |
| Band shifted +2 levels | 0.017 | 0.095 | 0.12 | 0.40 | 17.2 % | 0.950 | 84.9 % |
| Band shifted ±4 levels | 0.031–0.047 | 0.19–0.20 | 0.22 | 0.73–1.20 | 17–27 % | 0.900–0.961 | 62–88 % |
| Filter not rounded | **0.002** | 0.016 | 0.03 | 0.06 | **0.0 %** | 0.937 | 79.9 % |
| σ = 15 instead of 40 | 0.035 | 0.219 | 0.71 | 0.83 | 17.7 % | 0.915 | 77.4 % |
| σ = 80 instead of 40 | 0.047 | 0.351 | 0.62 | 1.05 | 16.5 % | 0.974 | 83.9 % |
| Per-node reference (+5/−6) | 0.152 | 0.722 | 0.99 | 3.41 | 38.0 % | 0.885 | 71.2 % |
| Per-node reference, band of 40 | 0.663 | 0.980 | 0.99 | 15.61 | 95.4 % | 0.296 | 6.0 % |

What comes out of it:

**1. The width of the band dominates everything else.** Going from 40 levels to the 11 of
the RGB flow moves the saturation by 0.105 on average and 0.61 at the 95th percentile: an
order of magnitude more than any other decision. Being wrong by a mere factor of two (20
instead of 40) still moves it by 0.062. It is by far the number that has to be measured
well.

**2. The position of the band weighs as much as the ratchet.** Shifting it by **two gray
levels** — what separates two neighbouring clicks when marking the reference by hand —
changes the mean saturation by 0.02 and moves the fraction of nodes that end up saturated
by ten points, from 75 % to 85 %. Those two ends currently come from picking two points on
an image and adding four constants written in the code. They should come from a measurement,
not from a click.

**3. The ratchet moves the saturation and leaves the moisture exactly as it was: 0.000.**
That is the inconsistency between the two fields, quantified: up to 0.53 of difference in
saturation without the moisture changing by a thousandth. A threshold at 0.8 doubles the
effect of the ratchet alone (0.022 against 0.008) and takes the nodes declared saturated
from 80 % to 93 %.

**4. Rounding the filter does not matter here**: 0.002 on average and not one node above
0.05. And yet in the RGB case that same rounding was *all* of the residual. The difference
is the band width: with 40 levels one gray level is 2.5 % of the scale; with 11 it is 9 %.
The rounding was never the problem, it was a symptom of too narrow a band.

**5. σ is not a cosmetic parameter.** Changing it to 15 or to 80 moves the result by 0.035
to 0.047, as much as being 30–50 % wrong about the band width. It deserves to be justified
by the grain size and the image scale, not picked by eye.

**6. The per-node reference depends entirely on the reference image being in a known
state.** With the only candidate available in that case — an image where the registration
points were marked, not a photograph of dry soil — the result falls apart: 0.66 of mean
difference and 95 % of the nodes affected. That does not say the per-node method is worse
(it cancels the fixed texture pattern, which is its advantage); it says the method
**requires an experimental condition that is written down nowhere**.

### What was changed with those measurements in hand

* **The incremental policy now acts on the normalised gray**, not on the saturation, so both
  fields come from the same value and cannot contradict each other. Checked on the RGB case:
  the saturation comes out **exactly as before** (maximum difference 0 at every step, so no
  compatibility switch is needed) and the 579 nodes declared saturated go from a median
  moisture of 1.14 % to 24.03 %, that of the saturated soil.
* **The ratchet is off by default** (`INCREMENTAL = 0`) and its threshold is 0.95. It is a
  hypothesis about the test rather than a measurement, and with it on the drying branch
  cannot be measured at all.
* **A warning when the band is narrow**, below 20 gray levels, where one level weighs more
  than 5 % of the scale. The RGB case, with 11, triggers it.
* **A quality mark per node**: measured, no data, or sitting at one end of the band, where
  the value is a bound and not a measurement. The analysis reports the breakdown when it
  finishes. On the paper case the answer is uncomfortable and worth having in plain sight:
  **82 % of the values with data sit at the wet end**, which makes them an "at least" and not
  a measurement. Only 18 % are measurements.

### Priorities this sets

1. The width of the band, from a column test, per soil.
2. The position of the band, measured and not clicked.
3. σ, justified.
4. Keeping moisture and saturation consistent, and keeping the ratchet out of the
   measurement.
5. The rounding, which takes care of itself as soon as the band is wide.

All of this is measured on the SWIR case, with a band of 40 levels. In the RGB flow, with
11, every sensitivity is larger.

## 5. The boundary correction (`ICONTOUR`)

The difficulty in measuring this on real data is that at the boundary there is no truth to
compare against: that is exactly where PIVlab does not measure. The way out is to **hide
nodes that do have a measurement and that sit next to the boundary** — their neighbourhood
is as incomplete as that of a real boundary node — rebuild them with each method and compare
against what PIVlab measured. It is the same situation the correction has to solve, but with
an answer.

The baseline is not correcting, which leaves the node at zero: being wrong by the whole
velocity it should have had. The deciding column is **in what fraction of the nodes the
rebuild lands closer to the real value than that zero does**.

| Case | 1 · neighbour average | 3 · extrapolation | 2 · particle average |
|---|---|---|---|
| Dam failure, stage 1, step 10, IVERSION=1 | **80 %** · 0.52× | 69 % · 0.66× | 29 % · 1.48× |
| Dam failure, stage 1, step 10, IVERSION=2 | **80 %** · 0.52× | 69 % · 0.66× | 30 % · 1.46× |
| Dam failure, stage 2, step 18 | **76 %** · 0.61× | 56 % · 0.74× | does not rebuild |
| Centrifuge, step 75 | **96 %** · 0.05× | 94 % · 0.10× | 92 % · 0.41× |
| Slope RGB, step 80 | **98 %** · 0.13× | 97 % · 0.15× | 50 % · 0.75× |

(the second number is the median error, in multiples of the typical boundary velocity;
without correction it is 1.00× by construction)

**The neighbour average wins on all five.** It beats leaving the zero between 76 % and 98 %
of the time and cuts the median error to between a half and a twentieth. It also always
rebuilds: it never ran out of neighbours in any case tried.

**Extrapolation comes second**, very close where the field is smooth and clearly behind
where the motion is abrupt (76 % against 56 % on stage 2 of the dam failure). That makes
sense: continuing a slope amplifies whatever noise there is at the boundary.

**The particle average is the worst and often does harm.** On the dam failure it does worse
than not correcting on 70 % of the nodes, and on stage 2 it rebuilds nothing. The reason is
structural: it averages velocities of particles that were themselves interpolated with the
boundary nodes set to zero. It is circular, and feeds back the very error it was meant to
remove. The program warns when `ICONTOUR=2` is selected.

The benefit also has a clear dependence: the smoother the field near the boundary, the more
the correction gains. On the centrifuge the median error drops to 0.05× of what not
correcting cost; on the dam failure, where the motion is abrupt, it stays at 0.52×.

## 6. Rotation, measured against a known answer

A rigid rotation does not deform the material, but accumulating linear strain increments
reports strain for it. The synthetic `rotation` cases pin that down: 50 steps of 1° leave
εxx = εyy = n·(cos Δθ − 1) and an equivalent shear strain of about 0.005 where there is no
deformation at all. Nothing in the results told the reader that.

The curl of the velocity field does, and it comes nearly free: shear and rotation are the
symmetric and the antisymmetric halves of the same velocity gradient, which the solver
already evaluates to get the strains. Two results were added from it, `Vorticity` (s⁻¹) and
`Rotation` (degrees accumulated).

Measured on the cases where the answer is known:

| Case | True value | Measured | Error |
|---|---|---|---|
| Simple shear at 0.05 s⁻¹ for 2 s (`shear-block`) | ω = −0.050000 s⁻¹, −2.864789° | identical to every digit | 0 |
| Rigid rotation, 1°/s for 50 s (`rotation_1P`), boundary rebuilt | 50° | 49.997°, the same on every particle | **0.003°** |
| the same case, boundary not corrected | 50° | 49.12° on average, spread 48.2 to 50.0 | 0.9° |

So on the very case where the strain claims a deformation that is not there, the rotation
lands within three thousandths of a degree of the truth. The second and third rows are also
one more measurement of what the boundary correction is worth.

The accumulated rotation inherits the first-order nature of the formulation, which is why it
was checked against a finite 50° rotation rather than assumed; it is carried in the `.REC`
(extension version 3), so a restarted analysis does not start turning from zero again.

Neither result is written in `--legacy-compat` mode: the original Fortran had no such blocks,
and the regression suite compares the whole file against it.

### The rotation artifact is not shear: it is an apparent contraction

Breaking the spurious strain of that case into its components turned out to matter more than
the magnitude. On a material that does not deform at all:

| Component | Value |
|---|---|
| εxx | −0.0076152 |
| εyy | −0.0076152 |
| γxy | **0.0000000** |
| `Vol_strain` | **−0.0152305** |

The two normal strains are **equal** and the shear is **exactly zero**. So the artifact is an
apparent isotropic in-plane contraction, and what it shows up in is **`Vol_strain`: a volume
loss of 1.5 % that never happened**. `Equi_strain` reports 0.0051 only because the deviatoric
invariant includes εzz = 0, which makes an in-plane isotropic contraction deviatoric in 3D.

This is worth knowing before reading a volumetric map over a rotating zone. In soil
mechanics a contraction means densification and pore-pressure build-up, so a rotating block
can be mistaken for a compacting one. The ratio of rotation to shear is what tells them
apart: on this case the deviatoric shear is zero to machine precision, which makes
`Vorticity_num` infinite and `Rot_angle` exactly 90° — "this is a rotation", stated as
plainly as it can be.

### How the deformation splits, on a real case

`Vorticity_num` and `Rot_angle` are the kinematic vorticity number of the accumulated
deformation and its arctangent. Measured on the dam-break case at the last step, keeping only
the particles that actually deformed:

| How much the material deformed | Particles | Median vorticity number | Quartiles |
|---|---|---|---|
| `Equi_strain` < 0.01 (practically intact) | 553 | 1.31 | 0.60 – 2.44 |
| 0.01 – 0.025 | 1971 | 0.60 | 0.27 – 1.19 |
| 0.025 – 0.05 | 1301 | 0.51 | 0.24 – 0.94 |
| 0.05 – 0.10 | 1059 | 0.64 | 0.30 – 0.98 |
| over 0.10 (well deformed) | 1032 | 0.60 | 0.30 – 0.88 |

The first row is the warning: where nothing happened the number is a ratio of two near-zero
quantities and means nothing, and it is *also* where it looks most dramatic. From the second
row on it settles around 0.5–0.6, between pure and simple shear, with half the scatter.

Read regionally rather than particle by particle — even in well-deformed material the
interquartile range is 0.30 to 0.88 — the map does separate the upstream face, which slides
with little rotation, from the downstream material, which rotates as it collapses.

## 7. The finite strain, measured against the incremental one

The strain the original computes is a sum of linear increments, one per step. The `Finite_*`
results differentiate once instead, over the accumulated displacement, from the deformation
gradient `F = I + ∂u/∂X`. Both are published; these are the measurements that justify having
the second.

| | incremental | finite | truth |
|---|---|---|---|
| **`rotation_1P`** — turns 50°, deforms nothing | | | |
| equivalent shear | 0.005077 | **0.000000** | 0 |
| volumetric / area change | −0.015230 | **0.000000** | 0 |
| rotation | −49.9975 | **−50.0000** | 50 |
| **`shear-block`** — simple shear γ = 0.1 | | | |
| engineering shear | 0.100000 | 0.100000 | 0.100000 |
| equivalent shear | 0.057735 | **0.057831** | 0.057831 |

The rotation case is the sharp one: the Green-Lagrange strain is *identically* zero for any
rigid motion, so the 1.5 % of apparent contraction the incremental measure reports simply
does not arise. It is not a better approximation, it is a property of the measure.

On `shear-block` the two differ and the finite one is exact. For engineering shear γ the
closed form is `E_xx = 0`, `E_yy = γ²/2`, `2E_xy = γ`; putting γ = 0.1 through the same
equivalent-shear formula gives 0.057831. The γ²/2 is a real term the linear theory drops.

On the real dam-break case, measured against the velocity filters that were tried and
rejected (see `experiments/` on the `noise-reduction` branch):

| | peak kept | shear band kept | quiet-zone noise cut | scatter |
|---|---|---|---|---|
| incremental | 100 % | 100 % | 0 % | 0.03765 |
| **finite, same velocities** | 144 % | **85 %** | **21 %** | **0.02802** |
| incremental + a 3×3 average of the velocities | 63 % | 63 % | 21 % | 0.01726 |

It removes as much noise as a mild spatial filter and keeps the feature the filter destroys,
because it is not smoothing anything: it computes a more correct quantity from the same
velocities. The spurious volumetric strain falls by 84 %.

One figure there is observed and not explained: the peak comes out 44 % *higher*, and the
finite correction accounts for only about 4 % at these strains. The rest is most likely the
incremental sum losing deformation along a path that rotates. It should be isolated on a
synthetic case with large shear before anyone relies on it.

The cost is four more blocks in the `.POST.RES`, about 40 % more file, and a few per cent of
run time.

## 8. Still open

* **How the dry–saturated band of each soil is measured.** It is the most sensitive number
  of the whole moisture method and currently comes from two clicks and four constants.
* **Measuring the drying branch** (without the incremental policy) against a test that
  actually dries.
* **The SWIR flow with two real cameras.** The registration machinery is in place and was
  validated by recovering a homography from the data, but no pair of genuinely simultaneous
  visible and SWIR images has been available: in the collections checked the SWIR files were
  byte-for-byte copies of the visible ones.
* **Particle mass from `S_DENSITY` and `POROSITY`**, which are read but not yet used (the
  mass is 1).
* **A `.PAR` with named keys**, so the order stops mattering at all. See
  [`ARCHITECTURE.md`](ARCHITECTURE.md), P2.
