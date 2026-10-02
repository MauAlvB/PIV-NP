# Results

## Phase 1 — spatial smoothing of the velocity

Verdict first: **do not put any of these in PIV-NP.** The experiment found one thing worth
keeping (how a filter must treat the edge) and one thing worth abandoning (the filtering
itself). What it mainly bought is a reason to go after phase 2 instead.

Every filter is applied to the velocity field before it is differentiated, through a source
that wraps the PIVlab one, so the solver is untouched.

### 1. Bias, on the cases whose answer is known

| Filter | γxy (truth 0.100000) | rotation (truth 50.000) | strain noise | quiet zone |
|---|---|---|---|---|
| none | 0.100000 | 49.997 | 0.03885 | 0.01096 |
| median 3×3 | 0.087399 | 43.201 | 0.02361 | 0.01302 |
| average 3×3 | 0.087399 | 43.201 | 0.01920 | 0.00863 |
| gaussian σ=1.0 | 0.086914 | 42.509 | 0.01689 | 0.00755 |
| gaussian σ=1.5 | 0.077148 | 35.815 | 0.01195 | 0.00680 |
| savitzky-golay 5×5 | 0.086914 | 42.519 | 0.02265 | 0.00995 |
| **edge: median 3×3** | **0.100000** | **49.997** | 0.02380 | 0.01308 |
| **edge: median 5×5** | **0.100000** | **49.997** | 0.01551 | 0.01001 |
| **edge: average 3×3** | **0.100000** | **49.997** | 0.01940 | 0.00871 |
| **edge: gaussian σ=1.5** | 0.099960 | 49.967 | 0.01205 | 0.00748 |

**All the bias came from the boundary.** Dropping the particles near the edge of the
material makes it vanish entirely:

| Filter | all | 1 cell in | 2 cells in | 3 cells in |
|---|---|---|---|---|
| average 3×3 | 0.092857 | 0.091667 | **0.100000** | **0.100000** |
| savitzky-golay 5×5 | 0.090739 | 0.091313 | 0.096780 | **0.100000** |
| gaussian σ=1.5 | 0.080900 | 0.082350 | 0.088399 | 0.092488 |

In the interior the filters are exact whatever their width, and the damaged band is exactly
as wide as the kernel reaches. At the edge a symmetric kernel has no neighbours on one side,
so renormalising over the side it does have leans the result inwards and flattens the very
gradient being measured — the opposite of what `ICONTOUR` does.

Continuing the field outwards by linear extrapolation before filtering removes the bias
completely. **Both kinds of edge have to be handled**: the gaps inside the grid, and the
border of the grid itself, which is where the material ends in `shear-block`. Missing the
second is why the first attempt changed nothing at all.

**This part is a real finding and it survives.** Any future filtering anywhere in PIV-NP has
to extend the field past the boundary before it touches it, or it will quietly cost between
13 % and 28 % of the strain.

### 2. Curvature, where the synthetic cases stop being easy

Both cases with a known answer have a **linear** velocity field, which a symmetric kernel
reproduces exactly. `horizontal_1P` has a strongly curved one (0, 0.1, 0.3, 0.9, 2.4 cm/s
from column to column), and there the real cost appears:

| Filter | peak kept | contrast kept |
|---|---|---|
| edge: median 3×3, 5×5 | **100 %** | **100 %** |
| edge: savitzky-golay 5×5 | 89 % | 66 % |
| edge: average 3×3 | 84 % | 53 % |
| edge: gaussian σ=1.5 | 77 % | 43 % |
| gaussian σ=1.5, naive | 44 % | 10 % |

The median looked like the answer: unbiased, and exact on the curved case too.

### 3. The real case, which disagrees

| Filter | peak kept | strongest 5 % kept | quiet-zone noise cut |
|---|---|---|---|
| edge: savitzky-golay 5×5 | 75 % | 70 % | 6 % |
| edge: median 3×3 | 74 % | 70 % | **−19 %** |
| gaussian σ=0.7 | 72 % | 73 % | 26 % |
| edge: average 3×3 | 63 % | 63 % | 21 % |
| edge: median 5×5 | 52 % | 48 % | 9 % |
| edge: gaussian σ=1.5 | 40 % | 45 % | 32 % |

**Every one of them removes more signal than noise.** The best trade on offer keeps 73 % of
the shear band to remove 26 % of the quiet-zone noise, and the shear band is visibly a
coherent diagonal feature, not scatter. The median, which was exact on the synthetic cases,
keeps only half the peak here and makes the quiet zone *worse*.

The reason the median falls apart is width: it removes anything narrower than about half its
window, the synthetic columns are many cells wide, and a real shear band is one or two.

### A metric of ours that was wrong

The first version of `measure.py` took the quiet zone as the lowest quarter **of each
filtered field**, so every filter was compared against its own quartile. That makes "it
scaled the whole map down" look identical to "it removed noise", and it flattered the
filters badly: `edge: median 5×5` appeared to cut the quiet zone by 48 % when, measured on
the same particles throughout, it cuts 9 %. The mask is now fixed from the unfiltered run.

### What phase 1 concludes

* **Keep**: the edge treatment. It is a prerequisite for any filtering, and without it a
  filter silently costs 13–28 % of the strain.
* **Drop**: spatial smoothing of the velocity as a way to reduce noise. On real data the
  noise is at the same spatial scale as the signal, so a local kernel cannot separate them;
  it only scales the map down. The synthetic cases, being linear or piecewise-constant over
  wide regions, cannot show this and made the filters look far better than they are.
* **Why phase 2 is the better bet**: the redundancy in this data is in **time**, not in
  space. The analysis has 20 steps of the same material. Differentiating once over the
  accumulated displacement, instead of 20 times and adding up, uses that redundancy without
  touching the spatial resolution — and it removes the rigid-rotation artifact at the same
  time.

## Phase 3 — the images, and why none of this was ever going to work

This phase was supposed to ask whether preprocessing the photographs would give PIV a better
field. It answered that, and then answered a larger question that makes the rest of the
investigation moot.

A small PIV was written for it (`piv.py`): window cross-correlation by FFT with a Gaussian
three-point sub-pixel fit. Not a competitor to PIVlab — no multi-pass, no window deformation
— but the *same* algorithm applied to two versions of an image, so that the difference
between them is the preprocessing and nothing else. It was checked first against known
shifts of a real photograph and recovers them to **0.07–0.14 px**, the normal range for a
sub-pixel PIV. That number matters later.

### Preprocessing: already done, and more makes it worse

The `Masked_*.jpg` the analysis ran on carry **52 % more local contrast** than the camera
originals, so they had already been enhanced. Running the same PIV on both:

| Image | scatter | outliers |
|---|---|---|
| original, untouched | 0.5972 px | 3.9 % |
| **as delivered (enhanced)** | **0.4716 px** | **1.9 %** |
| original + local contrast | 1.7518 px | 7.7 % |
| original + high-pass σ=8 | 1.4490 px | 6.4 % |
| original + anti-blocking | 0.5754 px | 3.9 % |

The enhancement already applied is worth 21 % of the scatter and half the outliers. Every
further step tried makes things worse, some by a factor of three.

The reason is worth keeping: these are photographs of **soil texture**, not of tracer
particles. Soil carries strong large-scale texture — grains, shadows — that correlates well,
and a high-pass throws exactly that away and keeps the fine detail, which on a 177 KB JPEG
is mostly compression noise. PIV preprocessing recipes are written for particle images;
borrowing them for soil is counterproductive.

One measurable defect remains, and it is small: the JPEG 8×8 grid is 37 % stronger than the
detail around it, which is structure locked to the pixels that does not move with the soil.
Damping it buys about 4 % of scatter on the originals and nothing on the delivered images.
For future tests the lesson is to mask and enhance **once**, from the camera file, and save
lossless.

### The measurement that ends the investigation

| | |
|---|---|
| Image scale | 0.54 mm/px |
| **Typical displacement between two photographs** | **0.163 px** |
| Scatter between neighbouring vectors | 0.129 px |
| **Signal to noise** | **1.3 : 1** |
| What sub-pixel PIV can resolve | ~0.1 px |

**The soil moves about a sixth of a pixel per frame, and PIV resolves about a tenth.** The
scatter is not noise left behind by something: it is the resolution limit of the method, and
the signal sits barely above it.

That explains every result above at once. No filter helped because there was nothing to
separate. The quiet zones look like noise because they are moving *below* what can be
measured at all. Spatial smoothing removed signal because at that level the signal has no
spatial redundancy left to exploit.

### What to do instead

Not a better filter. A **bigger displacement**, which raises the signal while the 0.1 px
floor stays where it is:

1. **Correlate frames further apart.** Nothing in PIV-NP has to change: it is a choice when
   exporting from PIVlab, and the `.PAR` already checks that `DT` matches the interval the
   files were exported with. Going from every frame to every fifth would take the
   displacement from 0.16 px to about 0.8 px and the signal-to-noise from 1.3 to roughly 6.
   The classic quarter-rule limit for a 32 px window is 8 px, so there is a lot of headroom —
   the current data uses about 2 % of what the window could carry. The cost is temporal
   resolution, and the limit is that the soil must not deform so much between frames that
   the windows stop looking alike.
2. **Photograph at a finer scale.** More pixels per millimetre is more pixels of
   displacement for the same motion.
3. **Phase 2 still stands**, and for the same reason: differentiating once over twenty steps
   of accumulated displacement uses the signal twenty times over against a floor that does
   not grow.

The cheapest of the three by far is the first, and it is available on data already recorded.

## Phase 2 — not smoothing, but differentiating once

**This is the one worth putting in PIV-NP.**

The particles start on a regular lattice, so the accumulated displacement can be
differentiated with respect to the initial coordinates directly, giving the deformation
gradient `F = I + du/dX`. From it the Green-Lagrange strain `E = (F'F − I)/2`, the rotation
from the polar decomposition (closed form for 2×2, so no LAPACK), and the true area change
`det(F) − 1`. One differentiation over the whole analysis instead of one per step.

### On the cases whose answer is known

| | incremental | finite |
|---|---|---|
| **`rotation_1P`** — turns 50°, does not deform | | |
| equivalent shear | 0.005077 | **0.000000** |
| volumetric | −0.015230 | **0.000000** |
| rotation | −49.9975 | **−50.0000** |
| **`shear-block`** — simple shear γ = 0.1 | | |
| engineering shear | 0.100000 | **0.100000** |
| equivalent shear | 0.057735 | 0.057831 |
| volumetric | 0 | 0 |

The rotation artifact is gone outright: zero strain, zero volume change, and the rotation
exact to four decimals. That is not a better approximation, it is a property of the measure
— `E` is identically zero for any rigid motion.

The difference on `shear-block` is not an error either. For engineering shear γ the exact
Green-Lagrange strain is `E_xx = 0`, `E_yy = γ²/2`, `2E_xy = γ`. Putting γ = 0.1 through the
same equivalent-shear formula gives **0.057831**, which is what came out. The γ²/2 is a real
term the linear theory drops.

### On the real case, against the filters of phase 1

| method | peak kept | band kept | quiet-zone noise cut | scatter |
|---|---|---|---|---|
| incremental (as shipped) | 100 % | 100 % | 0 % | 0.03765 |
| **finite, same velocities** | **144 %** | **85 %** | **21 %** | **0.02802** |
| incremental + edge: average 3×3 | 63 % | 63 % | 21 % | 0.01726 |
| incremental + edge: median 5×5 | 52 % | 48 % | 9 % | 0.01362 |
| incremental + edge: gaussian σ=1.5 | 40 % | 45 % | 32 % | 0.00944 |

It removes the same noise as the 3×3 average and keeps the feature the average destroys. It
is not trading signal for noise, because it is not smoothing anything: it computes a more
correct quantity from the same velocities. The spurious volumetric strain drops by 84 % as
well, from −0.00286 to −0.00047.

**One thing is observed and not explained.** The peak comes out 44 % *higher*, not lower.
The finite correction itself only accounts for about 4 % at these strains — checked by hand
on simple shear — so the rest is most likely the incremental sum losing deformation along a
path that rotates, which is a known effect. It should be isolated on a synthetic case with
large shear before anyone leans on it.

### What it would take in PIV-NP

The arithmetic is small and the inputs already exist: `initial_position` and `displacement`
are kept per particle, and the initial lattice is regular by construction. What needs
deciding is not the arithmetic but the contract:

* whether the finite measure **replaces** the incremental one or is published beside it.
  Replacing changes every existing result, which is a break; publishing both costs two more
  blocks and lets them be compared on real work before anything is retired.
* `--legacy-compat` has to keep the incremental one exactly, as it does for the other
  additions.
* lost particles and gaps must stay out of the differences, or a particle that stopped moving
  when it left the mesh drags its neighbours. The prototype already does this.

## Phase 4 — moisture

Not started.

---

## Phase 5 — PIV from the photographs, inside PIV-NP

Verdict first: **it works, within the 10 % agreed, and it is an entry point rather than a
replacement for PIVlab.** The displacements agree to 8 %, the amount of strain to 7 %, and
the two agree on *where* the strain is only in broad strokes. Which is the honest version of
what was asked for: something a newcomer can run from photographs alone, and then decide
whether to go and learn PIVlab.

Run it with `python experiments/end_to_end_piv.py` (needs `PIVNP_TEST_IMAGES`; add
`--fresh` to ignore the cached runs).

### What was compared, after two wrong attempts

Both wrong attempts are worth recording, because each made the implementation look worse
than it was and nearly sent the work off after the wrong thing.

1. **Accumulating with `nan_to_num`.** The first comparison added the twenty step fields
   together after turning every `NaN` into zero. About 45 % of each field has no data, in
   *different places* in the two, so this quietly compared a sum of twenty measurements
   against a sum of eleven. It reported 19 % and the real figure was never that bad. (The
   script that did it is not in the repository; it would print a believable wrong number.)
2. **Judging the strain particle by particle.** The second comparison asked for the strain
   of each particle to match within 10 %. The strain is a difference between neighbouring
   vectors, and at 0.16 px of movement per step it is a small difference between two noisy
   numbers; two runs of PIVlab with different validation settings would not agree either.
   It reported 105 % while the median strain of the whole field was 2 % apart.

What is compared now: the **displacement particle by particle**, which is a fair question
and the strict one, and the **strain over the field** — how much of it there is, and whether
it sits in the same places.

### The result

| | from PIVlab | from the photographs | difference |
|---|---|---|---|
| displacement x, per particle | 2.437 mm | 2.243 mm | 7.6 % |
| displacement y, per particle | 3.105 mm | 2.959 mm | 5.8 % |
| equivalent strain, field | 0.08325 | 0.08187 | 1.7 % |
| shear strain, field | 0.06075 | 0.06416 | 5.6 % |
| volumetric strain, field | 0.05874 | 0.06290 | 7.1 % |

Particle by particle, the two put 60 % of particles within 10 % of each other and 91 %
within 25 %; the median particle ends 0.29 mm from where PIVlab put it, having travelled
about 3.9 mm.

### Where it does *not* agree, and what that means

The two agree on how much strain there is far better than on where it is. Grouping the
particles into patches and correlating:

| patch | particles averaged | equivalent strain | shear |
|---|---|---|---|
| none | 1 | 0.53 | 0.49 |
| 10 mm | 5 | 0.58 | 0.42 |
| 20 mm | 16 | 0.60 | 0.36 |
| 40 mm | 48 | 0.61 | 0.34 |
| 80 mm | 140 | 0.83 | 0.76 |

If the difference were point-to-point noise it would average away as the patches grow and
the correlation would climb steadily. **It does not.** It is flat from 1 to 48 particles and
only climbs at 140, which says the disagreement is *spatially structured* at scales up to
about 40 mm — the two measure genuinely different fields at that scale and agree on the
coarse pattern. Not yet run down; candidates are the different coverage (42 % of our field
missing against PIVlab's 49 %, in different places) and whatever PIVlab's own validation
does that this does not.

So: good enough to see where a slope failed and roughly how hard, not good enough to publish
a strain pattern from without checking it against PIVlab.

### Two bugs found by the tests, not by the comparison

* **The window pushed past the edge.** The second pass takes each window of the second image
  from where the first pass said the soil went. Near the border of the photograph that
  offset points outside the image, and the window is pulled back in — but the code added
  back the offset it had *asked for*, so a true −1.7 px came back as −3.7. It hid itself:
  the vectors were wrong enough that the outlier test discarded them, and a discarded row
  looks like the edge of the material rather than like a mistake. Fixing it kept **100 % of
  the windows of a synthetic pair instead of 93 %** — the bug had been eating the border of
  every field. It changes nothing on the dam-break case, which moves 0.16 px per step: the
  integer offset is zero almost everywhere, so the clipping never triggers. It would bite
  any faster test.
* **Smoothing that flattened the gradient at the edge.** Phase 1 had already measured that a
  kernel renormalised over the neighbours that exist leans inwards and flattens the very
  gradient being measured, and that continuing the field outwards first removes the reason
  for it. The first version of the PIV smoothing renormalised anyway. On a field with a
  known gradient and 42 % of its points missing:

  | | typical point | worst point | full grid |
  |---|---|---|---|
  | continued outwards | 0.012 | 0.147 | exact |
  | renormalised | 0.060 | 0.196 | 0.154 |

  Five times better at the typical point, so the shipping version continues the field
  outwards. **The test that should have caught this checked only `[2:-2, 2:-2]`** — it
  passed over a smoother that was wrong at every edge, which is the lesson worth keeping.

### Why the field is smoothed at all

The raw field is 1.45× rougher than PIVlab's between neighbouring windows (0.082 px against
0.055), and the strain gap was 1.38× — essentially all of it. PIVlab's exported field is
smoothed by its own post-processing, so comparing a raw field to it compares two different
things. Sweeping the width:

| σ (grid points) | roughness vs PIVlab | distance to PIVlab's vectors |
|---|---|---|
| 0 | 1.46× | 0.119 px |
| 0.4 | 1.29× | 0.109 |
| **0.6** | **0.84×** | **0.092** |
| 1.0 | 0.46× | 0.101 |
| 2.0 | 0.18× | 0.130 |

σ = 0.6 both matches PIVlab's roughness and *minimises* the distance to its vectors, which
is the independent evidence that what is removed is noise and not soil — had it only matched
the roughness, this would just be fitting to PIVlab. It is the `SMOOTH` key of the `.PIV`
file and `SMOOTH = 0` turns it off.

---

## Phase 6 — why the built-in PIV and PIVlab disagree on the strain pattern

Verdict first: **the question was half wrong, and the half that survives has a well-supported
answer that has not been implemented.** Phase 5 reported that the two agree on how much
strain there is but only broadly on where, and read a correlation of 0.5 per particle as
evidence of something structured. That reading had an untested assumption in it — that 0.5 is
low — and measuring the thing it was compared against changes the conclusion.

Run it with `python experiments/piv_reproducibility.py` (needs `PIVNP_TEST_IMAGES`).

### 1. The baseline that was missing: how well does a measurement agree with itself?

The twenty steps split into the odd ones and the even ones, each run as a complete ten-step
analysis. The two sample the same physical event, interleaved, so they should give the same
strain pattern at about half the magnitude. How well they agree is the ceiling; nothing can
beat it.

| correlation of the equivalent strain | per particle | 20 mm patches | 40 mm | 80 mm |
|---|---|---|---|---|
| PIVlab against itself | 0.569 | 0.889 | 0.952 | 0.940 |
| **the built-in PIV against itself** | **0.731** | 0.884 | 0.960 | **0.978** |
| the built-in PIV against PIVlab | 0.528 | 0.482 | 0.608 | 0.828 |

Two things follow, and the first overturns phase 5:

* **0.52 per particle is not low.** It is the reproducibility floor: PIVlab manages 0.569
  against its own data. The strain of a single particle on this test is simply not a
  reproducible quantity, for anybody, because it is a difference between neighbouring
  vectors at 0.16 px of movement per step. Reading 0.52 as a defect was wrong.
* **Our field is the more reproducible of the two**, not the less — 0.73 against 0.57 per
  particle, 0.98 against 0.94 over 80 mm patches.

What survives is narrower and still real: our agreement with PIVlab climbs more slowly with
averaging than PIVlab's agreement with itself (0.83 against 0.94 at 80 mm). Two internally
consistent measurements that differ from each other differ *systematically*. The claim in
phase 5 that the disagreement is structured "up to about 40 mm" was also too strong: the
per-step difference field decorrelates within 2 to 4 grid points, which is 4 to 9 mm.

### 2. Four candidates, measured and eliminated

| candidate | how it was tested | verdict |
|---|---|---|
| peak locking | correlate the disagreement with the fractional part of the displacement | **no**: r = 0.00 to 0.03 over five steps |
| different coverage | we measure 13.5 % of the grid PIVlab rejects; reject it too and compare | **no**, and it makes things *worse* |
| smoothing length-scale | sweep our smoothing and watch the shear-rate pattern agreement | **no**: it peaks at our own default |
| grid registration | compare the two sets of grid coordinates directly | **no**, but see below |

The coverage test is worth keeping because it nearly convinced. The correlation peak ratio
separates the two populations cleanly — median 1.85 where both measure, **1.14 where only we
do** — so PIVlab is plainly rejecting the windows without a clear peak and we are keeping
them. Thresholding at 1.10 matches its coverage almost exactly (50.6 % against 51.0 %). But
rejecting them *degrades* the agreement, from 0.61 to 0.49 at 40 mm patches and 0.83 to 0.73
at 80 mm, and on the synthetic case it doubles the scatter of a shear that should be
constant. Removing a vector leaves a hole, the strain of its neighbours is then taken across
that hole, and an interpolated vector is worse than a measured one with noise in it. So the
extra coverage is not the cause, and `MIN_PEAK_RATIO` was implemented, measured, and **taken
back out**.

The smoothing sweep is the other informative negative. Correlating the shear rate of the two
fields: 0.596 with no smoothing, **0.619 at σ = 0.6**, 0.535 at 1.0, 0.376 at 1.5, 0.185 at
3.0. The default was chosen in phase 5 for a different reason and turns out to maximise this
too, so no filter can close the gap.

Grid registration found a real but small flaw: our grid sits **0.26 mm** from PIVlab's in
both axes, a constant offset of about half a pixel, on a grid step of 4.33 mm. It comes from
taking the window centre at `(window-1)/2` where PIVlab takes `window/2`. It cannot explain
the disagreement — a translated field has the same strain pattern — but it means our
coordinates are half a pixel from where they claim to be, which matters for registering
against the moisture images. Worth fixing on its own account.

### 3. What the evidence does point at, and the measurement that shows it

The second pass offsets each window of the second image by what the first pass found, and
**rounds that offset to a whole pixel**. On a test moving 0.16 px per step the rounded offset
is zero nearly everywhere, so the second pass is the first one again with a smaller window.
Measured: at 0.16 px, one pass and two passes give *identical* numbers — bias −0.0088 px,
scatter 0.0394 px.

PIVlab's multi-pass does shift by the fraction, interpolating the image. So the two are not
the same estimator, which fits every observation above: a systematically different field,
each internally consistent, differing in a way no filter and no validation threshold can
reconcile.

A prototype that shifts by the fraction, against displacements known exactly:

| true displacement | bias, rounded | scatter, rounded | bias, fractional | scatter, fractional |
|---|---|---|---|---|
| 0.05 px | −0.0039 | 0.0123 | −0.0003 | **0.0026** |
| 0.10 | −0.0066 | 0.0250 | −0.0003 | **0.0051** |
| 0.16 | −0.0088 | 0.0394 | −0.0001 | **0.0077** |
| 0.25 | −0.0095 | 0.0623 | 0.0003 | **0.0110** |
| 0.35 | −0.0096 | 0.0842 | 0.0006 | **0.0135** |
| 0.50 | 0.0054 | 0.1090 | 0.0219 | 0.1000 |
| 0.65 | 0.0094 | 0.0853 | −0.0007 | **0.0163** |
| 1.30 | −0.0094 | 0.0743 | 0.0006 | **0.0129** |
| 2.30 | −0.0094 | 0.0743 | 0.0005 | **0.0125** |
| 4.70 | 0.0094 | 0.0747 | −0.0006 | **0.0134** |

**Bias about fifteen times smaller, scatter about five times smaller**, at every displacement
except the half-pixel case where the two tie. This is the largest improvement found anywhere
in this folder, and it lands exactly where this project's data lives: slow tests whose
displacement is a fraction of a pixel.

Two honest qualifications. The prototype applies one shift to the whole image, which is only
valid because the synthetic field is uniform; a real implementation needs the fraction per
window, which is more code and more cost, and on a field with real shear the gain will be
smaller. And the prototype had two bugs of its own before it gave this — the image was
shifted the wrong way, doubling the displacement, and the final outlier rejection was missing
so the scatter was 0.65 px and meaningless. The first table it produced argued confidently
against the idea.

### 4. Where this leaves the open question

Answered enough to act on, not closed. The remaining disagreement is consistent with the two
being different estimators, and the way to find out is to implement the fractional offset and
re-run the comparison. If agreement moves towards PIVlab's own reproducibility, that was the
cause.

Not done here, because it changes every number the built-in PIV publishes and wants its own
round of validation.

---

## Phase 7 — reading the windows between pixels, implemented and measured

Verdict first: **implemented, better by every measure of a single step, and off by default,
because on the one case with a known accumulated answer it makes the answer worse.** Phase 6
proposed this as the leading improvement and the leading explanation of the disagreement with
PIVlab. Half of that held up.

`SUBPIXEL_OFFSET = 1` in the `.PIV`. `python experiments/compare_fields.py` draws the
comparison.

### What it does to a single step, where it is unambiguous

Known shifts, 16 px windows, two passes:

| true | bias, rounded | scatter, rounded | bias, between pixels | scatter |
|---|---|---|---|---|
| 0.05 px | −0.0039 | 0.0123 | −0.0005 | **0.0023** |
| 0.16 | −0.0088 | 0.0394 | −0.0011 | **0.0072** |
| 0.35 | −0.0096 | 0.0842 | −0.0019 | **0.0163** |
| 0.50 | 0.0054 | 0.1090 | 0.0009 | **0.0214** |
| 2.30 | −0.0094 | 0.0743 | −0.0017 | **0.0176** |

And peak locking, the error that swings with where the displacement falls between two pixels,
is gone:

| shift | rounded | between pixels |
|---|---|---|
| 1.00 | 0.0000 | 0.0225 |
| 1.25 | 0.0784 | 0.0254 |
| 1.50 | **0.1439** | **0.0211** |
| 1.75 | 0.0804 | 0.0182 |
| 2.00 | 0.0000 | 0.0226 |

Rounding is exact at a whole pixel and worst at a half; interpolating is flat. Worst case
across the fractions, 0.144 against 0.027.

### Why it is off anyway

On `examples/piv-from-images`, where the shear was put in and is known:

| | shear (truth 0.15625) | scatter | volumetric (truth 0) | vertical (truth 0) |
|---|---|---|---|---|
| rounded, as shipped | 0.15450 (**−1.1 %**) | 0.0626 | −0.00911 | 0.062 mm |
| between pixels | 0.14305 (**−8.4 %**) | 0.0587 | **−0.00189** | **0.034 mm** |

Five times less of the artifact that should not exist, and eight times more error in the
quantity that should. That is not a trade worth making by default on a code whose value is
that its numbers can be trusted.

What it is **not**: both estimators read the gradient of each single step to within 1.3 %
(cubic −1.3 %, −0.1 %, −0.8 % at steps 1, 5, 10; rounded +0.5 %, +0.9 %, 0.0 %). So the
estimator is fine and the loss happens over the ten steps of accumulation. Nor is it the
interpolation order: bilinear gave −8.6 % and cubic −8.4 %, which is why the code is cubic
(better on every other count) but also why that was not the cause. Nor is it concentrated at
the edge of the block — by depth band the error runs −12.5 %, −5.2 %, −12.2 %, +5.2 %,
scattered rather than concentrated.

Unresolved. The pinning test is
`test_reading_between_pixels_still_accumulates_worse_than_rounding`, which fails the day it
stops being true — which is the day the default should change.

### Against PIVlab on the real case, which muddies it further

| | median \|u\| | equivalent strain | within 10 % of PIVlab, displacement |
|---|---|---|---|
| PIVlab | 0.596 mm | 0.03046 | — |
| as shipped | 0.573 | 0.02957 | 57 % of particles |
| between pixels | 0.668 | 0.02871 | **71 %** |

So on the real test, reading between pixels agrees with PIVlab *better* — 71 % of particles
within 10 % against 57 % — while on the synthetic test it gets the known answer *worse*. Both
are true and they pull opposite ways. The figures show why it is not noise: the difference
maps are coherent regions, not speckle, and the two estimators bracket PIVlab — the rounded
one reads low over the body of the slope and the interpolating one reads high.

That PIVlab sits between them is consistent with phase 6's reading that the three are
different estimators rather than one being wrong, and it is a reason to trust the synthetic
case over the agreement with PIVlab when they disagree: only one of the two has an answer.

### Two bugs of my own, caught by tests rather than by reasoning

* The first version reserved a pixel of room for the interpolation by pulling the *window*
  in, which charged the last row and column of the grid a whole pixel of offset they never
  asked for, and that then cost the 0.12 px a one-pixel shift costs. The test on two
  identical photographs caught it. Clamping the neighbouring *sample* instead is correct.
* `floor` of a displacement of −1e-16 — which is what two identical photographs produce from
  round-off — is −1, with a fraction of 0.999…, so a window with nothing to resample got
  resampled. Snapping a fraction within 1e-9 of a whole pixel fixes it, and identical
  photographs measure exactly nothing again.

And one bug that was already in `main`, exposed by the change shifting memory around:
`gid_writer._finite` memoised on `id(p.displacement)`. An `id` is only unique while the
object lives, so the cache answered one analysis with another's results — an array of the
wrong length, between two cases in the same test session. It holds the array now.

---

## Phase 8 — Slope_RGB, a test that actually moves

Verdict first: **the built-in PIV handles it, and the interrogation window is what decides
whether the answer is worth anything.** Everything measured about it until now was measured
on the dam break, where the soil moves 0.16 px between photographs. This slope moves up to
19.5 px, which is the regime most tests are in, and it exposes something the dam break could
not.

`python experiments/slope_rgb.py`. The case is read-only; everything runs on a copy.

### The window, and a wrong answer that looked plausible

The case's grid steps by 50 px. A 64 px window at 0.219 overlap steps by 50, lands on
PIVlab's 60 x 35 points exactly, and runs without complaint. The result:

| | median \|u\| | equivalent strain | shear | rotation |
|---|---|---|---|---|
| PIVlab | 1261 mm | 0.357 | 0.144 | 9.8° |
| built-in PIV, 64 px window | 946 mm | **0.699** | **0.511** | **19.8°** |

Twice the strain and twice the rotation. Nothing in the run said anything was wrong.

What is wrong is the window. Measured over four steps, on the same grid, against PIVlab:

| window | overlap | measured | our roughness | PIVlab's | median gap |
|---|---|---|---|---|---|
| 64 | 0.219 | 64 % | **0.257** | 0.096 | 0.869 px |
| 100 | 0.500 | 74 % | 0.162 | 0.096 | 0.666 |
| 128 | 0.609 | 80 % | 0.141 | 0.096 | 0.560 |
| 160 | 0.688 | 86 % | 0.115 | 0.096 | 0.518 |
| **200** | 0.750 | **89 %** | **0.091** | 0.096 | 0.486 |

At 64 px the field is 2.7 times rougher than PIVlab's, and the strain is a difference between
neighbouring vectors, so a field 2.7 times rougher is a strain about twice as large. That is
the whole of it. PIVlab reaching a 50 px step at its usual half overlap would have used a
100 px window; 200 px at 0.75 overlap is what matches what it actually produced.

**The lesson generalises past this case.** The grid step is set by `WINDOW x (1 - OVERLAP)`,
so any step can be reached with a small window and little overlap or a large one with a lot.
Those two are not equivalent and nothing in the settings says so: the first is noise and the
second is a measurement. The window has to be chosen for how far and how unevenly the soil
moves, and the overlap then follows from the grid wanted. The README and the guide now say
this where the key is described.

At 200 px one pass is enough: a second changes the roughness from 0.096 to 0.095 and costs
5.7 times the time, because 19.5 px is a tenth of that window rather than a third of it.

### Two things this case forced into the code

* **The window no longer has to be a power of two.** Nothing in the transform needs it, and
  the rule was rejecting 200 px -- the only window that works here.
* **Correlating in bands**, since a 200 px window over a 4-megapixel photograph wants
  gigabytes of correlation planes at once. The result is bit-identical; the test forces the
  smallest band and requires equality rather than closeness.

And a third that is pure speed: the source was decoding every photograph twice, because step
*n* reads images *n* and *n+1* and step *n+1* reads *n+1* again. Remembering the last image
took five steps of this case from 74 s to 41 s.

---

## A note on this machine

Two more things on this machine, found in phase 7 and both worked around in
`compare_fields.py`:

* **matplotlib cannot save a figure.** It imports, it plots, and then it takes the
  interpreter down inside `savefig` with the same `0xc06d007f` — PNG and SVG alike, so it is
  the font and raster layer rather than the format. PIL writes PNGs perfectly well, so the
  figures are drawn with numpy and PIL.
* **Three analyses in one process is one too many.** Each holds several hundred megabytes of
  window stacks; freeing them between runs is not enough. A process per analysis is reliable.

`numpy` here crashes (`0xc06d007f`) on **any** matrix multiplication, even 3×3, and on
anything reaching LAPACK. Uninstalling the pip `scipy` removed one of the three BLAS
libraries but not the fault; the environment still has conda's MKL shims and Intel MKL
together. Switching the environment to OpenBLAS is the usual way out:

```
conda install -p <env> "libblas=*=*openblas"
```

PIV-NP itself is unaffected — it uses no linear algebra, which is why its 347 tests never
caught it — but the code in this folder works around it on purpose.
