# Where the project is

One page, kept current, for picking the work up again. It holds **no numbers of its own** on
purpose: everything measured lives in one authoritative place and this points at it, so there
is nothing here to go stale.

Last updated: 2026-10-01.

## Where each kind of information lives

| What | Where | Authoritative for |
|---|---|---|
| How to use the code | [`GUIDE.md`](GUIDE.md), [`README.md`](../README.md) | input files, options, results |
| What was checked and how much each decision moves | [`VALIDATION.md`](VALIDATION.md) | every measured claim about the shipped code |
| How it is laid out and how to extend it | [`ARCHITECTURE.md`](ARCHITECTURE.md) | the seams, the proposals, the roadmap |
| Visualisation | [`VISUALIZATION.md`](VISUALIZATION.md) | ParaView |
| Ideas that were tried, kept or rejected | [`experiments/`](../experiments/) | measurements of things not in the shipped path |
| Why a change was made | the commit messages | reasoning behind each step |

If a number appears in two of those, the one in `VALIDATION.md` wins for shipped behaviour
and the one in `experiments/results.md` wins for an idea that was only tried.

## Shipped, on `main`

The Python version is byte-for-byte identical to the compiled 2024 Fortran, which the
regression suite checks on every run. Beyond that: a single `.PAR` format with a converter
for the old dialects, moisture measured from the test photographs (`MOISTER=2`), the four
boundary corrections with the measured ranking, an interchangeable displacement source, two
runnable examples with their data, a getting-started guide, `Vorticity`, `Rotation`,
`Vorticity_num` and `Rot_angle`, and the finite strain from the deformation gradient
published beside the incremental one.

**PIV from the photographs**, so a case can be completed without PIVlab:
`pivnp <case> --source images`, configured by a `<case>.PIV` file. Inside the 10 % agreed
against PIVlab, and the numbers and — more to the point — the limits are in
[`VALIDATION.md`](VALIDATION.md#8-the-built-in-piv): it agrees on how much strain there is
far better than on where it is, and that is how it is described to users rather than being
left for them to find out. Its example, `examples/piv-from-images`, has photographs generated
with a known shear so that anyone can check the answer; the real test photographs are not in
the repository and are not ours to publish.

## What the noise investigation concluded

Three phases of measurement into where the noise comes from, with the conclusions in
[`experiments/results.md`](../experiments/results.md). Filtering the velocities and
preprocessing the images were both measured and **rejected**, because the soil moves about a
sixth of a pixel per frame while PIV resolves a tenth — there was no noise to separate from
signal. What survived and shipped is the finite strain from the deformation gradient, which
removes the rigid-rotation artifact outright and cuts the scatter without costing signal.

Kept here because a negative result is worth as much as a positive one: it is the reason not
to try these again.

## Decided, not yet done

* **Export PIVlab with frames further apart.** Costs nothing, needs no code change, and is
  the largest single improvement available to the data: the soil moves about a sixth of a
  pixel per frame while PIV resolves a tenth.
* **Two papers in parallel**, a software one and a methodological one, sharing the same
  validation data. The review protocol is written ([`review-protocol.md`](review-protocol.md))
  and the next move is the user's: nine citation exports, as set out in
  [`review-step-1.md`](review-step-1.md). Nothing else in the review can start until those
  arrive, so it is parked rather than in progress.

## Open questions

* The finite measure gives a peak 44 % higher than the incremental one and the finite
  correction only explains about 4 % of that. Isolate it on a synthetic case with large shear
  before relying on it.
* How the dry–saturated band of each soil is measured. It is the most sensitive number in the
  moisture method and currently comes from two clicks and four constants.
* **Why `SUBPIXEL_OFFSET` cannot be turned on yet.** Reading the windows between pixels is
  implemented and shipped, off by default. Per single step it is better by every measure —
  four to five times less bias and scatter, peak locking gone — and on the real case it
  agrees with PIVlab better than the default does, 71 % of particles within 10 % against
  57 %. But on the example whose answer is known it makes the accumulated shear worse, −8.4 %
  against −1.1 %, while cutting the artifacts that should be zero fivefold. Both estimators
  read each single step's gradient to within 1.3 %, so the loss is in the accumulation and
  nobody has found where; it is not the interpolation order and it is not the edge of the
  material. `test_reading_between_pixels_still_accumulates_worse_than_rounding` fails the day
  this is solved. See [`experiments/results.md`](../experiments/results.md), phase 7.
* The PIV grid sits 0.26 mm — half a pixel — from where its coordinates say, from taking the
  window centre at `(window-1)/2` instead of `window/2`. Harmless to the strain, not harmless
  to registering the grid against the moisture images.
* Whether the experimental data of the paper case can be published.
