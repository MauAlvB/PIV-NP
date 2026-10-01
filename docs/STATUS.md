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
| Work in progress, not shipped | `experiments/` on its branch | measurements of ideas being tried |
| Why a change was made | the commit messages | reasoning behind each step |

If a number appears in two of those, the one in `VALIDATION.md` wins for shipped behaviour
and the one in `experiments/results.md` wins for work in progress.

## Shipped, on `main`

The Python version is byte-for-byte identical to the compiled 2024 Fortran, which the
regression suite checks on every run. Beyond that: a single `.PAR` format with a converter
for the old dialects, moisture measured from the test photographs (`MOISTER=2`), the four
boundary corrections with the measured ranking, an interchangeable displacement source, a
runnable example with its data, a getting-started guide, `Vorticity`, `Rotation`,
`Vorticity_num` and `Rot_angle`, and the finite strain from the deformation gradient
published beside the incremental one.

## In progress, on `noise-reduction`

Three phases of measurement into where the noise comes from, with the conclusions in
[`experiments/results.md`](../experiments/results.md). The short version: filtering the
velocities and preprocessing the images were both measured and **rejected**, and the reason
is that the soil moves about a sixth of a pixel per frame while PIV resolves a tenth. What
survived is the finite strain from the deformation gradient, which removes the rigid-rotation
artifact outright and cuts the scatter without costing signal.

## Decided, not yet done

* **Export PIVlab with frames further apart.** Costs nothing, needs no code change, and is
  the largest single improvement available to the data: the soil moves about a sixth of a
  pixel per frame while PIV resolves a tenth.
* **A paper**, with the literature review built to PRISMA. Not started.

## Open questions

* The finite measure gives a peak 44 % higher than the incremental one and the finite
  correction only explains about 4 % of that. Isolate it on a synthetic case with large shear
  before relying on it.
* How the dry–saturated band of each soil is measured. It is the most sensitive number in the
  moisture method and currently comes from two clicks and four constants.
* Whether the experimental data of the paper case can be published.
