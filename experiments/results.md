# Results

## Phase 1 — spatial smoothing of the velocity

Every filter is applied to the velocity field before it is differentiated, through a source
that wraps the PIVlab one, so the solver is untouched. Bias is measured on the two cases
whose answer is known; noise on the real dam-break case.

| Filter | γxy (truth 0.100000) | rotation (truth 50.000) | strain noise | quiet zone |
|---|---|---|---|---|
| none | 0.100000 | 49.997 | 0.03885 | 0.01096 |
| median 3×3 | 0.087399 | 43.201 | 0.02361 | 0.00839 |
| average 3×3 | 0.087399 | 43.201 | 0.01920 | 0.00537 |
| average 3×3 twice | 0.081099 | 40.077 | 0.01440 | 0.00400 |
| gaussian σ=0.7 | 0.092834 | 46.039 | 0.02265 | 0.00662 |
| gaussian σ=1.0 | 0.086914 | 42.509 | 0.01689 | 0.00485 |
| gaussian σ=1.5 | 0.077148 | 35.815 | 0.01195 | 0.00327 |
| savitzky-golay 5×5 | 0.086914 | 42.519 | 0.02265 | 0.00647 |
| **edge: average 3×3** | **0.100000** | **49.997** | **0.01940** | 0.00546 |
| **edge: average twice** | 0.099976 | 49.992 | **0.01468** | 0.00413 |
| **edge: gaussian σ=1.0** | 0.099975 | 49.986 | 0.01716 | 0.00494 |
| **edge: gaussian σ=1.5** | 0.099960 | 49.967 | **0.01205** | 0.00348 |
| **edge: savitzky-golay 5×5** | 0.099980 | 49.993 | 0.02299 | 0.00667 |

*strain noise*: scatter of the equivalent strain between neighbouring particles in the real
case. *quiet zone*: mean strain of the quarter that deformed least, which should be zero.

### What came out of it

**1. Smoothing naively costs more than it gives.** Every filter in the first block cuts the
noise by 39 % to 69 % and pays for it with a bias of 13 % to 28 % on an answer we know
exactly. A code that measures strain cannot be 20 % wrong so that the map looks tidier.

**2. The bias is entirely a boundary effect.** Dropping the particles near the edge of the
material makes it vanish:

| Filter | all | 1 cell in | 2 cells in | 3 cells in |
|---|---|---|---|---|
| none | 0.100000 | 0.100000 | 0.100000 | 0.100000 |
| average 3×3 | 0.092857 | 0.091667 | **0.100000** | **0.100000** |
| savitzky-golay 5×5 | 0.090739 | 0.091313 | 0.096780 | **0.100000** |
| gaussian σ=1.5 | 0.080900 | 0.082350 | 0.088399 | 0.092488 |

The filters are exact in the interior, whatever their width, and the damaged band is exactly
as wide as the kernel reaches. At the edge a symmetric kernel has no neighbours on one side,
so renormalising over the side it does have leans the result inwards and flattens the very
gradient being measured. It is the opposite of what `ICONTOUR` does.

**3. Told what lies beyond the edge, the same filters stop biasing.** Continuing the field
outwards by linear extrapolation before filtering — past the gaps *and* past the border of
the grid, which is where the material ends in `shear-block` — leaves the known answers
intact and keeps the noise reduction:

| | naive | edge-aware | noise cut |
|---|---|---|---|
| average 3×3 | −12.6 % | **0.000 %** | 50 % |
| average twice | −18.9 % | −0.02 % | 62 % |
| gaussian σ=1.5 | −22.9 % | **−0.04 %** | 69 % |

A factor of 575 less bias, for the same noise reduction.

### What this does not yet prove

Both cases with a known answer have a **linear** velocity field (simple shear, and rigid
rotation is linear in x and y as well). A symmetric kernel reproduces a linear field exactly,
so part of why the edge-aware filters look perfect is that there is no curvature to flatten.
A real shear band is not linear, and that is where a wide kernel would cost something real.

The next thing to measure is therefore a field with curvature: the `horizontal_*` synthetic
cases, whose strain changes from one column of cells to the next, and the width of a real
shear band before and after filtering.

### Candidate to take further

`edge: average 3×3` is the conservative one — exactly unbiased on both cases, half the noise,
one line of kernel. `edge: gaussian σ=1.5` buys another 19 points of noise reduction for
0.04 % of bias, which is likely worth it but should be judged on the curvature test first.

## Phase 2 — not smoothing, but differentiating once

Not started.

## Phase 3 — vector validation upstream

Not started.

## Phase 4 — moisture

Not started.

---

## A note on this machine

`numpy` here crashes (`0xc06d007f`) on **any** matrix multiplication, even 3×3, and on
anything that reaches LAPACK (`lstsq`, `polyfit`, `corrcoef`). The environment carries three
BLAS libraries and four OpenMP runtimes at once, which is the usual cause; a `scipy`
installed with pip into a conda environment brings its own OpenBLAS and is the likely
culprit. PIV-NP itself is unaffected — it uses no linear algebra, which is why its 347 tests
never caught it — but the code in this folder works around it on purpose.
