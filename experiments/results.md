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

## Phase 2 — not smoothing, but differentiating once

Next. The strain from the deformation gradient of the accumulated displacement.

## Phase 3 — vector validation upstream

Not started. Worth checking whether the PIVlab exports carry the normalised median test
before anything else is attempted.

## Phase 4 — moisture

Not started.

---

## A note on this machine

`numpy` here crashes (`0xc06d007f`) on **any** matrix multiplication, even 3×3, and on
anything reaching LAPACK. Uninstalling the pip `scipy` removed one of the three BLAS
libraries but not the fault; the environment still has conda's MKL shims and Intel MKL
together. Switching the environment to OpenBLAS is the usual way out:

```
conda install -p <env> "libblas=*=*openblas"
```

PIV-NP itself is unaffected — it uses no linear algebra, which is why its 347 tests never
caught it — but the code in this folder works around it on purpose.
