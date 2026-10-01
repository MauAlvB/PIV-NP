# Noise reduction: work plan

Working notes of the `noise-reduction` branch. Nothing here is part of PIV-NP: these are
experiments that measure whether an idea is worth putting in it.

## What we are trying to fix

Measured on the dam-break case, step 10:

| | |
|---|---|
| Typical speed of the field | 0.086 mm/s |
| Point-to-neighbour scatter of the velocity | **0.045 mm/s, 52 % of the signal** |
| Amplification when differentiating (`dt`/cell) | **231** |
| Resulting strain noise, **per step** | **0.0104** |

The last figure is the problem. The band of "practically intact" material in that case has
an accumulated `Equi_strain` below 0.01, which is *less than one step of noise*. So in the
quiet zones the strain map is mostly noise, and the same goes for everything derived from
it: `Vol_strain`, `Vorticity`, and above all `Vorticity_num`, which divides one noisy number
by another.

The noise is already in the velocity field at short wavelength. It is not produced by
PIV-NP, which only differentiates what it is given — but differentiating is what turns 52 %
of scatter into a strain of the same size as the real one.

## The rule of the exercise

**A filter is judged by whether it moves an answer we know, not by how the map looks.**

That is the asset this repository has and most PIV work does not:

| Case | Known answer |
|---|---|
| `shear-block` | γxy = 0.100000, equivalent strain 0.057735, vorticity −0.0500, Wm = 1 |
| `rotation_1P` | rotation 50.000°, deviatoric shear exactly 0 |
| synthetic cases | strain fields computed by finite differences, independently |
| regression suite | byte for byte against the 2024 Fortran |
| dam-break SWIR | reproduces a published moisture analysis |

A filter that cuts the noise by 80 % but shifts γxy from 0.100 to 0.098 has broken the
measurement. Both numbers go in the same table, always.

## Plan

### Phase 0 — the yardstick (prerequisite)

A harness that, for a given filter, reports **bias** and **noise** side by side:

* bias: how much the known answers move, on `shear-block` and `rotation_1P`
* noise: the point-to-neighbour scatter of the strain on the real case, and the spurious
  strain left in the zone that should be still

### Phase 1 — spatial smoothing of the velocity (cheap, big effect)

Candidates, all applied to the velocity **before** it is differentiated, and all of them
having to respect the gaps (`NaN`) and the edge of the material:

| # | Filter | Why it is a candidate |
|---|---|---|
| 1.1 | 3×3 moving average, one and two passes | the floor: a plain average already cut the noise by 79 % |
| 1.2 | Gaussian, σ in grid cells | the same with a softer kernel and a tunable width |
| 1.3 | Savitzky–Golay (local polynomial) | gives the *derivative* smoothed, which is what is actually wanted |
| 1.4 | Penalised least squares (Garcia 2010) | the standard in PIV post-processing; picks its own smoothing by cross-validation and handles gaps natively |

### Phase 2 — not smoothing, but differentiating once

Computing the strain from the **deformation gradient of the accumulated displacement**
instead of summing per-step increments. Two gains at once: it removes the rigid-rotation
artifact (the 1.5 % of contraction that never happened) and it stops accumulating the
differentiation noise of every step. Already on the roadmap in `ARCHITECTURE.md`.

### Phase 3 — the input, not the output

Checking whether the PIVlab exports carry vector validation (the normalised median test).
If they do not, that is the cheapest correction available and it belongs upstream of
everything here.

### Phase 4 — moisture

A separate diagnosis: there the dominant error is **systematic, not random** (band width
0.105, badly chosen reference 0.152), so no filter reaches it. What is worth trying is
illumination correction and a multi-channel calibration. Last, because the first three
phases pay more.

## Results

Filled in as the experiments run. See `results.md`.
