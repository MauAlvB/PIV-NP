# Example: PIV straight from the photographs, with no PIVlab

A complete case with **no PIVlab files at all**. PIV-NP measures the displacements from the
photographs itself and then analyses them as usual.

```bash
pivnp examples/piv-from-images --source images
```

It takes about a second. If you have never used PIV-NP before, run
[`examples/shear-block`](../shear-block) first — that one is the walkthrough in
[`docs/GUIDE.md`](../../docs/GUIDE.md), and this one assumes you have seen it.

## What it is

A block of soil 192 × 128 px, photographed eleven times, sheared a little between each pair:
the base stays put and the top slides 2 px further each step, 20 px in all. At 0.5 mm per
pixel that is 10 mm of slide over a 64 mm block.

The photographs are **generated**, not photographed, and that is the point. The shear was put
in by [`make_example.py`](make_example.py), so you can check what comes out instead of
trusting it. Real test photographs are hundreds of megabytes and usually belong to someone,
which is why no repository ships them.

## What has to come out

The true answer is a simple shear of γxy = 20 px / 128 px = **0.15625**, the same on every
particle, with no change of volume and no vertical movement.

| Result | Expected | Measured | |
|---|---|---|---|
| `Total_strain` xy (γxy) | 0.15625 | **0.15450** | −1.1 % |
| `Equi_strain` | 0.09021 (= γ/√3) | **0.09307** | +3.2 % |
| `Vol_strain` | 0 | **−0.00911** | see below |
| `Displacement` x, at the top | +10.00 mm | **+9.27 mm** | |
| `Displacement` x, at the base | 0 | **+0.75 mm** | |
| `Displacement` y | 0 | **0.06 mm** | |
| `Vorticity` | −0.015625 s⁻¹ | **−0.01495** | −4.3 % |
| particles published | | **1440** of 1500 | |

So the shear itself comes back to about one per cent, which is as good as this gets, and the
quantities that should be zero are small rather than zero.

**Why `Vol_strain` is not zero.** A simple shear changes no volume, so −0.9 % of apparent
contraction is entirely an artifact — the one described in
[`docs/VALIDATION.md`](../../docs/VALIDATION.md), where scatter in a rotating field shows up
as a fake volumetric contraction rather than as fake shear. It is worth seeing on a case
whose answer you know, because on a real test it is indistinguishable from soil actually
compacting.

The finite-strain fields published beside the incremental ones reduce it but do not remove
it here:

| | γxy or 2·Exy | area change | equivalent |
|---|---|---|---|
| truth | 0.15625 | 0 | 0.09021 |
| incremental | 0.15450 | −0.00911 | 0.09307 |
| finite, from `det F` | 0.15340 | −0.00760 | 0.09269 |

Closer on all three, and by much less than the difference between the two measures on a real
test. On this case the deformation is small enough that the two should agree, and they do —
which is itself worth knowing.

The sixty missing particles are the ones that leave the grid as the block shears, which is
what happens in a real test too.

## The files

| File | What it is |
|---|---|
| `PIV-NP.TXT` | the case name, `shearphotos` |
| `shearphotos.PAR` | what to analyse: the grid, the 10 steps, the options |
| `shearphotos.PIV` | **where the photographs are and how to measure them** |
| `images/shear_001..011.png` | the eleven photographs, 220 KiB in all |
| `make_example.py` | regenerates the photographs, and prints the answer |

There is no `pivlab/` folder, and that is the difference from every other case.

## The `.PIV` file, which is the only new thing

```
IMAGES = images/shear_{n:03d}.png   ! {n} is the image number; {n:03d} makes it 001, 002...
SCALE  = 0.0005                     ! metres per pixel
WINDOW  = 16                        ! the patch of soil each vector describes
OVERLAP = 0.5                       ! so there is a vector every 8 px
PASSES  = 2                         ! 32 px to find the movement, then 16 to pin it down
REGION = 28, 28, 248, 164           ! where the material is, so the grid is not on background
MASK_BELOW = 20                     ! the background here is black
```

`IMAGES` and `SCALE` have no defaults and nothing can guess them. `SCALE` is the one that
matters most: every result is in metres **because of that number**, so if it is wrong
everything is wrong by the same factor and nothing in the output will look odd.

Every key, with its default, is listed in the
[README](../../README.md#input-files), and the walkthrough is
[`docs/GUIDE.md`](../../docs/GUIDE.md).

## Things worth trying

**Turn the smoothing off.** Set `SMOOTH = 0` in the `.PIV` and run it again. The slide at the
top of the block goes from 9.274 mm to 9.275 — no change worth the word — while the spread of
γxy across the block grows from 0.063 to 0.096, half again as noisy. That is the whole
argument for smoothing in one experiment: the strain is a difference between neighbouring
vectors, so sub-pixel scatter that hardly shows in the displacement dominates it.

**Make the window smaller.** `WINDOW = 8` gives four times as many vectors and trusts each
one less, and the `.PAR` will have to change — run it and PIV-NP will print the `BLOCK 2`
line it wants.

**Break it on purpose.** Multiply `SCALE` by ten. Every displacement and every length comes
out ten times larger, the strains come out unchanged, and nothing warns you. Worth seeing
once.

## How this compares with PIVlab

On this case there is nothing to compare against — the answer is known, which is better. On
a real test the built-in PIV agrees with PIVlab to about 8 % on displacement and 7 % on how
much strain there is, but only broadly on *where* the strain is. It is a way in, not a
replacement; the measurements and the limits are in
[`docs/VALIDATION.md`](../../docs/VALIDATION.md#8-the-built-in-piv).
