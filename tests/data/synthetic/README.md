# Synthetic validation cases

Velocity fields built by hand on small meshes, with a known solution. They check the
computation itself, not just that it has not changed with respect to earlier versions.
`cases.json` records the original (Spanish) name each case had in the group's collection,
so the provenance of the data can still be traced.

Every folder holds the PIVlab files (`datos (n).txt`), the `.PAR` **in the old format**
(block 3 with 3 values and no block 4, which also checks that it is still read) and, under
`expected/`, the `.POST.RES` of the analysis the group ran at the time.

| Case | Mesh | Imposed field | Expected result |
|---|---|---|---|
| `displacement_1P`, `_4P` | 4×4 and 5×5 cells of 0.5 m | 10 s at 2 cm/s along +x and 2 cm/s along −y on every node | u = (+0.20, −0.20) m on every particle, no strain |
| `shear_1P`, `_4P` | 5×4 cells | 10 s with a horizontal velocity of 2, 1.5, 1, 0.5 and 0 cm/s from top to bottom | γxy = 0.10 uniform, εxx = εyy = 0, εq = 0.1/√3 |
| `horizontal_uniform_1P`, `_4P` | 5×4 cells | 24 s of uniform stretching | εxx = 0.24 uniform, εyy = γxy = 0 |
| `horizontal_1P`, `_4P` | 5×4 cells | 15 s with a velocity of 0, 0.1, 0.3, 0.9 and 2.4 cm/s per column | a different εxx in each column of cells: 0.03, 0.06, 0.18… |
| `rotation_1P`, `_4P` | 6×6 cells | 50 s rotating at 1°/s around the central node | the solid does not deform, so the strain should be 0 |

In several cases some particles leave the mesh and stop being computed; that is why the
number of active particles at the end is smaller than at the start.

## About the rotation

A rigid-body rotation does not deform the material, so every strain that comes out is
error. Measuring the case shows that it has **two distinct origins**:

| | Equivalent shear strain (mean) |
|---|---|
| No contour correction | 0.0169 |
| Average of the neighbouring nodes | 0.0111 |
| Average of the surrounding particles | 0.0125 |
| Extrapolation from the interior | **0.0051, the same on every particle** |

The case has instants with 24 of its 49 nodes without data, so **two thirds of the error
came from the boundary**. Once those nodes are rebuilt by extrapolation, the remaining
error is uniform across the solid and matches what the theory predicts: accumulating
linear strain increments through a finite rotation leaves εxx = εyy = n·(cos Δθ − 1),
which for 50 steps of 1° gives an equivalent shear strain of 0.00508, against the 0.0051
measured.

That residual really is a limitation of the formulation, not of this implementation:
removing it would need a finite-strain measure (from the deformation gradient). The tests
pin both values down so that any change making them worse is caught.
