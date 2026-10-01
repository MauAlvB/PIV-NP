"""Strain from the deformation gradient: differentiating once instead of once per step.

PIV-NP builds its strain by adding a linear increment at every step, which is what the
original Fortran did. That has two costs, both measured in ``docs/VALIDATION.md``:

* it reports strain for a **rigid rotation**, which does not deform anything. On the
  ``rotation`` test cases, 50° of turning comes out as an equivalent shear of 0.005 and an
  apparent 1.5 % loss of volume that never happened;
* it adds one step's worth of differentiation noise at every step, and on real data that
  noise is the same size as the signal.

The particles start on a regular lattice, so the accumulated displacement can be
differentiated against the **initial** coordinates directly. That gives the deformation
gradient

    F = I + du/dX

and from it the Green-Lagrange strain ``E = (F'F - I)/2``, which is identically zero for any
rigid motion; the rotation of the polar decomposition, in closed form for 2x2; and the true
change of area ``det(F) - 1``.

These results are published **beside** the incremental ones rather than instead of them, so
that analyses made with either can be compared before anything is retired. Everything here is
computed from what the analysis already keeps -- ``initial_position`` and ``displacement`` --
so there is no extra state and nothing to carry across a restart.
"""

from __future__ import annotations

import numpy as np

from .state import Particles


def initial_lattice(initial_position: np.ndarray) -> tuple[np.ndarray, np.ndarray,
                                                           float, float]:
    """Row, column and spacing of the regular lattice the particles were created on.

    Worked out from the positions themselves rather than from the order particles are
    generated in, so it does not depend on that convention.
    """
    xs = np.unique(np.round(initial_position[:, 0], 9))
    ys = np.unique(np.round(initial_position[:, 1], 9))
    column = np.searchsorted(xs, np.round(initial_position[:, 0], 9))
    row = np.searchsorted(ys, np.round(initial_position[:, 1], 9))
    dx = float(np.median(np.diff(xs))) if xs.size > 1 else 1.0
    dy = float(np.median(np.diff(ys))) if ys.size > 1 else 1.0
    return row, column, dx, dy


def _derivative(field: np.ndarray, spacing: float, axis: int) -> np.ndarray:
    """Central difference where both neighbours are known, one-sided where only one is."""
    forward = np.roll(field, -1, axis=axis)
    backward = np.roll(field, 1, axis=axis)
    first = [slice(None)] * field.ndim
    first[axis] = 0
    backward[tuple(first)] = np.nan
    last = [slice(None)] * field.ndim
    last[axis] = -1
    forward[tuple(last)] = np.nan

    both = np.isfinite(forward) & np.isfinite(backward)
    out = np.where(both, (forward - backward) / (2.0 * spacing), np.nan)
    ahead = ~both & np.isfinite(forward) & np.isfinite(field)
    out = np.where(ahead, (forward - field) / spacing, out)
    behind = ~both & np.isfinite(backward) & np.isfinite(field)
    return np.where(behind, (field - backward) / spacing, out)


def deformation_gradient(particles: Particles,
                         usable: np.ndarray) -> tuple[np.ndarray, ...]:
    """The four components of ``F`` for every particle, as ``(F00, F01, F10, F11)``.

    ``usable`` marks the particles whose accumulated displacement means something: still
    inside the grid, and measured from the first step. The others are left out of the
    differences, so a particle frozen when it left the mesh does not drag its neighbours.
    Particles without enough neighbours come back as ``NaN``.
    """
    row, column, dx, dy = initial_lattice(particles.initial_position)
    shape = (int(row.max()) + 1, int(column.max()) + 1)

    gradients = []
    for component in (0, 1):
        laid_out = np.full(shape, np.nan)
        laid_out[row[usable], column[usable]] = particles.displacement[usable, component]
        gradients.append((_derivative(laid_out, dx, axis=1),    # d/dX
                          _derivative(laid_out, dy, axis=0)))   # d/dY

    (dudx, dudy), (dvdx, dvdy) = gradients
    return (1.0 + dudx[row, column], dudy[row, column],
            dvdx[row, column], 1.0 + dvdy[row, column])


def green_lagrange(f: tuple[np.ndarray, ...]) -> tuple[np.ndarray, ...]:
    """``E = (F'F - I)/2`` as ``(xx, yy, engineering xy)``, zero for any rigid motion."""
    f00, f01, f10, f11 = f
    xx = 0.5 * (f00 * f00 + f10 * f10 - 1.0)
    yy = 0.5 * (f01 * f01 + f11 * f11 - 1.0)
    xy = 0.5 * (f00 * f01 + f10 * f11)
    return xx, yy, 2.0 * xy


def rotation_degrees(f: tuple[np.ndarray, ...]) -> np.ndarray:
    """Rotation of the polar decomposition ``F = R U``, in degrees.

    For 2x2 the angle of ``R`` is ``atan2(F10 - F01, F00 + F11)``, which is exact for a
    rotation of any size -- unlike accumulating the curl, which is a small-angle sum.
    """
    f00, f01, f10, f11 = f
    return np.degrees(np.arctan2(f10 - f01, f00 + f11))


def area_change(f: tuple[np.ndarray, ...]) -> np.ndarray:
    """``det(F) - 1``: the true change of area, not the trace of a linear increment."""
    f00, f01, f10, f11 = f
    return f00 * f11 - f01 * f10 - 1.0


def equivalent_shear(strain: tuple[np.ndarray, ...], j2_threshold: float = 0.0) -> np.ndarray:
    """The same equivalent measure the incremental strain is published with.

    Written out rather than reusing ``solver.deviatoric_q`` because that one is compiled for
    scalars inside the particle loop; this works on the whole array at once.
    """
    xx, yy, xy = strain
    mean = (xx + yy) / 3.0
    dx, dy, dz = xx - mean, yy - mean, -mean
    j2 = (dx * dx + dy * dy + dz * dz) / 2.0 + (xy / 2.0) ** 2
    return 2.0 * np.sqrt(3.0 * np.where(j2 > j2_threshold, j2, 0.0)) / 3.0
