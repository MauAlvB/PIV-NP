"""Phase 2: the strain from the deformation gradient, differentiating once instead of N times.

PIV-NP builds the strain by adding up a linear increment at every step. That has two known
costs. It invents strain for a rigid rotation -- the synthetic rotation case turns 50 degrees
without deforming and still reports an equivalent strain of 0.005, as an apparent 1.5 %
contraction. And it adds one step's worth of differentiation noise twenty times, which phase
3 showed is the same size as the signal itself.

The alternative is to differentiate once, over the whole analysis. The particles start on a
regular lattice, so the accumulated displacement can be differentiated with respect to the
initial coordinates directly, giving the deformation gradient

    F = I + du/dX

from which everything follows exactly rather than incrementally:

* Green-Lagrange strain  E = (F'F - I)/2, zero for any rigid motion by construction
* rotation, from the polar decomposition, in closed form for 2x2 so no LAPACK is needed
* volume change, det(F) - 1, which is exact instead of the trace of a linear increment

Nothing here is wired into PIV-NP; it is computed from what a finished analysis leaves in
memory, so the two can be compared on the same run.
"""
from __future__ import annotations

import numpy as np

from pivnp.state import Particles


def lattice(initial: np.ndarray) -> tuple[np.ndarray, np.ndarray, float, float]:
    """Row and column of every particle on the lattice it started from.

    Worked out from the initial positions themselves rather than from the order particles
    are created in, so it does not depend on that convention.
    """
    xs = np.unique(np.round(initial[:, 0], 9))
    ys = np.unique(np.round(initial[:, 1], 9))
    col = np.searchsorted(xs, np.round(initial[:, 0], 9))
    row = np.searchsorted(ys, np.round(initial[:, 1], 9))
    dx = float(np.median(np.diff(xs))) if xs.size > 1 else 1.0
    dy = float(np.median(np.diff(ys))) if ys.size > 1 else 1.0
    return row, col, dx, dy


def _gradient(field: np.ndarray, spacing: float, axis: int) -> np.ndarray:
    """Central difference where both neighbours exist, one-sided where only one does."""
    out = np.full(field.shape, np.nan)
    forward = np.roll(field, -1, axis=axis)
    backward = np.roll(field, 1, axis=axis)
    edge = [slice(None)] * field.ndim
    edge[axis] = 0
    backward[tuple(edge)] = np.nan
    edge[axis] = -1
    forward[tuple(edge)] = np.nan

    both = np.isfinite(forward) & np.isfinite(backward)
    out = np.where(both, (forward - backward) / (2.0 * spacing), out)
    only_forward = ~both & np.isfinite(forward) & np.isfinite(field)
    out = np.where(only_forward, (forward - field) / spacing, out)
    only_backward = ~both & np.isfinite(backward) & np.isfinite(field)
    out = np.where(only_backward, (field - backward) / spacing, out)
    return out


def deformation_gradient(particles: Particles, usable: np.ndarray):
    """``F`` of every particle, from the accumulated displacement on the initial lattice.

    ``usable`` marks the particles whose displacement means something: the ones still inside
    the grid and measured from the start. The rest are left out of the differences, so a
    particle that stopped moving when it left the mesh does not drag its neighbours.
    """
    initial = particles.initial_position
    row, col, dx, dy = lattice(initial)
    shape = (row.max() + 1, col.max() + 1)

    fields = []
    for component in (0, 1):
        grid = np.full(shape, np.nan)
        grid[row[usable], col[usable]] = particles.displacement[usable, component]
        fields.append(grid)

    # du/dX and du/dY for each component; rows run along Y, columns along X
    dudx = _gradient(fields[0], dx, axis=1)
    dudy = _gradient(fields[0], dy, axis=0)
    dvdx = _gradient(fields[1], dx, axis=1)
    dvdy = _gradient(fields[1], dy, axis=0)

    f00 = 1.0 + dudx[row, col]
    f01 = dudy[row, col]
    f10 = dvdx[row, col]
    f11 = 1.0 + dvdy[row, col]
    return f00, f01, f10, f11


def green_lagrange(f00, f01, f10, f11):
    """``E = (F'F - I)/2``: exactly zero for any rigid motion, rotation included."""
    exx = 0.5 * (f00 * f00 + f10 * f10 - 1.0)
    eyy = 0.5 * (f01 * f01 + f11 * f11 - 1.0)
    exy = 0.5 * (f00 * f01 + f10 * f11)
    return exx, eyy, 2.0 * exy          # the third one is the engineering shear


def rotation_degrees(f00, f01, f10, f11):
    """Rotation of the polar decomposition, in closed form for 2x2.

    ``F = R U``; the angle of ``R`` is ``atan2(F10 - F01, F00 + F11)``. For a pure rotation
    it gives the angle exactly, whatever its size, with no small-angle assumption.
    """
    return np.degrees(np.arctan2(f10 - f01, f00 + f11))


def volume_change(f00, f01, f10, f11):
    """``det(F) - 1``: the true change of area, not the trace of a linear increment."""
    return f00 * f11 - f01 * f10 - 1.0


def equivalent_shear(exx, eyy, gxy):
    """The same measure PIV-NP publishes, so the two can be put side by side."""
    mean = (exx + eyy) / 3.0
    dx, dy, dz = exx - mean, eyy - mean, -mean
    j2 = (dx * dx + dy * dy + dz * dz) / 2.0 + (gxy / 2.0) ** 2
    return 2.0 * np.sqrt(3.0 * np.maximum(j2, 0.0)) / 3.0
