"""Update of the particles at every step (formerly the ``SOLMOV`` and ``INVAR2`` subroutines).

Every particle is computed independently, so the loops run in parallel with
``numba.prange``. The expressions keep the operation order of the original in order to
reproduce its rounding.
"""

from __future__ import annotations

import math

import numpy as np
from numba import njit, prange

from .config import CaseConfig
from .constants import (
    GRAVITY,
    J2_THRESHOLD,
    LEGACY_GRAVITY_STEP,
    LEGACY_J2_THRESHOLD,
    MACHINE_EPSILON,
    NODE_SIGN_X,
    NODE_SIGN_Y,
)
from .mesh import Grid, cell_node
from .nodal import grid_locate
from .particles import update_lost_flags
from .state import Nodes, Particles

#: Values of ``Particles.nan_initial``.
ACTIVE, NAN_AT_START = 0, 1


@njit(cache=True, nogil=True, inline="always")
def deviatoric_q(sx, sy, sz, sxy, threshold=J2_THRESHOLD):
    """Deviatoric stress (or strain) q = sqrt(3·J2); 0 if J2 <= ``threshold`` (``INVAR2``)."""
    mean = (sx + sy + sz) / 3.0
    dx, dy, dz = sx - mean, sy - mean, sz - mean
    j2 = (dx * dx + dy * dy + dz * dz) / 2.0 + sxy * sxy
    return math.sqrt(3.0 * j2) if j2 > threshold else 0.0


@njit(parallel=True, cache=True)
def _advance_kernel(lost, cells, cell_x, cell_y, n_cols, dx, dy, dt,
                    momentum, momentum_increment, nodal_mass, node_is_nan,
                    nodal_moisture, nodal_saturation,
                    position, velocity, acceleration, increment, displacement, step_disp,
                    moisture, saturation, nan_initial, nan_step,
                    step, restart, with_moisture, mesh_version, legacy_restart_displacement):
    for i in prange(position.shape[0]):
        if not lost[i]:
            cell = cells[i]
            xx = 2 * (position[i, 0] - (cell_x[i] + dx / 2.0)) / dx
            yy = 2 * (position[i, 1] - (cell_y[i] + dy / 2.0)) / dy
            for k in range(2):
                increment[i, k] = 0.0
                velocity[i, k] = 0.0
                acceleration[i, k] = 0.0
                step_disp[i, k] = 0.0
            nan_step[i] = 0
            if step == 1 and not restart:
                nan_initial[i] = 0
            if step == 1 and legacy_restart_displacement:
                # The original added the accumulated displacement to the "instantaneous" one
                # on the first step after a restart.
                step_disp[i, 0] = displacement[i, 0]
                step_disp[i, 1] = displacement[i, 1]
            nan_corners = 0
            moisture[i] = 0.0
            saturation[i] = 0.0

            for j in range(4):
                fn = (1.0 + xx * NODE_SIGN_X[j]) * (1.0 + yy * NODE_SIGN_Y[j]) / 4.0
                node = cell_node(cell, n_cols, j)
                if nodal_mass[node] >= MACHINE_EPSILON:
                    f1 = fn
                    f2 = fn * dt
                else:
                    f1 = 0.0
                    f2 = 0.0
                for k in range(2):
                    velocity[i, k] = velocity[i, k] + momentum[node, k] * f1
                    acceleration[i, k] = acceleration[i, k] + momentum_increment[node, k] * f1 / dt
                    increment[i, k] = increment[i, k] + momentum[node, k] * f2
                    step_disp[i, k] = step_disp[i, k] + momentum[node, k] * f2
                    displacement[i, k] = displacement[i, k] + momentum[node, k] * f2
                if with_moisture:
                    moisture[i] = moisture[i] + nodal_moisture[node] * f1
                    saturation[i] = saturation[i] + nodal_saturation[node] * f1
                if mesh_version == 1:
                    nan_corners += node_is_nan[node]
                elif node_is_nan[cell] == 1:  # staggered cell centred on a node without data
                    if step == 1 and not restart:
                        nan_initial[i] = NAN_AT_START
                    nan_step[i] = 1

            if mesh_version == 1 and nan_corners == 4:  # all 4 nodes without data
                if step == 1 and not restart:
                    nan_initial[i] = NAN_AT_START
                else:
                    nan_step[i] = 1

        if nan_initial[i] == NAN_AT_START:  # no data from the start: it does not move
            increment[i, 0] = 0.0
            increment[i, 1] = 0.0
            displacement[i, 0] = 0.0
            displacement[i, 1] = 0.0
            lost[i] = True


def advance_particles(particles: Particles, nodes: Nodes, grid: Grid, config: CaseConfig,
                      step: int, legacy_compat: bool = False) -> None:
    """Interpolate particle velocity, acceleration and displacement from the nodes.

    First half of ``SOLMOV``. Marks as ``NAN_AT_START`` the particles that on the first step
    sit in cells without any data (outside the material).
    """
    p = particles
    n = p.position.shape[0]
    cells, cell_x, cell_y = grid_locate(p.position, grid)
    update_lost_flags(p.lost[:n], cells, step)
    _advance_kernel(
        p.lost, cells, cell_x, cell_y, grid.n_cols, grid.dx, grid.dy, config.dt,
        nodes.momentum, nodes.momentum_increment, nodes.mass, nodes.is_nan,
        nodes.moisture, nodes.saturation,
        p.position, p.velocity, p.acceleration, p.position_increment, p.displacement,
        p.step_displacement, p.moisture, p.saturation, p.nan_initial, p.nan_step,
        step, config.restart, config.moisture, config.mesh_version,
        config.restart and legacy_compat,
    )


@njit(parallel=True, cache=True)
def _strain_kernel(lost, cells, n_cols, dx, dy, dt, momentum, nodal_mass,
                   position, increment, velocity, mass, strain, strain_inc,
                   vol_strain, vol_strain_inc, potential, kinetic, total, moisture,
                   eq_strain, eq_strain_inc, vorticity, rotation,
                   legacy_divide_by_mass, gravity, j2_threshold):
    for i in prange(position.shape[0]):
        if lost[i]:
            continue
        cell = cells[i]
        d0 = 0.0
        d1 = 0.0
        d2 = 0.0
        d3 = 0.0
        for j in range(4):
            node = cell_node(cell, n_cols, j)
            # Derivatives of the shape functions at the centre of the element: the strain is
            # uniform within each cell.
            dndx = NODE_SIGN_X[j] * 0.5 / dx
            dndy = NODE_SIGN_Y[j] * 0.5 / dy
            if nodal_mass[node] >= MACHINE_EPSILON:
                # The original divided by the nodal mass, which with IVERSION=2 is about
                # NPC² and leaves the strains scaled as the velocity divided by NPC². Only
                # reproduced in compatibility mode.
                if legacy_divide_by_mass:
                    fx = dndx * dt / nodal_mass[node]
                    fy = dndy * dt / nodal_mass[node]
                else:
                    fx = dndx * dt
                    fy = dndy * dt
            else:
                fx = 0.0
                fy = 0.0
            d0 = d0 + momentum[node, 0] * fx
            d1 = d1 + momentum[node, 1] * fy
            d2 = d2 + momentum[node, 0] * fy + momentum[node, 1] * fx
            # The antisymmetric counterpart of d2: the same two derivatives with a minus,
            # which is the curl. Nothing else to look up, so it comes practically free.
            d3 = d3 + momentum[node, 1] * fx - momentum[node, 0] * fy
        if abs(d0) < MACHINE_EPSILON and d0 != 0.0:
            d0 = 0.0
        if abs(d1) < MACHINE_EPSILON and d1 != 0.0:
            d1 = 0.0
        if abs(d2) < MACHINE_EPSILON and d2 != 0.0:
            d2 = 0.0
        if abs(d3) < MACHINE_EPSILON and d3 != 0.0:
            d3 = 0.0
        # d3 is the curl times dt. Half of it is the rotation of the material element over
        # the step, which accumulates; the curl itself is reported as it stands.
        vorticity[i] = d3 / dt
        rotation[i] = rotation[i] + 0.5 * d3 * 180.0 / np.pi
        strain_inc[i, 0] = d0
        strain_inc[i, 1] = d1
        strain_inc[i, 2] = d2
        strain[i, 0] = strain[i, 0] + d0
        strain[i, 1] = strain[i, 1] + d1
        strain[i, 2] = strain[i, 2] + d2

        vol_strain[i] = strain[i, 0] + strain[i, 1]
        vol_strain_inc[i] = d0 + d1
        potential[i] = mass[i] * gravity * position[i, 1]
        kinetic[i, 0] = 0.5 * mass[i] * velocity[i, 0] * velocity[i, 0]
        kinetic[i, 1] = 0.5 * mass[i] * velocity[i, 1] * velocity[i, 1]
        total[i] = potential[i] + kinetic[i, 0] + kinetic[i, 1]
        if moisture[i] < 0.0:
            moisture[i] = 0.0

        position[i, 0] = position[i, 0] + increment[i, 0]
        position[i, 1] = position[i, 1] + increment[i, 1]

        # Equivalent shear strain: 2/3·q (with the engineering shear strain / 2)
        eq_strain[i] = 2.0 * deviatoric_q(strain[i, 0], strain[i, 1], strain[i, 3],
                                          strain[i, 2] / 2.0, j2_threshold) / 3.0
        eq_strain_inc[i] = 2.0 * deviatoric_q(d0, d1, 0.0, d2 / 2.0, j2_threshold) / 3.0


def update_strains(particles: Particles, nodes: Nodes, grid: Grid, config: CaseConfig,
                   step: int, legacy_compat: bool = False) -> None:
    """Strains, energies and new position of every particle (second half of SOLMOV)."""
    p = particles
    n = p.position.shape[0]
    cells, _, _ = grid_locate(p.position, grid)
    update_lost_flags(p.lost[:n], cells, step)
    _strain_kernel(
        p.lost, cells, grid.n_cols, grid.dx, grid.dy, config.dt, nodes.momentum, nodes.mass,
        p.position, p.position_increment, p.velocity, p.mass, p.strain, p.strain_increment,
        p.vol_strain, p.vol_strain_increment, p.potential_energy, p.kinetic_energy,
        p.total_energy, p.moisture, p.eq_strain, p.eq_strain_increment,
        p.vorticity, p.rotation, legacy_compat,
        LEGACY_GRAVITY_STEP if legacy_compat else GRAVITY,
        LEGACY_J2_THRESHOLD if legacy_compat else J2_THRESHOLD,
    )


#: Below this, a strain or a rotation is numerical noise rather than deformation.
_NOTHING_HAPPENED = 1e-12


def _turn_and_shear(particles: Particles) -> tuple[np.ndarray, np.ndarray]:
    """The two magnitudes the kinematics are judged by, from the accumulated deformation.

    ``turn`` is twice the rotation in radians, which is what the curl integrates to, and
    ``shear`` is the in-plane deviatoric strain. For simple shear the two are equal, which
    is what makes their ratio 1 there.
    """
    strain = particles.strain
    turn = 2.0 * np.abs(np.radians(particles.rotation))
    shear = np.hypot(strain[:, 0] - strain[:, 1], strain[:, 2])
    return turn, shear


def vorticity_number(particles: Particles) -> np.ndarray:
    """Kinematic vorticity number of the deformation accumulated so far.

    ``0`` pure shear (it deforms without turning), ``1`` simple shear (a shear band), and
    the larger it grows the more the rotation dominates. ``NaN`` where it is not defined:
    a rigid rotation has no deviatoric strain to divide by, and a particle that neither
    turned nor deformed has no kinematics to describe.

    Taking the ratio of the accumulated quantities is how the number is estimated from
    finite strain in deformed rocks. It equals the time average of the instantaneous
    vorticity number when the deformation is steady; when it is not, it describes the
    deformation as a whole. The instantaneous one can be had from ``Vorticity`` and
    ``Inc_strain``, which are both published.
    """
    turn, shear = _turn_and_shear(particles)
    out = np.full(turn.shape, np.nan)
    usable = shear > _NOTHING_HAPPENED
    out[usable] = turn[usable] / shear[usable]
    return out


def rotation_angle(particles: Particles) -> np.ndarray:
    """The same thing bounded, in degrees, so that it can be drawn without a singularity.

    ``0`` pure shear, ``45`` simple shear, ``90`` rigid rotation. It is the arctangent of
    :func:`vorticity_number`, so one converts into the other: ``Wm = tan(angle)``. ``NaN``
    where nothing happened at all.
    """
    turn, shear = _turn_and_shear(particles)
    out = np.degrees(np.arctan2(turn, shear))
    out[(turn <= _NOTHING_HAPPENED) & (shear <= _NOTHING_HAPPENED)] = np.nan
    return out


def output_mask(particles: Particles, grid: Grid, step: int) -> np.ndarray:
    """Particles located in the grid when printing (``IDONDE /= -1`` in IMPRES_GiD)."""
    n = particles.position.shape[0]
    cells, _, _ = grid_locate(particles.position, grid)
    update_lost_flags(particles.lost[:n], cells, step)
    return ~particles.lost[:n]


@njit(parallel=True, cache=True)
def _count_nan_kernel(lost, cells, n_cols, node_is_nan, mesh_version, out):
    for i in prange(lost.size):
        if lost[i]:
            out[i] = 0
            continue
        if mesh_version == 1:  # nodes of the element without data (0 to 4)
            total = 0
            for j in range(4):
                total += node_is_nan[cell_node(cells[i], n_cols, j)]
            out[i] = total
        else:  # the PIVlab point at the centre of the element has no data (0 or 1)
            out[i] = node_is_nan[cells[i]]


def count_nan_nodes(particles: Particles, nodes: Nodes, grid: Grid, mesh_version: int,
                    step: int) -> np.ndarray:
    """How much data is missing around each particle (the ``NaNs`` result)."""
    n = particles.position.shape[0]
    cells, _, _ = grid_locate(particles.position, grid)
    update_lost_flags(particles.lost[:n], cells, step)
    _count_nan_kernel(particles.lost[:n], cells, grid.n_cols, nodes.is_nan, mesh_version,
                      particles.missing_data)
    return particles.missing_data
