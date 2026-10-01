"""State of the simulation: particle and node arrays.

It replaces the ``COMMON /PARTICULAS/`` and ``COMMON /NODOS/`` blocks of the original.
``NamedTuple`` is used so that Numba can take them straight into the kernels.
"""

from __future__ import annotations

from typing import NamedTuple

import numpy as np


class Particles(NamedTuple):
    """Per-particle variables. ``n`` = number of particles (NP)."""

    position: np.ndarray  # XP (n, 2)
    initial_position: np.ndarray  # position when the analysis started (the GiD mesh)
    displacement: np.ndarray  # UP (n, 2): accumulated displacement
    step_displacement: np.ndarray  # UPO (n, 2): displacement of this step
    position_increment: np.ndarray  # XPP (n, 2)
    velocity: np.ndarray  # VP (n, 2)
    acceleration: np.ndarray  # ACEP (n, 2)
    strain: np.ndarray  # EPS (n, 4): xx, yy, xy (engineering), zz (always 0)
    strain_increment: np.ndarray  # DEPS (n, 3)
    eq_strain: np.ndarray  # EPSEQ: accumulated equivalent shear strain
    eq_strain_increment: np.ndarray  # EPSEQ2
    vol_strain: np.ndarray  # EVOLUMETRIC
    vol_strain_increment: np.ndarray  # EVOL_IN
    #: Curl of the velocity field at the particle, dv/dx - du/dy, in 1/s. Positive
    #: counter-clockwise. The antisymmetric half of the same gradient the strains come from,
    #: so it costs almost nothing and it is what tells rotation apart from shear.
    vorticity: np.ndarray
    #: Rotation accumulated by the particle, in degrees, from adding vorticity*dt/2.
    rotation: np.ndarray
    mass: np.ndarray  # AMP (always 1)
    potential_energy: np.ndarray  # E_Potential
    kinetic_energy: np.ndarray  # E_Cinetic_x, E_Cinetic_y (n, 2)
    total_energy: np.ndarray  # E_Total
    moisture: np.ndarray  # SMOIST_NP
    saturation: np.ndarray  # SATURA_NP
    nan_initial: np.ndarray  # NaN_P: 0 active, 1 without data on step 1
    nan_step: np.ndarray  # NaN_P2: without data on the current step (informative only)
    missing_data: np.ndarray  # the "NaNs" result: missing data around the particle
    lost: np.ndarray  # IDONDE == -1 (size max(n, PIVlab nodes))

    @classmethod
    def zeros(cls, n: int, n_lost: int) -> Particles:
        def vec(k: int = 0) -> np.ndarray:
            return np.zeros((n, k) if k else n, dtype=np.float64)

        return cls(
            position=vec(2), initial_position=vec(2), displacement=vec(2),
            step_displacement=vec(2),
            position_increment=vec(2), velocity=vec(2), acceleration=vec(2),
            strain=vec(4), strain_increment=vec(3), eq_strain=vec(), eq_strain_increment=vec(),
            vol_strain=vec(), vol_strain_increment=vec(),
            vorticity=vec(), rotation=vec(), mass=np.ones(n),
            potential_energy=vec(), kinetic_energy=vec(2), total_energy=vec(),
            moisture=vec(), saturation=vec(),
            nan_initial=np.zeros(n, dtype=np.int8), nan_step=np.zeros(n, dtype=np.int8),
            missing_data=np.zeros(n, dtype=np.int64),
            lost=np.zeros(max(n, n_lost), dtype=np.bool_),
        )


class Nodes(NamedTuple):
    """Nodal variables.

    The measured fields (``velocity``, ``is_nan``…) are indexed with the PIV-NP node
    numbering (= the PIVlab points). The computation fields (``momentum``, ``mass``…) are
    indexed with the nodes of the particle grid, which is the same one when IVERSION = 1 and
    the staggered grid when IVERSION = 2.
    """

    velocity: np.ndarray  # VEL_X_NODO, VEL_Y_NODO (nn, 2)
    previous_velocity: np.ndarray  # VEL_X_NODO_V, VEL_Y_NODO_V (nn, 2)
    is_nan: np.ndarray  # Node_NaN: PIVlab measured no velocity at the node
    filled: np.ndarray  # velocity rebuilt by the boundary correction
    moisture_measured: np.ndarray  # SMOISTURE_N2
    saturation_measured: np.ndarray  # SATURATION_N2
    momentum: np.ndarray  # PV (nm, 2)
    momentum_increment: np.ndarray  # APV (nm, 2)
    mass: np.ndarray  # AM (nm)
    moisture: np.ndarray  # SMOISTURE_N (nm)
    saturation: np.ndarray  # SATURATION_N (nm)
    active_count: np.ndarray  # ICOUNT_NO_NAN (nm)

    @classmethod
    def zeros(cls, n_measured: int, n_mesh: int) -> Nodes:
        return cls(
            velocity=np.zeros((n_measured, 2)),
            previous_velocity=np.zeros((n_measured, 2)),
            is_nan=np.zeros(n_measured, dtype=np.int8),
            filled=np.zeros(n_measured, dtype=np.bool_),
            moisture_measured=np.zeros(n_measured),
            saturation_measured=np.zeros(n_measured),
            momentum=np.zeros((n_mesh, 2)),
            momentum_increment=np.zeros((n_mesh, 2)),
            mass=np.zeros(n_mesh),
            moisture=np.zeros(n_mesh),
            saturation=np.zeros(n_mesh),
            active_count=np.zeros(n_mesh, dtype=np.int64),
        )
