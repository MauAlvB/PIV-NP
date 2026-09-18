"""Estado de la simulación: arrays de partículas y de nodos.

Sustituye a los bloques ``COMMON /PARTICULAS/`` y ``COMMON /NODOS/`` del original.
Se usan ``NamedTuple`` para que Numba pueda recibirlos directamente en los kernels.
"""

from __future__ import annotations

from typing import NamedTuple

import numpy as np


class Particles(NamedTuple):
    """Variables por partícula numérica. ``n`` = número de partículas (NP)."""

    position: np.ndarray  # XP (n, 2)
    displacement: np.ndarray  # UP (n, 2): desplazamiento acumulado
    step_displacement: np.ndarray  # UPO (n, 2): desplazamiento del paso
    position_increment: np.ndarray  # XPP (n, 2)
    velocity: np.ndarray  # VP (n, 2)
    acceleration: np.ndarray  # ACEP (n, 2)
    strain: np.ndarray  # EPS (n, 4): xx, yy, xy (ingenieril), zz (siempre 0)
    strain_increment: np.ndarray  # DEPS (n, 3)
    eq_strain: np.ndarray  # EPSEQ: deformación de corte equivalente acumulada
    eq_strain_increment: np.ndarray  # EPSEQ2
    vol_strain: np.ndarray  # EVOLUMETRIC
    vol_strain_increment: np.ndarray  # EVOL_IN
    mass: np.ndarray  # AMP (siempre 1)
    potential_energy: np.ndarray  # E_Potential
    kinetic_energy: np.ndarray  # E_Cinetic_x, E_Cinetic_y (n, 2)
    total_energy: np.ndarray  # E_Total
    moisture: np.ndarray  # SMOIST_NP
    saturation: np.ndarray  # SATURA_NP
    nan_initial: np.ndarray  # NaN_P: 0 activa, 1 sin datos en el paso 1, 2 seguimiento PTV
    nan_step: np.ndarray  # NaN_P2: sin datos en el paso actual (solo informativo)
    lost: np.ndarray  # IDONDE == -1 (tamaño max(n, nodos PIVlab), ver H-10)

    @classmethod
    def zeros(cls, n: int, n_lost: int) -> Particles:
        def vec(k: int = 0) -> np.ndarray:
            return np.zeros((n, k) if k else n, dtype=np.float64)

        return cls(
            position=vec(2), displacement=vec(2), step_displacement=vec(2),
            position_increment=vec(2), velocity=vec(2), acceleration=vec(2),
            strain=vec(4), strain_increment=vec(3), eq_strain=vec(), eq_strain_increment=vec(),
            vol_strain=vec(), vol_strain_increment=vec(), mass=np.ones(n),
            potential_energy=vec(), kinetic_energy=vec(2), total_energy=vec(),
            moisture=vec(), saturation=vec(),
            nan_initial=np.zeros(n, dtype=np.int8), nan_step=np.zeros(n, dtype=np.int8),
            lost=np.zeros(max(n, n_lost), dtype=np.bool_),
        )


class Nodes(NamedTuple):
    """Variables nodales.

    Los campos medidos (``velocity``, ``is_nan``...) se indexan con la numeración de los
    nodos PIV-NP (= puntos PIVlab). Los campos de cálculo (``momentum``, ``mass``...) se
    indexan con los nodos de la malla de partículas, que coincide con la anterior si
    IVERSION = 1 y es la malla desplazada si IVERSION = 2.
    """

    velocity: np.ndarray  # VEL_X_NODO, VEL_Y_NODO (nn, 2)
    previous_velocity: np.ndarray  # VEL_X_NODO_V, VEL_Y_NODO_V (nn, 2)
    is_nan: np.ndarray  # Node_NaN: PIVlab no midió velocidad en el nodo
    filled: np.ndarray  # velocidad reconstruida por la corrección de contorno
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
