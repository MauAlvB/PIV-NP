"""Numerical constants.

Several literals of the original Fortran code were written without the ``d0`` suffix (say
``9.81`` or ``1.E-10``), so Fortran stored them in single precision (REAL*4) before promoting
them to REAL*8, losing precision from the seventh digit on. Here they are used in double
precision; the rounded versions are kept with the ``LEGACY_`` prefix for the compatibility
mode.
"""

from __future__ import annotations

import numpy as np


def as_fortran_real4(value: float) -> float:
    """Return ``value`` rounded to single precision, like a Fortran REAL literal."""
    return float(np.float32(value))


#: ``EPSILON(1.d0)``: threshold below which a nodal mass or a strain increment counts as zero.
MACHINE_EPSILON: float = float(np.finfo(np.float64).eps)

#: Gravity [m/s²] for the potential energy.
GRAVITY: float = 9.81
LEGACY_GRAVITY_STEP: float = as_fortran_real4(9.81)  # ``9.81`` (REAL*4) in SOLMOV

#: Threshold of the second invariant J2 below which q = 0 (``1.E-10`` in INVAR2).
J2_THRESHOLD: float = 1.0e-10
LEGACY_J2_THRESHOLD: float = as_fortran_real4(1.0e-10)

#: Local coordinates (in [-1, 1]) of the NPC rows/columns of particles of each cell.
#: NPC = 2 and 3 are uniform subdivisions (centres of the subcells); NPC = 4, 5 and 6 are
#: Gauss-Legendre points.
PARTICLE_LOCAL_COORDS: dict[int, tuple[float, ...]] = {
    1: (0.0,),
    2: (-0.5, 0.5),
    3: (-2 / 3, 0.0, 2 / 3),
    4: (-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526),
    5: (-0.9061798459386640, -0.5384693101056831, 0.0, 0.5384693101056831,
        0.9061798459386640),
    6: (-0.9324695142031521, -0.6612093864662645, -0.2386191860831969,
        0.2386191860831969, 0.6612093864662645, 0.9324695142031521),
}

#: The same coordinates as the original stored them (REAL*4 literals).
LEGACY_PARTICLE_LOCAL_COORDS: dict[int, tuple[float, ...]] = {
    1: (0.0,),
    2: (-0.5, 0.5),
    3: (as_fortran_real4(-0.66666666666667), 0.0, as_fortran_real4(0.66666666666667)),
    4: tuple(
        as_fortran_real4(v)
        for v in (-0.861136311594053, -0.339981043584856, 0.339981043584856, 0.861136311594053)
    ),
    5: tuple(
        as_fortran_real4(v)
        for v in (-0.906179845988664, -0.538469310105683, 0.0, 0.538469310105683,
                  0.906179845988664)
    ),
    6: tuple(
        as_fortran_real4(v)
        for v in (-0.932469514203152, -0.661209386466265, -0.238619186083197,
                  0.238619186083197, 0.661209386466265, 0.932469514203152)
    ),
}

#: Sign of every local node of the bilinear element (order: bottom-left, bottom-right,
#: top-left, top-right), the same as XN/YN in the original.
NODE_SIGN_X = np.array([-1.0, 1.0, -1.0, 1.0])
NODE_SIGN_Y = np.array([-1.0, -1.0, 1.0, 1.0])
