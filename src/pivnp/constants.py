"""Constantes numéricas.

Varios literales del código Fortran original se escribieron sin sufijo ``d0`` (por ejemplo
``9.81`` o ``1.E-10``), así que Fortran los almacenaba en simple precisión (REAL*4) antes de
promoverlos a REAL*8, perdiendo precisión a partir de la séptima cifra. Aquí se usan
en doble precisión; las versiones redondeadas se conservan con el sufijo ``LEGACY_`` para el
modo compatibilidad.
"""

from __future__ import annotations

import numpy as np


def as_fortran_real4(value: float) -> float:
    """Devuelve ``value`` redondeado a simple precisión, como un literal REAL de Fortran."""
    return float(np.float32(value))


#: ``EPSILON(1.d0)``: umbral para considerar nula la masa nodal o un incremento de deformación.
MACHINE_EPSILON: float = float(np.finfo(np.float64).eps)

#: Gravedad [m/s²] para la energía potencial.
GRAVITY: float = 9.81
LEGACY_GRAVITY_STEP: float = as_fortran_real4(9.81)  # ``9.81`` (REAL*4) en SOLMOV

#: Umbral del segundo invariante J2 por debajo del cual q = 0 (``1.E-10`` en INVAR2).
J2_THRESHOLD: float = 1.0e-10
LEGACY_J2_THRESHOLD: float = as_fortran_real4(1.0e-10)

#: Coordenadas locales (en [-1, 1]) de las NPC filas/columnas de partículas de cada celda.
#: NPC = 2 y 3 son subdivisiones uniformes (centros de las subceldas); NPC = 4, 5 y 6 son
#: puntos de Gauss-Legendre.
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

#: Las mismas coordenadas tal como las guardaba el original (literales REAL*4).
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

#: Signo de cada nodo local del elemento bilineal (orden: abajo-izq, abajo-der, arriba-izq,
#: arriba-der), igual que XN/YN en el original.
NODE_SIGN_X = np.array([-1.0, 1.0, -1.0, 1.0])
NODE_SIGN_Y = np.array([-1.0, -1.0, 1.0, 1.0])
