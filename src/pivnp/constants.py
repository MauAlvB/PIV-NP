"""Constantes numéricas heredadas del código Fortran original.

Varios literales del código original se escribieron sin sufijo ``d0`` (por ejemplo
``9.81`` o ``1.E-10``). Fortran los almacena en simple precisión (REAL*4) y luego los
promueve a REAL*8, lo que introduce un pequeño error de redondeo. Para que esta versión
reproduzca los resultados del original bit a bit, esos valores se guardan aquí con el
mismo redondeo. Ver ``docs/HALLAZGOS.md`` (H-07).
"""

from __future__ import annotations

import numpy as np


def as_fortran_real4(value: float) -> float:
    """Devuelve ``value`` redondeado a simple precisión, como un literal REAL de Fortran."""
    return float(np.float32(value))


#: ``EPSILON(1.d0)``: umbral para considerar nula la masa nodal o un incremento de deformación.
MACHINE_EPSILON: float = float(np.finfo(np.float64).eps)

#: Gravedad usada al inicializar la energía potencial (``9.81d0`` en PIVLAB_DATA).
GRAVITY_INITIAL: float = 9.81

#: Gravedad usada en cada paso (``9.81`` en SOLMOV, simple precisión).
GRAVITY_STEP: float = as_fortran_real4(9.81)

#: Umbral del segundo invariante J2 por debajo del cual q = 0 (``1.E-10`` en INVAR2).
J2_THRESHOLD: float = as_fortran_real4(1.0e-10)

#: Coordenadas locales (en [-1, 1]) de las partículas dentro de una celda, según NPC
#: (partículas por lado). NPC=2 y 3 son subdivisiones uniformes; NPC=4..6 son puntos de
#: Gauss-Legendre (ver H-06). En el original son literales REAL*4 (ver H-07).
PARTICLE_LOCAL_COORDS: dict[int, tuple[float, ...]] = {
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
