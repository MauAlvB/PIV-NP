"""Traducción literal (Python puro, índices base 1) de fragmentos del Fortran original.

Sirve como oráculo en las pruebas unitarias: las funciones optimizadas deben dar lo mismo
que estos bucles, escritos tal como están en ``legacy/MainCodePIV-NP.for``.
"""

from __future__ import annotations

import math

import numpy as np


def mesh_arrays(nc: int, nfil: int, axc: float, ayc: float, staggered: bool = False):
    """NCF, NNF, NNFA, XF, YF de PIVLAB_DATA (malla 1 o malla desplazada 2)."""
    nch = nc // nfil
    if staggered:
        nfil, nch = nfil + 1, nch + 1
        nc = nch * nfil
        xf0, yf0 = -axc / 2, -ayc / 2
    else:
        xf0, yf0 = 0.0, 0.0
    ncf = [0] * (nfil + 2)
    nnf = [0] * (nfil + 1)
    nnfa = [0] * (nfil + 1)
    ncf[1], nnf[1], nnfa[1] = 1, 1, 1 + nch + 1
    for i in range(2, nfil + 1):
        ncf[i] = ncf[i - 1] + nch
        nnf[i] = nnf[i - 1] + nch + 1
        nnfa[i] = nnfa[i - 1] + nch + 1
    ncf[nfil + 1] = nc + 1
    yf = [0.0] + [yf0 + (i - 1) * ayc for i in range(1, nfil + 2)]
    return nch, nfil, ncf, nnf, nnfa, xf0, yf


def ucelda(x: float, y: float, arrays, axc: float):
    """Búsqueda lineal de UCELDA: devuelve (INDC, IN(1:4)) en base 1, o None."""
    nch, nfil, ncf, nnf, nnfa, xf0, yf = arrays
    ifi = 0
    for j in range(1, nfil + 1):
        if yf[j] <= y < yf[j + 1]:
            ifi = j
            break
    if ifi == 0:
        return None
    for ic in range(ncf[ifi], ncf[ifi + 1]):
        xc = xf0 + axc * (ic - ncf[ifi])
        if xc <= x < xc + axc:
            n1 = nnf[ifi] + ic - ncf[ifi]
            n3 = nnfa[ifi] + ic - ncf[ifi]
            return ic, (n1, n1 + 1, n3, n3 + 1)
    return None


def generate_particles(nc: int, nfil: int, npc: int, axc: float, ayc: float, gauss,
                       staggered: bool = False) -> np.ndarray:
    """Bucle de generación de partículas de PIVLAB_DATA."""
    nch, nfil, ncf, _, _, xf0, yf = mesh_arrays(nc, nfil, axc, ayc, staggered)
    axp, ayp = axc / npc, ayc / npc
    xp = []
    for i in range(1, nfil + 1):
        for j in range(ncf[i], ncf[i + 1]):
            for ii in range(1, npc + 1):
                for jj in range(1, npc + 1):
                    if npc <= 10:
                        x = xf0 + (j - ncf[i]) * axc + axc / 2.0 + gauss[jj - 1] * axc / 2.0
                        y = yf[i] + ayc / 2.0 + gauss[ii - 1] * ayc / 2.0
                    else:
                        x = xf0 + (j - ncf[i]) * axc + (jj - 1) * axp + axp / 2.0
                        y = yf[i] + (ii - 1) * ayp + ayp / 2.0
                    xp.append((x, y))
    return np.array(xp)


def iconectividad(nch: int, nfil: int) -> dict[int, int]:
    """ICONECTIVIDAD(punto PIVlab) = nodo PIV-NP, base 1."""
    conn = {}
    for i in range(1, nch + 2):
        for j in range(1, nfil + 2):
            conn[i * (nfil + 2 - j) + (j - 1) * (i - 1)] = i + (j - 1) * (nch + 1)
    return conn


def invar2(sigx: float, sigy: float, sigz: float, sigxy: float, threshold: float) -> float:
    sigm = (sigx + sigy + sigz) / 3.0
    sx, sy, sz = sigx - sigm, sigy - sigm, sigz - sigm
    rj2 = (sx * sx + sy * sy + sz * sz) / 2.0 + sigxy * sigxy
    return math.sqrt(3.0 * rj2) if rj2 > threshold else 0.0
