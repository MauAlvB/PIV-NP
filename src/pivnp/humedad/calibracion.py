"""Curva de calibración: del gris normalizado a la saturación y la humedad.

La curva sale de ensayos de laboratorio y cambia con cada suelo, así que es un dato de
entrada más. El archivo es un CSV con tres columnas::

    gris,saturacion,humedad
    100,0.0079,0.0890
    99.5,0.0074,0.1084
    ...

* ``gris``: valor de gris normalizado entre la referencia saturada (0) y la seca (100).
* ``saturacion``: grado de saturación, entre 0 y 1.
* ``humedad``: humedad en %.

Entre los puntos de la tabla se interpola igual que el código MATLAB original (``pchip``),
de modo que los resultados son comparables con los análisis anteriores. Fuera del rango de
la tabla **no se extrapola**: se recorta al extremo más cercano y se cuentan los valores
recortados, porque extrapolar era el origen de humedades imposibles (del orden de −1700 %).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


def _pendientes_pchip(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Pendientes de Fritsch-Carlson, las mismas que usa ``pchip`` de MATLAB."""
    h = np.diff(x)
    delta = np.diff(y) / h
    d = np.zeros_like(y)

    interior = slice(1, len(y) - 1)
    mismo_signo = delta[:-1] * delta[1:] > 0
    w1 = 2 * h[1:] + h[:-1]
    w2 = h[1:] + 2 * h[:-1]
    with np.errstate(divide="ignore", invalid="ignore"):
        media = (w1 + w2) / (w1 / delta[:-1] + w2 / delta[1:])
    d[interior] = np.where(mismo_signo, media, 0.0)

    for extremo, (i, j) in ((0, (0, 1)), (len(y) - 1, (-1, -2))):
        pendiente = ((2 * h[i] + h[j]) * delta[i] - h[i] * delta[j]) / (h[i] + h[j])
        if pendiente * delta[i] <= 0:
            pendiente = 0.0
        elif delta[i] * delta[j] < 0 and abs(pendiente) > abs(3 * delta[i]):
            pendiente = 3 * delta[i]
        d[extremo] = pendiente
    return d


def interpolar_pchip(x: np.ndarray, y: np.ndarray, consulta: np.ndarray) -> np.ndarray:
    """Interpolación cúbica monótona por tramos, compatible con ``interp1(...,'pchip')``.

    ``x`` debe estar ordenado de menor a mayor. Los valores de ``consulta`` fuera del rango
    se evalúan con el tramo del extremo correspondiente (quien quiera recortar debe hacerlo
    antes; :class:`Calibracion` lo hace).
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    consulta = np.asarray(consulta, dtype=np.float64)
    if x.ndim != 1 or x.size < 2 or y.shape != x.shape:
        raise ValueError("x e y deben ser vectores del mismo tamaño y con 2 puntos o más")
    if np.any(np.diff(x) <= 0):
        raise ValueError("x debe estar ordenado de menor a mayor y sin valores repetidos")

    d = _pendientes_pchip(x, y)
    i = np.clip(np.searchsorted(x, consulta) - 1, 0, x.size - 2)
    h = x[i + 1] - x[i]
    t = (consulta - x[i]) / h
    t2, t3 = t * t, t * t * t
    return ((2 * t3 - 3 * t2 + 1) * y[i]
            + (t3 - 2 * t2 + t) * h * d[i]
            + (-2 * t3 + 3 * t2) * y[i + 1]
            + (t3 - t2) * h * d[i + 1])


@dataclass(frozen=True)
class Resultado:
    """Saturación y humedad de un conjunto de valores de gris."""

    saturacion: np.ndarray
    humedad: np.ndarray
    recortados: int  # valores que caían fuera de la tabla de calibración


@dataclass(frozen=True)
class Calibracion:
    """Tabla de calibración de un suelo, ordenada por gris creciente."""

    gris: np.ndarray
    saturacion: np.ndarray
    humedad: np.ndarray
    origen: str = "<memoria>"

    def __post_init__(self) -> None:
        if self.gris.size < 2:
            raise ValueError(f"{self.origen}: la calibración necesita al menos 2 puntos")
        if not (self.gris.shape == self.saturacion.shape == self.humedad.shape):
            raise ValueError(f"{self.origen}: las tres columnas deben tener el mismo tamaño")
        if np.any(np.diff(self.gris) <= 0):
            raise ValueError(f"{self.origen}: la columna de gris debe crecer sin repeticiones")
        if np.any(~np.isfinite(self.gris)) or np.any(~np.isfinite(self.saturacion)) \
                or np.any(~np.isfinite(self.humedad)):
            raise ValueError(f"{self.origen}: la calibración no admite valores vacíos o NaN")

    @property
    def rango(self) -> tuple[float, float]:
        return float(self.gris[0]), float(self.gris[-1])

    @classmethod
    def desde_csv(cls, ruta: Path) -> Calibracion:
        """Lee la calibración de un CSV con columnas ``gris,saturacion,humedad``.

        Se ignoran las líneas en blanco y las que empiezan por ``#``, y se admite que la
        tabla venga en cualquier orden.
        """
        ruta = Path(ruta)
        filas = []
        for numero, linea in enumerate(ruta.read_text(encoding="utf-8-sig").splitlines(), 1):
            limpia = linea.strip()
            if not limpia or limpia.startswith("#"):
                continue
            partes = [p.strip() for p in limpia.replace(";", ",").split(",")]
            if len(partes) < 3:
                raise ValueError(f"{ruta}:{numero}: se esperaban 3 columnas y hay {len(partes)}")
            if numero == 1 or not _es_numero(partes[0]):
                continue  # cabecera
            try:
                filas.append([float(p) for p in partes[:3]])
            except ValueError:
                raise ValueError(f"{ruta}:{numero}: valores no numéricos ({limpia!r})") from None
        if not filas:
            raise ValueError(f"{ruta}: no contiene ninguna fila de datos")

        tabla = np.array(filas, dtype=np.float64)
        orden = np.argsort(tabla[:, 0])
        return cls(tabla[orden, 0], tabla[orden, 1], tabla[orden, 2], origen=str(ruta))

    def evaluar(self, gris: np.ndarray) -> Resultado:
        """Saturación y humedad de cada valor de gris, recortando fuera de la tabla.

        Los valores no finitos (nodos sin dato) se propagan como NaN y no cuentan como
        recortados.
        """
        gris = np.asarray(gris, dtype=np.float64)
        validos = np.isfinite(gris)
        minimo, maximo = self.rango
        fuera = validos & ((gris < minimo) | (gris > maximo))
        acotado = np.clip(gris, minimo, maximo)

        saturacion = np.full(gris.shape, np.nan)
        humedad = np.full(gris.shape, np.nan)
        if validos.any():
            saturacion[validos] = interpolar_pchip(self.gris, self.saturacion, acotado[validos])
            humedad[validos] = interpolar_pchip(self.gris, self.humedad, acotado[validos])
        return Resultado(saturacion, humedad, int(fuera.sum()))


def _es_numero(texto: str) -> bool:
    try:
        float(texto)
    except ValueError:
        return False
    return True
