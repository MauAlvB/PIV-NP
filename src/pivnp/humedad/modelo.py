"""Del gris medido a la saturación y la humedad.

Tres pasos, los mismos del código MATLAB original:

1. **Referencias**: a partir del gris de una imagen de referencia se define, nodo a nodo, el
   gris del suelo seco y el del saturado, sumando y restando sendos desplazamientos.
2. **Normalización**: el gris de cada instante se lleva a una escala de 0 (saturado) a 100
   (seco), recortando los negativos a 0.
3. **Calibración**: ese valor normalizado se convierte en saturación y humedad con la curva
   del suelo.

Sobre esos tres pasos actúa la **política incremental**: la saturación de un nodo no puede
bajar respecto al instante anterior, y al superar el umbral se fija en 1. Refleja un frente
de humedecimiento que avanza; para medir también el secado habrá que desactivarla.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .calibracion import Calibracion


@dataclass(frozen=True)
class Referencias:
    """Gris del suelo seco y del saturado en cada nodo."""

    seco: np.ndarray
    saturado: np.ndarray

    def __post_init__(self) -> None:
        if self.seco.shape != self.saturado.shape:
            raise ValueError("las dos referencias deben tener el mismo tamaño")


def referencias_por_desplazamiento(gris_referencia: np.ndarray, desplazamiento_seco: float,
                                   desplazamiento_saturado: float) -> Referencias:
    """Referencias a partir de una sola imagen, sumando y restando un valor fijo.

    Es lo que hacía el MATLAB (+5 y −6 niveles de gris). Los dos números fijan toda la
    escala de saturación, así que conviene revisarlos con datos de laboratorio.
    """
    if desplazamiento_seco <= desplazamiento_saturado:
        raise ValueError("el desplazamiento seco debe ser mayor que el saturado")
    gris = np.asarray(gris_referencia, dtype=np.float64)
    return Referencias(gris + desplazamiento_seco, gris + desplazamiento_saturado)


def normalizar(gris: np.ndarray, referencias: Referencias) -> np.ndarray:
    """Lleva el gris a la escala 0 (saturado) - 100 (seco), recortando los negativos."""
    recorrido = referencias.seco - referencias.saturado
    with np.errstate(divide="ignore", invalid="ignore"):
        normalizado = (np.asarray(gris, dtype=np.float64) - referencias.saturado) / recorrido
    normalizado = np.where(np.isfinite(normalizado), normalizado * 100.0, np.nan)
    return np.where(normalizado < 0.0, 0.0, normalizado)


@dataclass
class Estado:
    """Resultado de un instante."""

    saturacion: np.ndarray
    humedad: np.ndarray
    recortados: int  # nodos fuera del rango de la calibración


class ModeloHumedad:
    """Convierte el gris de cada instante en saturación y humedad.

    Guarda la saturación alcanzada por cada nodo, que es lo que permite la política
    incremental.
    """

    def __init__(self, calibracion: Calibracion, umbral_saturacion: float = 0.8,
                 incremental: bool = True) -> None:
        if not 0 < umbral_saturacion <= 1:
            raise ValueError(f"el umbral debe estar en (0, 1] y vale {umbral_saturacion}")
        self.calibracion = calibracion
        self.umbral = float(umbral_saturacion)
        self.incremental = bool(incremental)
        self.alcanzada: np.ndarray | None = None

    def reiniciar(self) -> None:
        self.alcanzada = None

    def evaluar(self, gris_normalizado: np.ndarray) -> Estado:
        """Saturación y humedad del instante, aplicando la política incremental."""
        resultado = self.calibracion.evaluar(gris_normalizado)
        saturacion = resultado.saturacion.copy()

        if not self.incremental:
            return Estado(saturacion, resultado.humedad, resultado.recortados)

        if self.alcanzada is None:
            self.alcanzada = np.full(saturacion.shape, np.nan)
        elif self.alcanzada.shape != saturacion.shape:
            raise ValueError("el número de nodos ha cambiado entre instantes")

        medido = np.isfinite(saturacion)
        previa = np.where(np.isfinite(self.alcanzada), self.alcanzada, -np.inf)
        # al superar el umbral el nodo se da por saturado y ya no vuelve atrás
        nueva = np.where(saturacion >= self.umbral, 1.0, np.maximum(saturacion, previa))
        self.alcanzada = np.where(medido, nueva, self.alcanzada)
        return Estado(np.where(medido, nueva, np.nan), resultado.humedad,
                      resultado.recortados)
