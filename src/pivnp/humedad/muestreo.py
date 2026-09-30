"""Del sistema de la malla PIV al píxel de la imagen, y muestreo del gris.

Los nodos de la malla vienen en metros en los archivos de PIVlab; el factor de conversión
está en la cabecera de esos mismos archivos. El código MATLAB redondeaba las coordenadas al
píxel más cercano, y aquí se hace igual para que los valores coincidan.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Registro:
    """Transformación de la malla PIV a la imagen de humedad.

    Con una sola cámara es la identidad. Con dos (por ejemplo visible e infrarroja) son la
    escala y el desplazamiento que hacen coincidir ambas imágenes.
    """

    escala_x: float = 1.0
    origen_x: float = 0.0
    escala_y: float = 1.0
    origen_y: float = 0.0


SIN_REGISTRO = Registro()


def coordenadas_en_pixeles(x_m: np.ndarray, y_m: np.ndarray, metros_por_pixel: float,
                           registro: Registro = SIN_REGISTRO) -> tuple[np.ndarray, np.ndarray]:
    """Columna y fila (base 0) de cada nodo dentro de la imagen.

    Se redondea al píxel más cercano, como hacía el MATLAB, y se resta 1 porque allí los
    índices empiezan en 1.
    """
    if metros_por_pixel <= 0:
        raise ValueError(f"el factor de conversión debe ser positivo y vale {metros_por_pixel}")
    columna = np.floor(np.asarray(x_m, dtype=np.float64) / metros_por_pixel + 0.5)
    fila = np.floor(np.asarray(y_m, dtype=np.float64) / metros_por_pixel + 0.5)
    columna = np.floor(columna * registro.escala_x + registro.origen_x + 0.5)
    fila = np.floor(fila * registro.escala_y + registro.origen_y + 0.5)
    return columna.astype(np.int64) - 1, fila.astype(np.int64) - 1


def muestrear(imagen: np.ndarray, columna: np.ndarray, fila: np.ndarray,
              con_dato: np.ndarray | None = None) -> np.ndarray:
    """Gris de la imagen en cada nodo; NaN donde no hay dato o el nodo cae fuera.

    ``con_dato`` marca los nodos que PIVlab sí midió; el resto se descarta, igual que hacía
    el MATLAB con la máscara del archivo ``.mat``.
    """
    alto, ancho = imagen.shape
    dentro = (columna >= 0) & (columna < ancho) & (fila >= 0) & (fila < alto)
    utiles = dentro if con_dato is None else dentro & np.asarray(con_dato, dtype=bool)

    gris = np.full(columna.shape, np.nan)
    gris[utiles] = imagen[fila[utiles], columna[utiles]]
    return gris


def fuera_de_la_imagen(imagen: np.ndarray, columna: np.ndarray, fila: np.ndarray) -> int:
    """Cuántos nodos caen fuera de la imagen (indica un registro mal ajustado)."""
    alto, ancho = imagen.shape
    return int((~((columna >= 0) & (columna < ancho) & (fila >= 0) & (fila < alto))).sum())
