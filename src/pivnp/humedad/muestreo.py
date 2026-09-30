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

    Con una sola cámara es la identidad. Con dos (por ejemplo visible e infrarroja) es la
    homografía que hace coincidir ambas imágenes, la misma que el código MATLAB construía
    marcando puntos a mano sobre las dos::

        | columna |   | escala_x       inclinacion_xy  origen_x |   | x |
        | fila    | ~ | inclinacion_yx escala_y        origen_y | · | y |
        | 1       |   | perspectiva_x  perspectiva_y   1        |   | 1 |

    Con las dos cámaras en el mismo sitio basta la escala y el origen; los otros cuatro
    valores hacen falta cuando miran desde ángulos distintos.
    """

    escala_x: float = 1.0
    origen_x: float = 0.0
    escala_y: float = 1.0
    origen_y: float = 0.0
    inclinacion_xy: float = 0.0
    inclinacion_yx: float = 0.0
    perspectiva_x: float = 0.0
    perspectiva_y: float = 0.0

    @property
    def es_identidad(self) -> bool:
        return (self.escala_x, self.escala_y) == (1.0, 1.0) and not any(
            (self.origen_x, self.origen_y, self.inclinacion_xy, self.inclinacion_yx,
             self.perspectiva_x, self.perspectiva_y))


SIN_REGISTRO = Registro()


def coordenadas_en_pixeles(x_m: np.ndarray, y_m: np.ndarray, metros_por_pixel: float,
                           registro: Registro = SIN_REGISTRO) -> tuple[np.ndarray, np.ndarray]:
    """Columna y fila (base 0) de cada nodo dentro de la imagen.

    Se pasa de metros a píxeles, se aplica el registro y se redondea al píxel más cercano,
    como hacía el MATLAB; se resta 1 porque allí los índices empiezan en 1.
    """
    if metros_por_pixel <= 0:
        raise ValueError(f"el factor de conversión debe ser positivo y vale {metros_por_pixel}")
    x = np.asarray(x_m, dtype=np.float64) / metros_por_pixel
    y = np.asarray(y_m, dtype=np.float64) / metros_por_pixel
    if not registro.es_identidad:
        peso = registro.perspectiva_x * x + registro.perspectiva_y * y + 1.0
        peso = np.where(np.abs(peso) < 1e-12, np.nan, peso)
        x, y = ((registro.escala_x * x + registro.inclinacion_xy * y + registro.origen_x) / peso,
                (registro.inclinacion_yx * x + registro.escala_y * y + registro.origen_y) / peso)
    columna = np.floor(np.nan_to_num(x, nan=-1e9) + 0.5)
    fila = np.floor(np.nan_to_num(y, nan=-1e9) + 0.5)
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
