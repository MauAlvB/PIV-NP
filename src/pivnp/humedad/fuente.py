"""Orquestación: de las imágenes del ensayo a la humedad y la saturación de cada nodo.

Reúne las piezas anteriores y entrega, instante a instante, los dos campos que el análisis
necesita. Puede además escribir los archivos ``Moist_<n>.TXT``, que es el formato con el que
trabajaba el código MATLAB y sirve para comparar resultados.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .calibracion import Calibracion
from .configuracion import ConfiguracionHumedad
from .imagenes import desenfoque_gaussiano, leer_imagen
from .modelo import Estado, ModeloHumedad, Referencias, normalizar, referencias_por_desplazamiento
from .muestreo import Registro, coordenadas_en_pixeles, fuera_de_la_imagen, muestrear

CABECERA_MOIST = "x_m,y_m,moisture,saturation_degree"


@dataclass(frozen=True)
class Malla:
    """Nodos de la malla PIV, en metros, y el factor de conversión a píxeles."""

    x_m: np.ndarray
    y_m: np.ndarray
    metros_por_pixel: float

    def __post_init__(self) -> None:
        if self.x_m.shape != self.y_m.shape:
            raise ValueError("x e y deben tener el mismo número de nodos")
        if self.metros_por_pixel <= 0:
            raise ValueError("el factor de conversión debe ser positivo")


class FuenteHumedad:
    """Calcula humedad y saturación desde las imágenes, instante a instante."""

    def __init__(self, directorio: Path, configuracion: ConfiguracionHumedad,
                 calibracion: Calibracion, malla: Malla) -> None:
        self.directorio = Path(directorio)
        self.configuracion = configuracion
        self.malla = malla
        self.modelo = ModeloHumedad(calibracion, configuracion.umbral_saturacion,
                                    configuracion.incremental)
        registro = Registro(configuracion.escala_x, configuracion.origen_x,
                            configuracion.escala_y, configuracion.origen_y)
        self.columna, self.fila = coordenadas_en_pixeles(malla.x_m, malla.y_m,
                                                         malla.metros_por_pixel, registro)
        self.referencias: Referencias | None = None
        self.nodos_fuera = 0

    # --- preparación ---------------------------------------------------------------------
    def preparar(self, con_dato: np.ndarray) -> Referencias:
        """Calcula las referencias seca y saturada a partir de la imagen de referencia."""
        ruta = self.directorio / self.configuracion.referencia_seca
        filtrada = self._imagen_filtrada(ruta)
        self.nodos_fuera = fuera_de_la_imagen(filtrada, self.columna, self.fila)
        gris = muestrear(filtrada, self.columna, self.fila, con_dato)
        self.referencias = referencias_por_desplazamiento(
            gris, self.configuracion.desplazamiento_seco,
            self.configuracion.desplazamiento_saturado)
        return self.referencias

    # --- cálculo -------------------------------------------------------------------------
    def gris_normalizado(self, paso: int, con_dato: np.ndarray) -> np.ndarray:
        """Gris del instante, en la escala 0 (saturado) - 100 (seco)."""
        if self.referencias is None:
            raise RuntimeError("hay que llamar a preparar() antes del primer instante")
        filtrada = self._imagen_filtrada(self.configuracion.ruta_imagen(self.directorio, paso))
        gris = muestrear(filtrada, self.columna, self.fila, con_dato)
        return normalizar(gris, self.referencias)

    def _imagen_filtrada(self, ruta: Path) -> np.ndarray:
        imagen = leer_imagen(ruta, self.configuracion.canal)
        return desenfoque_gaussiano(imagen, self.configuracion.sigma,
                                    self.configuracion.redondeo_legado)

    def instante(self, paso: int, con_dato: np.ndarray) -> Estado:
        """Humedad y saturación de un instante."""
        estado = self.modelo.evaluar(self.gris_normalizado(paso, con_dato))
        if paso == 1 and self.configuracion.primer_instante == "legado":
            # El MATLAB dejaba la humedad a cero en el primer instante.
            estado = Estado(estado.saturacion, np.zeros_like(estado.humedad),
                            estado.recortados)
        return estado

    def reiniciar(self) -> None:
        self.modelo.reiniciar()


def escribir_moist(ruta: Path, malla: Malla, estado: Estado) -> None:
    """Escribe un archivo ``Moist_<n>.TXT`` con el formato del código MATLAB."""
    lineas = [CABECERA_MOIST]
    for x, y, humedad, saturacion in zip(malla.x_m, malla.y_m, estado.humedad,
                                         estado.saturacion, strict=True):
        lineas.append(f"{_numero(x, 10)},{_numero(y, 10)},"
                      f"{_numero(humedad, 17)},{_numero(saturacion, 15)}")
    Path(ruta).write_text("\n".join(lineas) + "\n", encoding="ascii")


def _numero(valor: float, cifras: int) -> str:
    """Formato de MATLAB: ``NaN`` para lo que falta y notación breve para el resto."""
    if not np.isfinite(valor):
        return "NaN"
    return f"{float(valor):.{cifras}g}"
