"""Orquestación: de las imágenes del ensayo a la humedad y la saturación de cada nodo.

Reúne las piezas anteriores y entrega, instante a instante, los dos campos que el análisis
necesita. Puede además escribir los archivos ``Moist_<n>.TXT``, que es el formato con el que
trabajaba el código MATLAB y sirve para comparar resultados.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .calibracion import Calibracion
from .configuracion import ConfiguracionHumedad, leer
from .imagenes import desenfoque_gaussiano, leer_imagen
from .modelo import (
    MEDIDO,
    SIN_DATO,
    TOPE_HUMEDO,
    TOPE_SECO,
    Estado,
    ModeloHumedad,
    Referencias,
    avisar_si_la_banda_es_estrecha,
    normalizar,
    referencias_globales,
    referencias_por_desplazamiento,
)
from .muestreo import coordenadas_en_pixeles, fuera_de_la_imagen, muestrear

CABECERA_MOIST = "x_m,y_m,moisture,saturation_degree"

log = logging.getLogger("pivnp")


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
        self.columna, self.fila = coordenadas_en_pixeles(malla.x_m, malla.y_m,
                                                         malla.metros_por_pixel,
                                                         configuracion.registro)
        self.referencias: Referencias | None = None
        self.nodos_fuera = 0
        #: Cuántos valores de cada marca de calidad se han publicado en todo el análisis.
        self.calidad: dict[int, int] = dict.fromkeys((MEDIDO, SIN_DATO, TOPE_HUMEDO,
                                                      TOPE_SECO), 0)

    # --- preparación ---------------------------------------------------------------------
    def preparar(self, con_dato: np.ndarray) -> Referencias:
        """Fija las referencias seca y saturada, que no cambian durante el ensayo.

        Con ``BANDA_SECA`` y ``BANDA_SATURADA`` en el ``.HUM`` son dos intensidades para toda
        la imagen y no hace falta imagen de referencia. Si no, se sacan de la imagen de
        referencia nodo a nodo, que es lo que hace el flujo RGB.
        """
        configuracion = self.configuracion
        if configuracion.banda_global:
            self.referencias = referencias_globales(self.columna.shape,
                                                    configuracion.banda_seca,
                                                    configuracion.banda_saturada)
        else:
            filtrada = self._imagen_filtrada(self.directorio / configuracion.referencia_seca)
            self.nodos_fuera = fuera_de_la_imagen(filtrada, self.columna, self.fila)
            gris = muestrear(filtrada, self.columna, self.fila, con_dato)
            self.referencias = referencias_por_desplazamiento(
                gris, configuracion.desplazamiento_seco,
                configuracion.desplazamiento_saturado)
        avisar_si_la_banda_es_estrecha(self.referencias, configuracion.origen)
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

    def resumen_de_calidad(self) -> str:
        """Cuánto de lo publicado es medida y cuánto es una cota, en todo el análisis."""
        total = sum(self.calidad.values())
        if not total:
            return "no se ha calculado ningún instante"
        con_dato = total - self.calidad[SIN_DATO]
        if not con_dato:
            return "ningún nodo con dato"
        partes = [f"{self.calidad[marca]} ({100 * self.calidad[marca] / con_dato:.0f} %) "
                  f"{nombre}"
                  for marca, nombre in ((MEDIDO, "medidos"),
                                        (TOPE_HUMEDO, "en el tope húmedo (son un 'al menos')"),
                                        (TOPE_SECO, "en el tope seco (son un 'como mucho')"))
                  if self.calidad[marca]]
        return f"de {con_dato} valores con dato: " + ", ".join(partes)

    def desde_gris(self, paso: int, normalizado: np.ndarray) -> Estado:
        """Humedad y saturación a partir del gris ya normalizado.

        Separado de la lectura de la imagen porque esta parte lleva memoria: la saturación no
        baja, así que los instantes tienen que pasar por aquí en orden. Lo de antes no la
        lleva y puede calcularse por adelantado en otro hilo.
        """
        estado = self.modelo.evaluar(normalizado)
        for marca in self.calidad:
            self.calidad[marca] += int((estado.calidad == marca).sum())
        if paso == 1 and self.configuracion.primer_instante == "legado":
            # El MATLAB dejaba la humedad a cero en el primer instante.
            estado = Estado(estado.saturacion, np.zeros_like(estado.humedad), estado.calidad)
        return estado

    def instante(self, paso: int, con_dato: np.ndarray) -> Estado:
        """Humedad y saturación de un instante."""
        return self.desde_gris(paso, self.gris_normalizado(paso, con_dato))

    def reiniciar(self) -> None:
        self.modelo.reiniciar()


def fuente_de_caso(case_dir: Path, case_name: str, malla: Malla,
                   con_dato: np.ndarray) -> FuenteHumedad:
    """Prepara la fuente de humedad de un caso a partir de su ``<caso>.HUM``.

    El archivo de calibración se busca junto al del caso, salvo que se dé una ruta absoluta.
    """
    from ..config import find_file  # aquí para no crear una dependencia circular

    case_dir = Path(case_dir)
    configuracion = leer(find_file(case_dir, f"{case_name}.HUM"))
    if configuracion.desconocidas:
        log.warning("%s: claves que no se reconocen y se ignoran: %s", configuracion.origen,
                    ", ".join(configuracion.desconocidas))
    ruta_calibracion = Path(configuracion.calibracion)
    if not ruta_calibracion.is_absolute():
        ruta_calibracion = find_file(case_dir, configuracion.calibracion)
    fuente = FuenteHumedad(case_dir, configuracion, Calibracion.desde_csv(ruta_calibracion),
                           malla)
    fuente.preparar(con_dato)
    if fuente.nodos_fuera:
        log.warning("%d nodos de la malla caen fuera de la imagen de humedad: revisa el "
                    "registro (ESCALA_X, ORIGEN_X, ESCALA_Y, ORIGEN_Y) del .HUM",
                    fuente.nodos_fuera)
    log.info("Humedad desde las imágenes: %s, canal %d, sigma %g, calibración %s",
             configuracion.patron_imagenes, configuracion.canal, configuracion.sigma,
             Path(ruta_calibracion).name)
    return fuente


def escribir_moist(ruta: Path, malla: Malla, humedad: np.ndarray,
                   saturacion: np.ndarray) -> None:
    """Escribe un archivo ``Moist_<n>.TXT`` con el formato del código MATLAB."""
    lineas = [CABECERA_MOIST]
    for x, y, h, s in zip(malla.x_m, malla.y_m, humedad, saturacion, strict=True):
        lineas.append(f"{_numero(x, 10)},{_numero(y, 10)},"
                      f"{_numero(h, 17)},{_numero(s, 15)}")
    Path(ruta).write_text("\n".join(lineas) + "\n", encoding="ascii")


def _numero(valor: float, cifras: int) -> str:
    """Formato de MATLAB: ``NaN`` para lo que falta y notación breve para el resto."""
    if not np.isfinite(valor):
        return "NaN"
    return f"{float(valor):.{cifras}g}"
