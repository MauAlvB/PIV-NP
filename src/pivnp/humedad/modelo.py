"""Del gris medido a la saturación y la humedad.

Tres pasos, los mismos del código MATLAB original:

1. **Referencias**: a partir del gris de una imagen de referencia se define, nodo a nodo, el
   gris del suelo seco y el del saturado, sumando y restando sendos desplazamientos.
2. **Normalización**: el gris de cada instante se lleva a una escala de 0 (saturado) a 100
   (seco), recortando los negativos a 0.
3. **Calibración**: ese valor normalizado se convierte en saturación y humedad con la curva
   del suelo.

Sobre esos tres pasos puede actuar la **política incremental**: el gris de un nodo no vuelve
a subir, y al pasar el umbral el nodo se da por saturado. Refleja un frente de humedecimiento
que avanza, y por eso **viene apagada**: es una hipótesis sobre el ensayo, no una medida, y
con ella puesta no se puede medir el secado.

Cuando se activa, actúa sobre el gris normalizado y no sobre la saturación, de modo que los
dos campos salen siempre del mismo valor. El código original la aplicaba solo a la saturación
y recalculaba la humedad entera en cada instante, con lo que los dos campos podían
contradecirse: en el caso de referencia acababa habiendo 634 nodos dados por saturados con la
humedad casi a cero.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np

from .calibracion import Calibracion

log = logging.getLogger("pivnp")


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


#: Por debajo de esta anchura, un nivel de gris pesa más del 5 % de la escala de saturación.
#: Con los 11 niveles del flujo RGB pesa el 9 %, y ahí el método se vuelve muy frágil.
ANCHO_MINIMO_DE_BANDA = 20.0


def avisar_si_la_banda_es_estrecha(referencias: Referencias, origen: str = "") -> float:
    """Comprueba la anchura de la banda, que es el número más sensible de todo el método.

    Devuelve la anchura mediana. Cuanto más estrecha, más pesa cada nivel de gris: medido
    sobre el caso del artículo, pasar de 40 niveles a 11 mueve la saturación 0,105 de media,
    un orden de magnitud más que cualquier otra decisión del cálculo.
    """
    ancho = referencias.seco - referencias.saturado
    finitos = ancho[np.isfinite(ancho)]
    if finitos.size == 0:
        return float("nan")
    mediana = float(np.median(finitos))
    if mediana < ANCHO_MINIMO_DE_BANDA:
        log.warning("%sla banda entre el suelo seco y el saturado son %.0f niveles de gris, "
                    "así que un nivel es el %.0f %% de la escala de saturación. Con una banda "
                    "tan estrecha el resultado depende mucho del ruido de la imagen; conviene "
                    "revisarla con el ensayo de calibración del suelo.",
                    f"{origen}: " if origen else "", mediana, 100.0 / mediana)
    return mediana


def referencias_globales(forma, gris_seco: float, gris_saturado: float) -> Referencias:
    """Referencias iguales en todos los nodos, medidas sobre el suelo del ensayo.

    Es lo que hace el flujo SWIR: en vez de sacar la referencia de la imagen nodo a nodo, se
    fijan dos intensidades para todo el material, la del suelo seco y la del saturado. La
    banda deja así de depender de la textura de cada punto.
    """
    if gris_seco <= gris_saturado:
        raise ValueError(f"el gris del suelo seco ({gris_seco}) debe ser mayor que el del "
                         f"saturado ({gris_saturado})")
    return Referencias(np.full(forma, float(gris_seco)), np.full(forma, float(gris_saturado)))


def normalizar(gris: np.ndarray, referencias: Referencias) -> np.ndarray:
    """Lleva el gris a la escala 0 (saturado) - 100 (seco), recortando los negativos."""
    recorrido = referencias.seco - referencias.saturado
    with np.errstate(divide="ignore", invalid="ignore"):
        normalizado = (np.asarray(gris, dtype=np.float64) - referencias.saturado) / recorrido
    normalizado = np.where(np.isfinite(normalizado), normalizado * 100.0, np.nan)
    return np.where(normalizado < 0.0, 0.0, normalizado)


#: Marcas de calidad de cada nodo, para poder leer los resultados sabiendo qué es medida.
MEDIDO = 0  #: el gris cae dentro de la banda: el valor es una medida
SIN_DATO = 1  #: el nodo no tiene dato (fuera de la imagen, o PIVlab no lo midió)
TOPE_HUMEDO = 2  #: el gris cae por debajo del extremo saturado: el valor es un "al menos"
TOPE_SECO = 3  #: el gris cae por encima del extremo seco: el valor es un "como mucho"


@dataclass
class Estado:
    """Resultado de un instante."""

    saturacion: np.ndarray
    humedad: np.ndarray
    #: Marca por nodo (:data:`MEDIDO`, :data:`SIN_DATO`, :data:`TOPE_HUMEDO`,
    #: :data:`TOPE_SECO`). Un nodo en un tope no es una medida, es una cota: conviene
    #: saberlo al leer un campo, porque suelen ser muchos.
    calidad: np.ndarray

    @property
    def recortados(self) -> int:
        """Nodos que caen fuera de la banda, y cuyo valor es por tanto una cota."""
        return int(np.isin(self.calidad, (TOPE_HUMEDO, TOPE_SECO)).sum())

    @property
    def medidos(self) -> int:
        return int((self.calidad == MEDIDO).sum())


class ModeloHumedad:
    """Convierte el gris de cada instante en saturación y humedad.

    La política incremental se aplica **sobre el gris normalizado**, no sobre la saturación:
    así los dos campos salen del mismo valor y no pueden contradecirse. El código original la
    aplicaba solo a la saturación y dejaba la humedad libre, de modo que un nodo podía quedar
    marcado como saturado con la humedad casi a cero.
    """

    def __init__(self, calibracion: Calibracion, umbral_saturacion: float = 0.95,
                 incremental: bool = False) -> None:
        if not 0 < umbral_saturacion <= 1:
            raise ValueError(f"el umbral debe estar en (0, 1] y vale {umbral_saturacion}")
        self.calibracion = calibracion
        self.umbral = float(umbral_saturacion)
        self.incremental = bool(incremental)
        #: Gris más bajo alcanzado por cada nodo: el suelo no se seca, así que no vuelve a subir.
        self.alcanzado: np.ndarray | None = None

    def reiniciar(self) -> None:
        self.alcanzado = None

    def evaluar(self, gris_normalizado: np.ndarray) -> Estado:
        """Saturación y humedad del instante, aplicando la política incremental."""
        gris = np.asarray(gris_normalizado, dtype=np.float64)
        medido = np.isfinite(gris)
        if self.incremental:
            gris = self._aplicar_trinquete(gris, medido)
        resultado = self.calibracion.evaluar(gris)
        return Estado(resultado.saturacion, resultado.humedad, self._calidad(gris, medido))

    # --- política incremental -------------------------------------------------------------
    def _aplicar_trinquete(self, gris: np.ndarray, medido: np.ndarray) -> np.ndarray:
        if self.alcanzado is None:
            self.alcanzado = np.full(gris.shape, np.inf)
        elif self.alcanzado.shape != gris.shape:
            raise ValueError("el número de nodos ha cambiado entre instantes")

        gris = np.where(medido, np.minimum(gris, self.alcanzado), gris)
        if self.umbral < 1.0:
            # Al pasar el umbral el nodo se da por saturado y ya no vuelve atrás.
            saturado = self.calibracion.evaluar(gris).saturacion >= self.umbral
            gris = np.where(medido & saturado, float(self.calibracion.rango[0]), gris)
        self.alcanzado = np.where(medido, gris, self.alcanzado)
        return gris

    def _calidad(self, gris: np.ndarray, medido: np.ndarray) -> np.ndarray:
        minimo, maximo = self.calibracion.rango
        calidad = np.full(gris.shape, SIN_DATO, dtype=np.int8)
        calidad[medido] = MEDIDO
        # el gris ya viene recortado por debajo en normalizar(), de ahí el <=
        calidad[medido & (gris <= minimo)] = TOPE_HUMEDO
        calidad[medido & (gris >= maximo)] = TOPE_SECO
        return calidad
