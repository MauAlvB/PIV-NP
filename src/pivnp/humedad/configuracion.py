"""Archivo de configuración de la medición de humedad (``<caso>.HUM``).

Cada valor va identificado por su nombre, el orden no importa y los comentarios empiezan por
``!`` como en Fortran::

    ! Configuración de la medición de humedad
    IMAGENES   = vis_{n}.jpg      ! {n} se sustituye por el número de instante
    CANAL      = gris             ! 1 rojo, 2 verde, 3 azul, 0 o "gris" para escala de grises
    SIGMA      = 40               ! radio de promediado, en píxeles
    REFERENCIA_SECA = ref2.jpg
    CALIBRACION = calibracion_slope_rgb.csv

Los valores que no aparezcan toman el valor por defecto, que es el del código MATLAB
original, para que los análisis sean comparables con los anteriores. Las dos excepciones son
``PRIMER_INSTANTE`` y ``REDONDEO_LEGADO``: el criterio nuevo es el de por defecto y el
antiguo se recupera poniéndolas a ``legado`` y a ``1``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from .imagenes import numero_de_canal
from .muestreo import Registro

#: Valor por defecto de cada clave. Son los del código MATLAB original salvo donde se acordó
#: cambiar de criterio: CANAL (el análisis de referencia se hizo en escala de grises, no con
#: el canal rojo), PRIMER_INSTANTE, REDONDEO_LEGADO, e INCREMENTAL con su UMBRAL_SATURACION,
#: que pasan a estar apagados porque son una hipótesis sobre el ensayo y no una medida.
PREDETERMINADOS: dict[str, str] = {
    "IMAGENES": "vis_{n}.jpg",
    "CANAL": "gris",
    "SIGMA": "40",
    "REFERENCIA_SECA": "ref2.jpg",
    "DESPLAZAMIENTO_SECO": "5",
    "DESPLAZAMIENTO_SATURADO": "-6",
    "BANDA_SECA": "",
    "BANDA_SATURADA": "",
    "CALIBRACION": "",
    "UMBRAL_SATURACION": "0.95",
    "INCREMENTAL": "0",
    "PRIMER_INSTANTE": "igual",
    "REDONDEO_LEGADO": "0",
    "ESCALA_X": "1.0",
    "ORIGEN_X": "0.0",
    "ESCALA_Y": "1.0",
    "ORIGEN_Y": "0.0",
    "INCLINACION_XY": "0.0",
    "INCLINACION_YX": "0.0",
    "PERSPECTIVA_X": "0.0",
    "PERSPECTIVA_Y": "0.0",
}

#: Qué hacer con el primer instante: como los demás, o como el MATLAB (humedad 0 y
#: saturación constante).
PRIMER_INSTANTE = ("igual", "legado")


class ConfiguracionError(ValueError):
    """El archivo de configuración de humedad falta o es incoherente."""


@dataclass(frozen=True)
class ConfiguracionHumedad:
    """Parámetros de la medición de humedad de un ensayo."""

    patron_imagenes: str
    canal: int
    sigma: float
    referencia_seca: str
    desplazamiento_seco: float
    desplazamiento_saturado: float
    #: Intensidades del suelo seco y del saturado, iguales para toda la imagen. Si se dan,
    #: sustituyen a la referencia por nodo; es lo que hace el flujo SWIR.
    banda_seca: float | None
    banda_saturada: float | None
    calibracion: str
    umbral_saturacion: float
    incremental: bool
    primer_instante: str
    redondeo_legado: bool
    escala_x: float
    origen_x: float
    escala_y: float
    origen_y: float
    inclinacion_xy: float = 0.0
    inclinacion_yx: float = 0.0
    perspectiva_x: float = 0.0
    perspectiva_y: float = 0.0
    origen: str = "<memoria>"
    desconocidas: tuple[str, ...] = field(default_factory=tuple)

    def ruta_imagen(self, directorio: Path, paso: int) -> Path:
        return Path(directorio) / self.patron_imagenes.format(n=paso)

    @property
    def banda_global(self) -> bool:
        """Si la banda seca-saturada se fija con dos intensidades para toda la imagen."""
        return self.banda_seca is not None and self.banda_saturada is not None

    @property
    def registro(self) -> Registro:
        """Transformación de la malla PIV a la imagen de humedad."""
        return Registro(self.escala_x, self.origen_x, self.escala_y, self.origen_y,
                        self.inclinacion_xy, self.inclinacion_yx,
                        self.perspectiva_x, self.perspectiva_y)

    @property
    def registro_es_identidad(self) -> bool:
        return self.registro.es_identidad

    def validar(self) -> None:
        errores = []
        if "{n}" not in self.patron_imagenes:
            errores.append("IMAGENES debe contener {n}, que se sustituye por el instante")
        if self.sigma <= 0:
            errores.append(f"SIGMA={self.sigma} debe ser positivo")
        if not 0 < self.umbral_saturacion <= 1:
            errores.append(f"UMBRAL_SATURACION={self.umbral_saturacion} debe estar en (0, 1]")
        if self.desplazamiento_seco <= self.desplazamiento_saturado:
            errores.append("DESPLAZAMIENTO_SECO debe ser mayor que DESPLAZAMIENTO_SATURADO")
        if (self.banda_seca is None) != (self.banda_saturada is None):
            errores.append("BANDA_SECA y BANDA_SATURADA van juntas: o las dos o ninguna")
        elif self.banda_global and self.banda_seca <= self.banda_saturada:
            errores.append(f"BANDA_SECA={self.banda_seca} debe ser mayor que "
                           f"BANDA_SATURADA={self.banda_saturada}")
        if self.primer_instante not in PRIMER_INSTANTE:
            errores.append(f"PRIMER_INSTANTE={self.primer_instante!r} debe ser "
                           f"{' o '.join(PRIMER_INSTANTE)}")
        if self.escala_x == 0 or self.escala_y == 0:
            errores.append("ESCALA_X y ESCALA_Y no pueden ser cero")
        if not self.calibracion:
            errores.append("falta CALIBRACION, el archivo con la curva del suelo")
        if errores:
            raise ConfiguracionError(f"{self.origen}: " + "; ".join(errores))


def analizar(texto: str, origen: str = "<memoria>") -> ConfiguracionHumedad:
    """Interpreta el contenido de un archivo ``.HUM``."""
    valores = dict(PREDETERMINADOS)
    desconocidas = []
    for numero, linea in enumerate(texto.splitlines(), 1):
        limpia = linea.split("!", 1)[0].strip()
        if not limpia:
            continue
        if "=" not in limpia:
            raise ConfiguracionError(f"{origen}:{numero}: se esperaba 'CLAVE = valor' y se "
                                     f"leyó {linea.strip()!r}")
        clave, valor = (parte.strip() for parte in limpia.split("=", 1))
        clave = clave.upper()
        if clave not in PREDETERMINADOS:
            desconocidas.append(clave)
        valores[clave] = valor

    def numero_real(clave: str) -> float:
        try:
            return float(valores[clave].replace(",", "."))
        except ValueError:
            raise ConfiguracionError(f"{origen}: {clave} debe ser un número y vale "
                                     f"{valores[clave]!r}") from None

    def afirmativo(clave: str) -> bool:
        return valores[clave].strip().lower() in ("1", "si", "sí", "true")

    def numero_o_nada(clave: str) -> float | None:
        return numero_real(clave) if valores[clave].strip() else None

    try:
        canal = numero_de_canal(valores["CANAL"] if not valores["CANAL"].lstrip("-").isdigit()
                                else int(valores["CANAL"]))
    except ValueError as error:
        raise ConfiguracionError(f"{origen}: {error}") from None

    configuracion = ConfiguracionHumedad(
        patron_imagenes=valores["IMAGENES"],
        canal=canal,
        sigma=numero_real("SIGMA"),
        referencia_seca=valores["REFERENCIA_SECA"],
        desplazamiento_seco=numero_real("DESPLAZAMIENTO_SECO"),
        desplazamiento_saturado=numero_real("DESPLAZAMIENTO_SATURADO"),
        banda_seca=numero_o_nada("BANDA_SECA"),
        banda_saturada=numero_o_nada("BANDA_SATURADA"),
        calibracion=valores["CALIBRACION"],
        umbral_saturacion=numero_real("UMBRAL_SATURACION"),
        incremental=afirmativo("INCREMENTAL"),
        primer_instante=valores["PRIMER_INSTANTE"].strip().lower(),
        redondeo_legado=afirmativo("REDONDEO_LEGADO"),
        escala_x=numero_real("ESCALA_X"),
        origen_x=numero_real("ORIGEN_X"),
        escala_y=numero_real("ESCALA_Y"),
        origen_y=numero_real("ORIGEN_Y"),
        inclinacion_xy=numero_real("INCLINACION_XY"),
        inclinacion_yx=numero_real("INCLINACION_YX"),
        perspectiva_x=numero_real("PERSPECTIVA_X"),
        perspectiva_y=numero_real("PERSPECTIVA_Y"),
        origen=origen,
        desconocidas=tuple(desconocidas),
    )
    configuracion.validar()
    return configuracion


def leer(ruta: Path) -> ConfiguracionHumedad:
    """Lee el archivo ``.HUM`` de un caso."""
    ruta = Path(ruta)
    if not ruta.exists():
        raise ConfiguracionError(f"{ruta}: no existe el archivo de configuración de humedad")
    return analizar(ruta.read_text(encoding="latin-1"), str(ruta))
