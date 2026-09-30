"""Lectura y validación de los datos del caso (``PIV-NP.TXT`` y ``<caso>.PAR``).

Hay un solo formato de ``.PAR``. La lectura es *list-directed* de Fortran: los valores se
separan por espacios, tabuladores o comas, y el programa original lo lee igual::

    <título del caso>
    BLOQUE 2: N_cel N_nod N_part_celda N_fil Ancho    Alto
              2006  2100  3            34    0.212115 0.212115
    BLOQUE 3: del_t total_steps impresion moister version pivlab contour rec track
              0.8   149         1         0       1       1      0       0   0
    BLOQUE 4: s_density porosity
              2650.0    0.4

Los nombres van encima de sus valores, así que el archivo se explica solo. El bloque 3 es,
por orden: el tiempo entre imágenes, el número de archivos PIVlab, cada cuántos pasos se
imprimen resultados, de dónde sale la humedad, la versión de la malla, el formato de los
archivos PIVlab, la corrección de contorno, si se continúa un análisis y las partículas de
seguimiento.

* MOISTER: 0 sin humedad, 1 leerla de los ``Moist_<n>.TXT``, 2 calcularla desde las imágenes
  del ensayo con la configuración de ``<caso>.HUM``.
* ITR (partículas de seguimiento PTV) debe ser 0: esa opción ya no está soportada.

El bloque 3 cambió de orden varias veces entre versiones del programa, y llegó a haber
dialectos con los mismos valores colocados de otra manera. Aquí no se leen: un ``.PAR`` de
una versión anterior se pasa una sola vez al formato de arriba con ``pivnp <directorio>
--convert-par`` (ver :mod:`pivnp.par_migrate`), que guarda el original al lado.
"""

from __future__ import annotations

import logging
import math
import re
import unicodedata
from dataclasses import dataclass
from pathlib import Path

log = logging.getLogger("pivnp")

CASE_INDEX_FILE = "PIV-NP.TXT"

_TOKEN_SEPARATORS = re.compile(r"[,\s]+")


class ConfigError(ValueError):
    """Datos de entrada ausentes o incoherentes."""


#: Máximo de partículas por lado de celda: con más, el original las colocaba todas en el
#: centro de la celda.
MAX_PARTICLES_PER_SIDE = 6


@dataclass(frozen=True)
class CaseConfig:
    """Parámetros del análisis. Entre paréntesis, el nombre de la variable en el Fortran."""

    title: str
    n_cells: int  # NC: celdas de la malla PIVlab
    n_nodes: int  # NN: nodos (puntos) PIVlab
    particles_per_side: int  # NPC: partículas por lado de celda (NPC² por celda)
    n_rows: int  # NFIL: filas de celdas
    cell_width: float  # AXC [m]
    cell_height: float  # AYC [m]
    dt: float  # DT: tiempo entre imágenes
    total_steps: int  # TOTAL_STEPS: número de archivos PIVlab
    print_every: int  # IMPPAS: pasos entre resultados impresos
    moisture: bool  # MOISTER >= 1: el análisis lleva humedad y saturación
    #: MOISTER = 2: la humedad se calcula desde las imágenes del ensayo, con la
    #: configuración de ``<caso>.HUM``, en vez de leerse de los ``Moist_<n>.TXT``.
    moisture_from_images: bool
    mesh_version: int  # IVERSION: 1 = malla PIV-NP = malla PIVlab; 2 = malla desplazada
    pivlab_format: int  # IPIVLAB: 1 = 4 columnas (PIVlab antiguo); otro = 5 columnas
    contour: int  # ICONTOUR: corrección de velocidades en el contorno (0 = no)
    restart: bool  # IREC: continuar desde <caso>.REC
    soil_density: float  # S_DENSITY [kg/m3] (no se usa en el cálculo)
    porosity: float  # POROSITY (no se usa en el cálculo)

    @property
    def n_cols(self) -> int:
        """NCH: celdas por fila."""
        return self.n_cells // self.n_rows

    @property
    def n_particles(self) -> int:
        """NP: partículas generadas en la malla."""
        cells = self.n_cells if self.mesh_version == 1 else (self.n_cols + 1) * (self.n_rows + 1)
        return cells * self.particles_per_side**2

    def validate(self) -> None:
        """Rechaza combinaciones con las que el original produce resultados sin sentido."""
        errors = []
        if self.n_rows < 1 or self.n_cells < 1:
            errors.append("NC y NFIL deben ser positivos")
        elif self.n_cells % self.n_rows:
            errors.append(f"NC={self.n_cells} no es múltiplo de NFIL={self.n_rows}")
        elif self.n_nodes != (self.n_cols + 1) * (self.n_rows + 1):
            errors.append(
                f"NN={self.n_nodes} no coincide con (NC/NFIL+1)*(NFIL+1)="
                f"{(self.n_cols + 1) * (self.n_rows + 1)}"
            )
        if not 1 <= self.particles_per_side <= MAX_PARTICLES_PER_SIDE:
            errors.append(f"NPC={self.particles_per_side} debe estar entre 1 y "
                          f"{MAX_PARTICLES_PER_SIDE}")
        if self.cell_width <= 0 or self.cell_height <= 0:
            errors.append("AXC y AYC deben ser positivos")
        if self.mesh_version not in (1, 2):
            errors.append(f"IVERSION={self.mesh_version} debe ser 1 o 2")
        if self.print_every < 1:
            errors.append("IMPPAS debe ser >= 1")
        if self.contour not in (0, 1, 2, 3):
            errors.append(f"ICONTOUR={self.contour} debe estar entre 0 y 3")
        if errors:
            raise ConfigError("; ".join(errors))


class _ListDirectedReader:
    """Imita las sentencias ``READ(u, 2000)`` (A80) y ``READ(u, *)`` de Fortran."""

    def __init__(self, lines: list[str], source: str) -> None:
        self._lines = lines
        self._pos = 0
        self.source = source

    def _next_line(self) -> str:
        if self._pos >= len(self._lines):
            raise ConfigError(f"{self.source}: fin de archivo inesperado (línea {self._pos + 1})")
        line = self._lines[self._pos]
        self._pos += 1
        return line

    def text(self) -> str:
        """Una línea de texto, truncada a 80 columnas como el formato A80 del original."""
        return self._next_line()[:80].rstrip()

    def comment(self) -> str:
        """Una línea de comentario, entera: aquí el original no lee nada, solo salta."""
        return self._next_line().rstrip()

    def values(self, count: int) -> list[str]:
        """Lee ``count`` valores; el resto de la última línea leída se descarta."""
        tokens: list[str] = []
        while len(tokens) < count:
            tokens.extend(t for t in _TOKEN_SEPARATORS.split(self._next_line().strip()) if t)
        return tokens[:count]

    def line_values(self) -> list[str]:
        """Todos los valores de la línea siguiente."""
        return [t for t in _TOKEN_SEPARATORS.split(self._next_line().strip()) if t]

    def at_end(self) -> bool:
        return self._pos >= len(self._lines) or not any(
            line.strip() for line in self._lines[self._pos:]
        )


def _to_int(token: str, name: str) -> int:
    try:
        return int(token)
    except ValueError:
        raise ConfigError(f"{name}: se esperaba un entero y se leyó {token!r}") from None


def _to_float(token: str, name: str) -> float:
    try:
        return float(token.replace("D", "E").replace("d", "e"))
    except ValueError:
        raise ConfigError(f"{name}: se esperaba un número y se leyó {token!r}") from None


#: Nombre de campo que corresponde a cada palabra de la línea de comentario del bloque 3.
#: Se busca como subcadena, sin acentos y en minúsculas, en el orden de esta lista.
_BLOCK3_NAMES: tuple[tuple[str, str], ...] = (
    ("del_t", "DT"), ("delt", "DT"), ("dt", "DT"),
    ("total", "TOTAL_STEPS"), ("steps", "TOTAL_STEPS"), ("pasos", "TOTAL_STEPS"),
    ("impres", "IMPPAS"), ("salto", "IMPPAS"), ("print", "IMPPAS"),
    ("moist", "MOISTER"), ("humedad", "MOISTER"),
    ("densi", "S_DENSITY"),
    ("porosi", "POROSITY"),
    ("pivlab", "IPIVLAB"),
    ("pivnp", "IVERSION"), ("version", "IVERSION"),
    ("contour", "ICONTOUR"), ("contorno", "ICONTOUR"),
    ("rec", "IREC"),
    ("ptv", "ITR"), ("ptr", "ITR"), ("track", "ITR"), ("segui", "ITR"),
)

#: Valor por defecto de lo que un ``.PAR`` puede no traer.
_BLOCK3_DEFAULTS = {"MOISTER": "0", "IVERSION": "1", "IPIVLAB": "1", "ICONTOUR": "0",
                    "IREC": "0", "ITR": "0"}

_PARENTHESES = re.compile(r"\([^)]*\)")


def _strip_accents(text: str) -> str:
    return "".join(c for c in unicodedata.normalize("NFD", text)
                   if unicodedata.category(c) != "Mn")


def names_in_header(header: str) -> list[str] | None:
    """Campos que nombra la línea de comentario del bloque 3, o ``None`` si no se entiende.

    Cada ``.PAR`` documenta su propio orden en esa línea (``del_t total_steps
    salto_impresion v.pivnp v.pivlab moister REC PTR``), que es la única forma fiable de
    saberlo: hay dialectos con los mismos ocho valores en distinto orden.
    """
    texto = _PARENTHESES.sub(" ", header)
    _, _, despues = texto.partition(":")  # quita el "BLOQUE 3:" de delante
    palabras = _strip_accents(despues or texto).lower().split()
    if not palabras:
        return None

    campos: list[str] = []
    for palabra in palabras:
        for clave, campo in _BLOCK3_NAMES:
            if clave in palabra:
                if campo in campos:  # un campo repetido delata que no es una lista de nombres
                    return None
                campos.append(campo)
                break
        else:
            return None  # una palabra que no se reconoce invalida toda la línea
    return campos


#: Orden de los valores del bloque 3 en el formato único.
ANALYSIS_ORDER = ("DT", "TOTAL_STEPS", "IMPPAS", "MOISTER", "IVERSION", "IPIVLAB", "ICONTOUR",
                  "IREC", "ITR")

_COMO_CONVERTIR = ("Si viene de una versión anterior del programa, pásalo al formato único "
                   "con:  pivnp <directorio> --convert-par")


def _analysis_block(reader: _ListDirectedReader, header: str) -> dict[str, str]:
    """Bloque 3 del formato único: nueve valores en el orden de :data:`ANALYSIS_ORDER`."""
    values = reader.line_values()
    if len(values) != len(ANALYSIS_ORDER):
        raise ConfigError(
            f"{reader.source}: el bloque 3 tiene {len(values)} valores y el formato único "
            f"tiene {len(ANALYSIS_ORDER)} ({' '.join(ANALYSIS_ORDER)}). {_COMO_CONVERTIR}")
    campos = names_in_header(header)
    if campos is not None and tuple(campos) != ANALYSIS_ORDER:
        raise ConfigError(
            f"{reader.source}: la línea de comentario del bloque 3 nombra los campos en otro "
            f"orden ({' '.join(campos)}). {_COMO_CONVERTIR}")
    return dict(zip(ANALYSIS_ORDER, values, strict=True))


@dataclass(frozen=True)
class RawPar:
    """Los valores de un ``.PAR``, tal cual vienen, ya identificados por su nombre.

    Guarda las cadenas y no los números, de modo que convertir un archivo al formato único
    sea reordenar lo que hay, sin volver a formatear ni redondear nada.
    """

    title: str
    geometry: list[str]  # NC NN NPC NFIL AXC AYC
    analysis: dict[str, str]  # lo que traiga el bloque 3, por nombre
    density: str
    porosity: str
    #: El ``.PAR`` está en el formato único actual. Solo el conversor lee los anteriores.
    current_dialect: bool = True

    def value(self, name: str) -> str:
        """Valor del bloque 3, o el que se asume cuando el archivo no lo trae."""
        if name in self.analysis:
            return self.analysis[name]
        return _BLOCK3_DEFAULTS[name]


def read_par_blocks(text: str, source: str = "<PAR>") -> RawPar:
    """Lee un ``.PAR`` en el formato único y devuelve sus valores por nombre.

    Los ``.PAR`` de versiones anteriores no se leen aquí: se convierten una vez con
    ``pivnp <directorio> --convert-par`` y a partir de ahí hay un solo formato de entrada.
    """
    reader = _ListDirectedReader(text.splitlines(), source)
    title = reader.text()
    reader.comment()
    geometry = reader.values(6)
    analysis = _analysis_block(reader, reader.comment())
    if reader.at_end():
        raise ConfigError(f"{source}: falta el bloque 4, con la densidad del suelo y la "
                          f"porosidad. {_COMO_CONVERTIR}")
    reader.comment()
    density, porosity = reader.values(2)
    return RawPar(title, geometry, analysis, density, porosity)


def parse_par(text: str, source: str = "<PAR>") -> CaseConfig:
    """Interpreta el contenido de un archivo ``.PAR`` (formato actual o anterior)."""
    return config_from_blocks(read_par_blocks(text, source), source)


def config_from_blocks(crudo: RawPar, source: str = "<PAR>") -> CaseConfig:
    """Valida los valores leídos y los convierte en la configuración del caso."""
    nc, nn, npc, nfil, axc, ayc = crudo.geometry
    dt, steps, imppas = (crudo.value("DT"), crudo.value("TOTAL_STEPS"), crudo.value("IMPPAS"))
    moister, iversion = crudo.value("MOISTER"), crudo.value("IVERSION")
    ipivlab, icontour = crudo.value("IPIVLAB"), crudo.value("ICONTOUR")
    irec, itr = crudo.value("IREC"), crudo.value("ITR")
    title, density, porosity = crudo.title, crudo.density, crudo.porosity
    dialecto_actual = crudo.current_dialect

    if _to_int(itr, "ITR") != 0:
        raise ConfigError(
            "ITR: las partículas de seguimiento PTV ya no están soportadas. Usa ITR=0 y "
            "quita el bloque 5 del .PAR"
        )

    irec_value = _to_int(irec, "IREC")
    if irec_value not in (0, 1):
        raise ConfigError(f"IREC={irec_value} debe ser 0 o 1")

    moister_value = _to_int(moister, "MOISTER")
    if not dialecto_actual and moister_value > 1:
        # En las versiones anteriores cualquier valor distinto de 0 activaba la lectura de
        # los archivos de humedad; el 2 de "calcular desde las imágenes" es nuevo.
        log.info("%s: MOISTER=%d en un .PAR de una versión anterior; se entiende como 1, "
                 "leer los Moist_<n>.TXT", source, moister_value)
        moister_value = 1
    if moister_value not in (0, 1, 2):
        raise ConfigError(f"MOISTER={moister_value} debe ser 0 (sin humedad), 1 (leerla de "
                          "los Moist_<n>.TXT) o 2 (calcularla desde las imágenes)")

    config = CaseConfig(
        title=title,
        n_cells=_to_int(nc, "NC"),
        n_nodes=_to_int(nn, "NN"),
        particles_per_side=_to_int(npc, "NPC"),
        n_rows=_to_int(nfil, "NFIL"),
        cell_width=_to_float(axc, "AXC"),
        cell_height=_to_float(ayc, "AYC"),
        dt=_to_float(dt, "DT"),
        # TOTAL_STEPS es REAL en el original: el bucle DO trunca su valor.
        total_steps=max(0, math.trunc(_to_float(steps, "TOTAL_STEPS"))),
        print_every=_to_int(imppas, "IMPPAS"),
        moisture=moister_value >= 1,
        moisture_from_images=moister_value == 2,
        mesh_version=_to_int(iversion, "IVERSION"),
        pivlab_format=_to_int(ipivlab, "IPIVLAB"),
        contour=_to_int(icontour, "ICONTOUR"),
        restart=irec_value == 1,
        soil_density=_to_float(density, "S_DENSITY"),
        porosity=_to_float(porosity, "POROSITY"),
    )
    config.validate()
    return config


def read_case_name(case_dir: Path) -> str:
    """Lee el nombre del caso desde ``PIV-NP.TXT`` (primer valor de la primera línea)."""
    path = find_file(case_dir, CASE_INDEX_FILE)
    first_line = path.read_text(encoding="latin-1").strip().splitlines()[0].strip()
    if first_line[:1] in ("'", '"'):
        return first_line[1:].split(first_line[0], 1)[0]
    return _TOKEN_SEPARATORS.split(first_line, 1)[0]


def load_case(case_dir: Path, case_name: str | None = None) -> tuple[str, CaseConfig]:
    """Devuelve el nombre del caso y su configuración validada."""
    case_dir = Path(case_dir)
    name = case_name or read_case_name(case_dir)
    par_path = find_file(case_dir, f"{name}.PAR")
    return name, parse_par(par_path.read_text(encoding="latin-1"), str(par_path))


def find_file(directory: Path, name: str) -> Path:
    """Busca ``name`` en ``directory`` sin distinguir mayúsculas (como Windows)."""
    candidate = Path(directory) / name
    if candidate.exists():
        return candidate
    lowered = name.lower()
    for entry in Path(directory).iterdir():
        if entry.name.lower() == lowered:
            return entry
    raise FileNotFoundError(candidate)
