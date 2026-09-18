"""Lectura y validación de los datos del caso (``PIV-NP.TXT`` y ``<caso>.PAR``).

Formato del archivo ``.PAR`` (lectura *list-directed* de Fortran: los valores se separan
por espacios, tabuladores o comas y pueden repartirse en varias líneas)::

    Bloque 1  título (una línea)
    Bloque 2  línea de comentario
              NC  NN  NPC  NFIL  AXC  AYC
    Bloque 3  línea de comentario
              DT  TOTAL_STEPS  IMPPAS  MOISTER  IVERSION  IPIVLAB  ICONTOUR  IREC  ITR
    Bloque 4  línea de comentario
              S_DENSITY  POROSITY
    Bloque 5  (solo si ITR != 0) línea de comentario
              DXT  DYT  PTVX1  PTVY1  PTVX2  PTVY2  PTVX3  PTVY3
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from pathlib import Path

CASE_INDEX_FILE = "PIV-NP.TXT"

_TOKEN_SEPARATORS = re.compile(r"[,\s]+")


class ConfigError(ValueError):
    """Datos de entrada ausentes o incoherentes."""


@dataclass(frozen=True)
class TrackingPoints:
    """Bloque 5: partículas adicionales para comparar con seguimiento PTV de laboratorio."""

    scale_x: float  # DXT
    scale_y: float  # DYT
    points: tuple[tuple[float, float], ...]  # (PTVXi, PTVYi), i = 1..3

    def positions(self) -> list[tuple[float, float]]:
        return [(self.scale_x * px, self.scale_y * py) for px, py in self.points]


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
    moisture: bool  # MOISTER: leer archivos Moist_<n>.TXT
    mesh_version: int  # IVERSION: 1 = malla PIV-NP = malla PIVlab; 2 = malla desplazada
    pivlab_format: int  # IPIVLAB: 1 = 4 columnas (PIVlab antiguo); otro = 5 columnas
    contour: int  # ICONTOUR: corrección de velocidades en el contorno (0 = no)
    restart: bool  # IREC: continuar desde <caso>.REC
    tracking: TrackingPoints | None  # ITR: partículas PTV adicionales
    soil_density: float  # S_DENSITY [kg/m3] (no se usa en el cálculo)
    porosity: float  # POROSITY (no se usa en el cálculo)

    @property
    def n_cols(self) -> int:
        """NCH: celdas por fila."""
        return self.n_cells // self.n_rows

    @property
    def n_base_particles(self) -> int:
        """NP0: partículas generadas en la malla (sin las de seguimiento PTV)."""
        cells = self.n_cells if self.mesh_version == 1 else (self.n_cols + 1) * (self.n_rows + 1)
        return cells * self.particles_per_side**2

    @property
    def n_particles(self) -> int:
        """NP: total de partículas, incluidas las de seguimiento PTV."""
        extra = len(self.tracking.points) if self.tracking else 0
        return self.n_base_particles + extra

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
        if self.particles_per_side < 1:
            errors.append("NPC debe ser >= 1")
        if self.cell_width <= 0 or self.cell_height <= 0:
            errors.append("AXC y AYC deben ser positivos")
        if self.mesh_version not in (1, 2):
            errors.append(f"IVERSION={self.mesh_version} debe ser 1 o 2")
        if self.print_every < 1:
            errors.append("IMPPAS debe ser >= 1")
        if errors:
            raise ConfigError("; ".join(errors))


class _ListDirectedReader:
    """Imita las sentencias ``READ(u, 2000)`` (A80) y ``READ(u, *)`` de Fortran."""

    def __init__(self, lines: list[str], source: str) -> None:
        self._lines = lines
        self._pos = 0
        self._source = source

    def _next_line(self) -> str:
        if self._pos >= len(self._lines):
            raise ConfigError(f"{self._source}: fin de archivo inesperado (línea {self._pos + 1})")
        line = self._lines[self._pos]
        self._pos += 1
        return line

    def text(self) -> str:
        return self._next_line()[:80].rstrip()

    def values(self, count: int) -> list[str]:
        """Lee ``count`` valores; el resto de la última línea leída se descarta."""
        tokens: list[str] = []
        while len(tokens) < count:
            tokens.extend(t for t in _TOKEN_SEPARATORS.split(self._next_line().strip()) if t)
        return tokens[:count]


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


def parse_par(text: str, source: str = "<PAR>") -> CaseConfig:
    """Interpreta el contenido de un archivo ``.PAR``."""
    reader = _ListDirectedReader(text.splitlines(), source)

    title = reader.text()
    reader.text()
    nc, nn, npc, nfil, axc, ayc = reader.values(6)
    reader.text()
    dt, steps, imppas, moister, iversion, ipivlab, icontour, irec, itr = reader.values(9)
    reader.text()
    density, porosity = reader.values(2)

    tracking = None
    if _to_int(itr, "ITR") != 0:
        reader.text()
        values = [_to_float(v, "PTV") for v in reader.values(8)]
        tracking = TrackingPoints(
            scale_x=values[0],
            scale_y=values[1],
            points=((values[2], values[3]), (values[4], values[5]), (values[6], values[7])),
        )

    irec_value = _to_int(irec, "IREC")
    if irec_value not in (0, 1):
        raise ConfigError(f"IREC={irec_value} debe ser 0 o 1")

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
        moisture=_to_int(moister, "MOISTER") == 1,
        mesh_version=_to_int(iversion, "IVERSION"),
        pivlab_format=_to_int(ipivlab, "IPIVLAB"),
        contour=_to_int(icontour, "ICONTOUR"),
        restart=irec_value == 1,
        tracking=tracking,
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
