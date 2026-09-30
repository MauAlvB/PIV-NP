"""Paso de los ``.PAR`` de cualquier versión al formato único.

El bloque 3 del ``.PAR`` ha cambiado de orden varias veces, y hay dialectos con los mismos
ocho valores significando cosas distintas. En vez de arrastrar esa ambigüedad en el lector,
los archivos se convierten una vez al formato actual::

    <título del caso>
    BLOQUE 2: N_cel N_nod N_part_celda N_fil Ancho Alto
              2006  2100  3            34    0.212115 0.212115
    BLOQUE 3: del_t total_steps impresion moister version pivlab contour rec track
              0.8   149          1         0       1       1      0       0   0
    BLOQUE 4: s_density porosity
              2650.0    0.4

Los valores se reordenan tal cual, sin volver a formatearlos, así que la conversión no
cambia ni un decimal. Las dos únicas excepciones se anotan siempre:

* ``IPIVLAB`` se escribe según las columnas que de verdad tienen los archivos PIVlab del
  caso. Los ``.PAR`` anteriores a 2024 no traen ese campo y el valor que se asumía (4
  columnas) descoloca la lectura de los archivos de 5.
* ``MOISTER=2`` en un ``.PAR`` anterior significaba "leer los archivos de humedad"; en el
  formato actual significa "calcularla desde las imágenes", así que se escribe 1.

El archivo original se conserva al lado, con el nombre ``<caso>.PAR.orig``.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path

from .config import (
    ConfigError,
    RawPar,
    _ListDirectedReader,
    config_from_blocks,
    find_file,
    names_in_header,
)
from .pivlab_io import VELOCITY_PATTERN, frame_interval_in_header

log = logging.getLogger("pivnp")

#: Sufijo del archivo original que se conserva al convertir.
BACKUP_SUFFIX = ".orig"

#: Nombres de los campos de cada bloque, en el orden del formato único.
GEOMETRY_NAMES = ("N_cel", "N_nod", "N_part_celda", "N_fil", "Ancho", "Alto")
ANALYSIS_NAMES = ("del_t", "total_steps", "impresion", "moister", "version", "pivlab",
                  "contour", "rec", "track")
ANALYSIS_FIELDS = ("DT", "TOTAL_STEPS", "IMPPAS", "MOISTER", "IVERSION", "IPIVLAB",
                   "ICONTOUR", "IREC", "ITR")
SOIL_NAMES = ("s_density", "porosity")

_SANGRIA = " " * len("BLOQUE 3: ")


@dataclass
class Conversion:
    """Qué se ha hecho con un ``.PAR``."""

    path: Path
    changed: bool = False
    backup: Path | None = None
    notes: list[str] = field(default_factory=list)
    error: str | None = None

    def __str__(self) -> str:
        if self.error:
            return f"{self.path.name}: NO se ha podido convertir - {self.error}"
        estado = "convertido" if self.changed else "ya estaba en el formato único"
        return f"{self.path.name}: {estado}" + "".join(f"\n    - {n}" for n in self.notes)


def _analysis_block_of_any_version(reader: _ListDirectedReader,
                                   header: str) -> tuple[dict[str, str], bool]:
    """Bloque 3 de cualquier versión: devuelve los valores por nombre y si ya era el actual.

    Primero se mira la línea de comentario, que es donde cada archivo nombra sus campos y la
    única forma fiable de saber el orden: hay dos dialectos con los mismos ocho valores
    colocados de otra manera. Si no se entiende, se recurre al número de valores.
    """
    values = reader.line_values()
    campos = names_in_header(header)
    if campos is not None and len(campos) == len(values):
        log.info("%s: bloque 3 leído por su cabecera (%s)", reader.source, " ".join(campos))
        return dict(zip(campos, values, strict=True)), "ICONTOUR" in campos

    orden: list[str]
    if len(values) >= 9:
        orden = ["DT", "TOTAL_STEPS", "IMPPAS", "MOISTER", "IVERSION", "IPIVLAB", "ICONTOUR",
                 "IREC", "ITR"]
        return dict(zip(orden, values[:9], strict=True)), True
    if len(values) == 8:
        log.info("%s: .PAR con densidad y porosidad en el bloque 3", reader.source)
        orden = ["DT", "TOTAL_STEPS", "IMPPAS", "MOISTER", "S_DENSITY", "POROSITY",
                 "IVERSION", "ITR"]
        return dict(zip(orden, values, strict=True)), False
    if len(values) in (3, 4):
        log.info("%s: .PAR en formato antiguo (%d valores en el bloque 3); se asumen "
                 "IVERSION=1, IPIVLAB=1, ICONTOUR=0, IREC=0 e ITR=0", reader.source, len(values))
        orden = ["DT", "TOTAL_STEPS", "IMPPAS", "MOISTER"]
        return dict(zip(orden, values, strict=False)), False
    raise ConfigError(
        f"{reader.source}: el bloque 3 tiene {len(values)} valores y su línea de comentario "
        "no dice qué es cada uno; no se reconoce como ninguna versión conocida del .PAR")


def read_any_par_blocks(text: str, source: str = "<PAR>") -> RawPar:
    """Lee un ``.PAR`` de cualquier versión. Solo lo usa la conversión."""
    reader = _ListDirectedReader(text.splitlines(), source)
    title = reader.text()
    reader.comment()
    geometry = reader.values(6)
    analysis, current = _analysis_block_of_any_version(reader, reader.comment())

    density, porosity = analysis.get("S_DENSITY"), analysis.get("POROSITY")
    if density is None:
        if reader.at_end():  # los .PAR antiguos no traen el bloque 4
            density = porosity = "0"
        else:
            reader.comment()
            density, porosity = reader.values(2)
    return RawPar(title, geometry, analysis, density, porosity, current)


def bloque(titulo: str, nombres: tuple[str, ...], valores: list[str]) -> list[str]:
    """Dos líneas alineadas: los nombres de los campos y sus valores debajo."""
    anchos = [max(len(n), len(v)) for n, v in zip(nombres, valores, strict=True)]
    cabecera = " ".join(n.ljust(a) for n, a in zip(nombres, anchos, strict=True))
    fila = " ".join(v.ljust(a) for v, a in zip(valores, anchos, strict=True))
    return [f"{titulo} {cabecera}".rstrip(), f"{_SANGRIA}{fila}".rstrip()]


def to_canonical(crudo: RawPar, moister: int, pivlab_format: int) -> str:
    """Texto del ``.PAR`` en el formato único, a partir de los valores leídos."""
    analisis = dict(crudo.analysis)
    analisis["MOISTER"] = str(moister)
    analisis["IPIVLAB"] = str(pivlab_format)
    valores = [analisis.get(campo, crudo.value(campo)) for campo in ANALYSIS_FIELDS]
    lineas = [crudo.title]
    lineas += bloque("BLOQUE 2:", GEOMETRY_NAMES, list(crudo.geometry))
    lineas += bloque("BLOQUE 3:", ANALYSIS_NAMES, valores)
    lineas += bloque("BLOQUE 4:", SOIL_NAMES, [crudo.density, crudo.porosity])
    return "\n".join(lineas) + "\n"


def columns_in_pivlab_files(case_dir: Path) -> int | None:
    """Columnas que traen los archivos PIVlab del caso, o ``None`` si no se encuentran."""
    try:
        ruta = find_file(case_dir, VELOCITY_PATTERN.format(step=1))
    except (FileNotFoundError, OSError):
        return None
    for linea in ruta.read_text(encoding="latin-1").splitlines()[3:]:
        if linea.strip():
            return len([t for t in linea.replace(",", " ").split() if t])
    return None


def convert_file(path: Path) -> Conversion:
    """Convierte un ``.PAR`` al formato único, guardando el original al lado."""
    path = Path(path)
    resultado = Conversion(path)
    texto = path.read_text(encoding="latin-1")
    try:
        crudo = read_any_par_blocks(texto, str(path))
        config = config_from_blocks(crudo, str(path))
    except (ConfigError, ValueError) as error:
        resultado.error = str(error)
        return resultado

    moister = int(config.moisture) + int(config.moisture_from_images)
    declarado = crudo.analysis.get("MOISTER")
    if declarado is not None and declarado.strip() != str(moister):
        resultado.notes.append(
            f"MOISTER={declarado.strip()} de una versión anterior se escribe como {moister}: "
            "allí significaba leer los archivos de humedad")

    pivlab = config.pivlab_format
    columnas = columns_in_pivlab_files(path.parent)
    if columnas is not None:
        medido = 1 if columnas <= 4 else 2
        if medido != pivlab:
            resultado.notes.append(
                f"IPIVLAB pasa de {pivlab} a {medido}: los archivos PIVlab del caso traen "
                f"{columnas} columnas")
        pivlab = medido
    elif "IPIVLAB" not in crudo.analysis:
        resultado.notes.append("IPIVLAB se deja en 1: no se han encontrado los archivos "
                               "PIVlab para comprobar cuántas columnas tienen")

    _avisar_del_intervalo(path.parent, config, resultado)

    nuevo = to_canonical(crudo, moister, pivlab)
    if nuevo == texto:
        return resultado

    copia = path.with_name(path.name + BACKUP_SUFFIX)
    if not copia.exists():
        copia.write_text(texto, encoding="latin-1")
        resultado.backup = copia
    else:
        resultado.notes.append(f"{copia.name} ya existía: se conserva el de la primera vez")
    path.write_text(nuevo, encoding="latin-1")
    resultado.changed = True
    return resultado


def _avisar_del_intervalo(case_dir: Path, config, resultado: Conversion) -> None:
    """El DT que no cuadra con PIVlab se avisa, pero no se toca: cambiaría los resultados."""
    try:
        intervalo = frame_interval_in_header(find_file(case_dir,
                                                       VELOCITY_PATTERN.format(step=1)))
    except (FileNotFoundError, OSError, ValueError):
        return
    if intervalo and abs(intervalo - config.dt) > 1e-3 * max(intervalo, config.dt):
        resultado.notes.append(
            f"DT={config.dt:g} no coincide con el intervalo {intervalo:g} s con el que se "
            "exportó desde PIVlab; se deja como está porque cambiarlo cambiaría los "
            "resultados")


def convert_tree(root: Path) -> list[Conversion]:
    """Convierte todos los ``.PAR`` que haya bajo un directorio."""
    root = Path(root)
    rutas = sorted(p for p in root.rglob("*.PAR") if p.is_file() and p.suffix == ".PAR")
    if not rutas and root.is_file():
        rutas = [root]
    return [convert_file(ruta) for ruta in rutas]
