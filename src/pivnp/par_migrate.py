"""Moving the ``.PAR`` files of any version to the single format.

Block 3 of the ``.PAR`` has changed order several times, and there are dialects carrying the
same eight values with different meanings. Instead of dragging that ambiguity through the
reader, the files are converted once to the current format::

    <case title>
    BLOCK 2: n_cells n_nodes n_part_cell n_rows width    height
             2006    2100    3           34     0.212115 0.212115
    BLOCK 3: dt  total_steps print_every moisture mesh_version pivlab contour restart track
             0.8 149         1           0        1            1      0       0       0
    BLOCK 4: s_density porosity
             2650.0    0.4
    ! ... legend with the meaning of every value, which is not read ...

The values are reordered as they are, without reformatting them, so the conversion does not
change a single decimal. The only two exceptions are always reported:

* ``IPIVLAB`` is written from the number of columns the PIVlab files of the case actually
  have. ``.PAR`` files older than 2024 do not carry that field, and the value that used to be
  assumed (4 columns) shifts every value when reading 5-column files.
* ``MOISTER=2`` in an older ``.PAR`` meant "read the moisture files"; in the current format it
  means "compute it from the images", so 1 is written instead.

The original file is kept next to it, named ``<case>.PAR.orig``.
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

#: Suffix of the original file kept when converting.
BACKUP_SUFFIX = ".orig"

#: Names of the fields of every block, in the order of the single format.
GEOMETRY_NAMES = ("n_cells", "n_nodes", "n_part_cell", "n_rows", "width", "height")
ANALYSIS_NAMES = ("dt", "total_steps", "print_every", "moisture", "mesh_version", "pivlab",
                  "contour", "restart", "track")
ANALYSIS_FIELDS = ("DT", "TOTAL_STEPS", "IMPPAS", "MOISTER", "IVERSION", "IPIVLAB",
                   "ICONTOUR", "IREC", "ITR")
SOIL_NAMES = ("s_density", "porosity")

_INDENT = " " * len("BLOCK 3: ")

#: Legend written at the end of the file. It is not read: it is there for whoever opens it.
LEGEND = """\
!------------------------------------------------------------------------------------------
! What every value means. Nothing below this line is read.
!
! BLOCK 2 - geometry of the PIVlab grid
!   n_cells       cells of the grid
!   n_nodes       nodes (points) of the grid; n_nodes = (n_cells/n_rows + 1) * (n_rows + 1)
!   n_part_cell   rows and columns of particles per cell, from 1 to 6 (its square per cell)
!   n_rows        rows of cells
!   width, height size of one cell, in metres
!
! BLOCK 3 - which analysis is run
!   dt            time between images, in seconds. It has to be the same interval the files
!                 were exported from PIVlab with; otherwise displacements come out at a
!                 different scale. The program checks the file header and warns.
!   total_steps   how many PIVlab files are processed
!   print_every   results are written at step 1 and then every 'print_every' steps
!   moisture      where the moisture comes from
!                   0 = no moisture is computed
!                   1 = it is read from the Moist_<n>.TXT files
!                   2 = it is computed from the test images, with the settings of the
!                       <case>.HUM file
!   mesh_version  where the velocities given by PIVlab sit
!                   1 = at the nodes of the grid
!                   2 = at the centre of each element (grid shifted half a cell)
!   pivlab        version of PIVlab the "Datos" files were exported with, whose format
!                 changes from one version to the next
!                   1 = 4-column export (x, y, u, v)
!                   2 = 5-column export (adds the vector type)
!                 The program looks at the file and warns if it does not match what is here.
!   contour       correction of the velocity at the boundary of the material
!                   0 = none
!                   1 = average of the neighbouring nodes that do have data
!                   2 = average of the particles surrounding the point
!                   3 = extrapolation from the interior towards the exterior
!   restart       0 = new analysis
!                 1 = continue from <case>.REC, for instance a second stage of the test
!   track         must be 0; PTV tracking particles are no longer supported
!
! BLOCK 4 - the soil
!   s_density     density of the soil, in kg/m3
!   porosity      initial porosity
!------------------------------------------------------------------------------------------
"""


@dataclass
class Conversion:
    """What was done with one ``.PAR``."""

    path: Path
    changed: bool = False
    backup: Path | None = None
    notes: list[str] = field(default_factory=list)
    error: str | None = None

    def __str__(self) -> str:
        if self.error:
            return f"{self.path.name}: could NOT be converted - {self.error}"
        state = "converted" if self.changed else "already in the single format"
        return f"{self.path.name}: {state}" + "".join(f"\n    - {n}" for n in self.notes)


def _analysis_block_of_any_version(reader: _ListDirectedReader,
                                   header: str) -> tuple[dict[str, str], bool]:
    """Block 3 of any version: returns the values by name and whether it was already current.

    The comment line is looked at first, because that is where every file names its fields and
    the only reliable way to know the order: there are two dialects carrying the same eight
    values arranged differently. If it is not understood, the number of values is used.
    """
    values = reader.line_values()
    fields = names_in_header(header)
    if fields is not None and len(fields) == len(values):
        log.info("%s: block 3 read from its header (%s)", reader.source, " ".join(fields))
        return dict(zip(fields, values, strict=True)), "ICONTOUR" in fields

    order: list[str]
    if len(values) >= 9:
        order = ["DT", "TOTAL_STEPS", "IMPPAS", "MOISTER", "IVERSION", "IPIVLAB", "ICONTOUR",
                 "IREC", "ITR"]
        return dict(zip(order, values[:9], strict=True)), True
    if len(values) == 8:
        log.info("%s: .PAR with density and porosity inside block 3", reader.source)
        order = ["DT", "TOTAL_STEPS", "IMPPAS", "MOISTER", "S_DENSITY", "POROSITY",
                 "IVERSION", "ITR"]
        return dict(zip(order, values, strict=True)), False
    if len(values) in (3, 4):
        log.info("%s: .PAR in the old format (%d values in block 3); IVERSION=1, IPIVLAB=1, "
                 "ICONTOUR=0, IREC=0 and ITR=0 are assumed", reader.source, len(values))
        order = ["DT", "TOTAL_STEPS", "IMPPAS", "MOISTER"]
        return dict(zip(order, values, strict=False)), False
    raise ConfigError(
        f"{reader.source}: block 3 carries {len(values)} values and its comment line does not "
        "say what each one is; it is not recognized as any known version of the .PAR")


def read_any_par_blocks(text: str, source: str = "<PAR>") -> RawPar:
    """Read a ``.PAR`` of any version. Only the conversion uses it."""
    reader = _ListDirectedReader(text.splitlines(), source)
    title = reader.text()
    reader.comment()
    geometry = reader.values(6)
    analysis, current = _analysis_block_of_any_version(reader, reader.comment())

    density, porosity = analysis.get("S_DENSITY"), analysis.get("POROSITY")
    if density is None:
        if reader.at_end():  # the old .PAR files do not carry block 4
            density = porosity = "0"
        else:
            reader.comment()
            density, porosity = reader.values(2)
    return RawPar(title, geometry, analysis, density, porosity, current)


def block(title: str, names: tuple[str, ...], values: list[str]) -> list[str]:
    """Two aligned lines: the names of the fields and their values underneath."""
    widths = [max(len(n), len(v)) for n, v in zip(names, values, strict=True)]
    header = " ".join(n.ljust(w) for n, w in zip(names, widths, strict=True))
    row = " ".join(v.ljust(w) for v, w in zip(values, widths, strict=True))
    return [f"{title} {header}".rstrip(), f"{_INDENT}{row}".rstrip()]


def to_canonical(raw: RawPar, moister: int, pivlab_format: int) -> str:
    """Text of the ``.PAR`` in the single format, from the values that were read."""
    analysis = dict(raw.analysis)
    analysis["MOISTER"] = str(moister)
    analysis["IPIVLAB"] = str(pivlab_format)
    values = [analysis.get(field_name, raw.value(field_name)) for field_name in ANALYSIS_FIELDS]
    lines = [raw.title]
    lines += block("BLOCK 2:", GEOMETRY_NAMES, list(raw.geometry))
    lines += block("BLOCK 3:", ANALYSIS_NAMES, values)
    lines += block("BLOCK 4:", SOIL_NAMES, [raw.density, raw.porosity])
    return "\n".join(lines) + "\n" + LEGEND


def columns_in_pivlab_files(case_dir: Path) -> int | None:
    """Columns the PIVlab files of the case carry, or ``None`` when they are not found."""
    try:
        path = find_file(case_dir, VELOCITY_PATTERN.format(step=1))
    except (FileNotFoundError, OSError):
        return None
    for line in path.read_text(encoding="latin-1").splitlines()[3:]:
        if line.strip():
            return len([t for t in line.replace(",", " ").split() if t])
    return None


def convert_file(path: Path) -> Conversion:
    """Convert a ``.PAR`` to the single format, keeping the original next to it."""
    path = Path(path)
    result = Conversion(path)
    text = path.read_text(encoding="latin-1")
    try:
        raw = read_any_par_blocks(text, str(path))
        config = config_from_blocks(raw, str(path))
    except (ConfigError, ValueError) as error:
        result.error = str(error)
        return result

    moister = int(config.moisture) + int(config.moisture_from_images)
    declared = raw.analysis.get("MOISTER")
    if declared is not None and declared.strip() != str(moister):
        result.notes.append(
            f"MOISTER={declared.strip()} from an earlier version is written as {moister}: "
            "there it meant reading the moisture files")

    pivlab = config.pivlab_format
    columns = columns_in_pivlab_files(path.parent)
    if columns is not None:
        measured = 1 if columns <= 4 else 2
        if measured != pivlab:
            result.notes.append(
                f"IPIVLAB goes from {pivlab} to {measured}: the PIVlab files of the case "
                f"carry {columns} columns")
        pivlab = measured
    elif "IPIVLAB" not in raw.analysis:
        result.notes.append("IPIVLAB is left at 1: the PIVlab files of the case were not "
                            "found, so their number of columns could not be checked")

    _warn_about_the_interval(path.parent, config, result)

    new_text = to_canonical(raw, moister, pivlab)
    if new_text == text:
        return result

    backup = path.with_name(path.name + BACKUP_SUFFIX)
    if not backup.exists():
        backup.write_text(text, encoding="latin-1")
        result.backup = backup
    else:
        result.notes.append(f"{backup.name} already existed: the first one is kept")
    path.write_text(new_text, encoding="latin-1")
    result.changed = True
    return result


def _warn_about_the_interval(case_dir: Path, config, result: Conversion) -> None:
    """A DT that does not match PIVlab is reported, not touched: it would change results."""
    try:
        interval = frame_interval_in_header(find_file(case_dir,
                                                      VELOCITY_PATTERN.format(step=1)))
    except (FileNotFoundError, OSError, ValueError):
        return
    if interval and abs(interval - config.dt) > 1e-3 * max(interval, config.dt):
        result.notes.append(
            f"DT={config.dt:g} does not match the interval of {interval:g} s the files were "
            "exported from PIVlab with; it is left as it is, because changing it would change "
            "the results")


def convert_tree(root: Path) -> list[Conversion]:
    """Convert every ``.PAR`` under a directory."""
    root = Path(root)
    paths = sorted(p for p in root.rglob("*.PAR") if p.is_file() and p.suffix == ".PAR")
    if not paths and root.is_file():
        paths = [root]
    return [convert_file(path) for path in paths]
