"""Reading and validation of the case input (``PIV-NP.TXT`` and ``<case>.PAR``).

There is a single ``.PAR`` format. Reading is Fortran *list-directed*: values are separated
by spaces, tabs or commas, and the original program reads it the same way::

    <case title>
    BLOCK 2: n_cells n_nodes n_part_cell n_rows width    height
             2006    2100    3           34     0.212115 0.212115
    BLOCK 3: dt  total_steps print_every moisture mesh_version pivlab contour restart track
             0.8 149         1           0        1            1      0       0       0
    BLOCK 4: s_density porosity
             2650.0    0.4

The names sit above their values, so the file explains itself. Block 3 is, in order: the time
between images, the number of PIVlab files, how often results are printed, where the moisture
comes from, the grid version, the format of the PIVlab files, the boundary correction, whether
an analysis is continued, and the tracking particles.

* MOISTURE: 0 no moisture, 1 read it from the ``Moist_<n>.TXT``, 2 compute it from the test
  images with the settings of ``<case>.HUM``.
* TRACK (PTV tracking particles) must be 0: that option is no longer supported.

Block 3 changed order several times between versions of the program, and there were even
dialects with the same values arranged differently. They are not read here: a ``.PAR`` from
an earlier version is moved once to the format above with ``pivnp <directory> --convert-par``
(see :mod:`pivnp.par_migrate`), which keeps the original next to it.
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
    """Input data that is missing or inconsistent."""


#: Maximum particles per cell side: with more, the original placed them all at the centre of
#: the cell.
MAX_PARTICLES_PER_SIDE = 6


@dataclass(frozen=True)
class CaseConfig:
    """Parameters of the analysis. In brackets, the name of the variable in the Fortran."""

    title: str
    n_cells: int  # NC: cells of the PIVlab grid
    n_nodes: int  # NN: PIVlab nodes (points)
    particles_per_side: int  # NPC: particles per cell side (NPC² per cell)
    n_rows: int  # NFIL: rows of cells
    cell_width: float  # AXC [m]
    cell_height: float  # AYC [m]
    dt: float  # DT: time between images
    total_steps: int  # TOTAL_STEPS: number of PIVlab files
    print_every: int  # IMPPAS: steps between printed results
    moisture: bool  # MOISTER >= 1: the analysis carries moisture and saturation
    #: MOISTER = 2: the moisture is computed from the test images, with the settings of
    #: ``<case>.HUM``, instead of being read from the ``Moist_<n>.TXT``.
    moisture_from_images: bool
    mesh_version: int  # IVERSION: 1 = PIV-NP grid = PIVlab grid; 2 = staggered grid
    pivlab_format: int  # IPIVLAB: 1 = 4 columns (older PIVlab); anything else = 5 columns
    contour: int  # ICONTOUR: velocity correction at the boundary (0 = none)
    restart: bool  # IREC: continue from <case>.REC
    soil_density: float  # S_DENSITY [kg/m3] (not used in the computation)
    porosity: float  # POROSITY (not used in the computation)

    @property
    def n_cols(self) -> int:
        """NCH: cells per row."""
        return self.n_cells // self.n_rows

    @property
    def n_particles(self) -> int:
        """NP: particles generated on the grid."""
        cells = self.n_cells if self.mesh_version == 1 else (self.n_cols + 1) * (self.n_rows + 1)
        return cells * self.particles_per_side**2

    def validate(self) -> None:
        """Reject combinations for which the original produces meaningless results."""
        errors = []
        if self.n_rows < 1 or self.n_cells < 1:
            errors.append("NC and NFIL must be positive")
        elif self.n_cells % self.n_rows:
            errors.append(f"NC={self.n_cells} is not a multiple of NFIL={self.n_rows}")
        elif self.n_nodes != (self.n_cols + 1) * (self.n_rows + 1):
            errors.append(
                f"NN={self.n_nodes} does not match (NC/NFIL+1)*(NFIL+1)="
                f"{(self.n_cols + 1) * (self.n_rows + 1)}"
            )
        if not 1 <= self.particles_per_side <= MAX_PARTICLES_PER_SIDE:
            errors.append(f"NPC={self.particles_per_side} must be between 1 and "
                          f"{MAX_PARTICLES_PER_SIDE}")
        if self.cell_width <= 0 or self.cell_height <= 0:
            errors.append("AXC and AYC must be positive")
        if self.mesh_version not in (1, 2):
            errors.append(f"IVERSION={self.mesh_version} must be 1 or 2")
        if self.print_every < 1:
            errors.append("IMPPAS must be >= 1")
        if self.contour not in (0, 1, 2, 3):
            errors.append(f"ICONTOUR={self.contour} must be between 0 and 3")
        if errors:
            raise ConfigError("; ".join(errors))


class _ListDirectedReader:
    """Mimics the Fortran ``READ(u, 2000)`` (A80) and ``READ(u, *)`` statements."""

    def __init__(self, lines: list[str], source: str) -> None:
        self._lines = lines
        self._pos = 0
        self.source = source

    def _next_line(self) -> str:
        if self._pos >= len(self._lines):
            raise ConfigError(f"{self.source}: unexpected end of file (line {self._pos + 1})")
        line = self._lines[self._pos]
        self._pos += 1
        return line

    def text(self) -> str:
        """One line of text, truncated to 80 columns like the A80 format of the original."""
        return self._next_line()[:80].rstrip()

    def comment(self) -> str:
        """One whole comment line: the original reads nothing here, it just skips."""
        return self._next_line().rstrip()

    def values(self, count: int) -> list[str]:
        """Read ``count`` values; the rest of the last line read is discarded."""
        tokens: list[str] = []
        while len(tokens) < count:
            tokens.extend(t for t in _TOKEN_SEPARATORS.split(self._next_line().strip()) if t)
        return tokens[:count]

    def line_values(self) -> list[str]:
        """All the values of the next line."""
        return [t for t in _TOKEN_SEPARATORS.split(self._next_line().strip()) if t]

    def at_end(self) -> bool:
        return self._pos >= len(self._lines) or not any(
            line.strip() for line in self._lines[self._pos:]
        )


def _to_int(token: str, name: str) -> int:
    try:
        return int(token)
    except ValueError:
        raise ConfigError(f"{name}: expected an integer and read {token!r}") from None


def _to_float(token: str, name: str) -> float:
    try:
        return float(token.replace("D", "E").replace("d", "e"))
    except ValueError:
        raise ConfigError(f"{name}: expected a number and read {token!r}") from None


#: Field name matching every word of the comment line of block 3. It is looked up as a
#: substring, without accents and in lower case, in the order of this list. Both the English
#: names of the single format and the Spanish ones of the earlier versions are recognized,
#: so that the converter can read any of them.
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
    ("restart", "IREC"), ("rec", "IREC"),
    ("ptv", "ITR"), ("ptr", "ITR"), ("track", "ITR"), ("segui", "ITR"),
)

#: Default value of what a ``.PAR`` may not carry.
_BLOCK3_DEFAULTS = {"MOISTER": "0", "IVERSION": "1", "IPIVLAB": "1", "ICONTOUR": "0",
                    "IREC": "0", "ITR": "0"}

_PARENTHESES = re.compile(r"\([^)]*\)")


def _strip_accents(text: str) -> str:
    return "".join(c for c in unicodedata.normalize("NFD", text)
                   if unicodedata.category(c) != "Mn")


def names_in_header(header: str) -> list[str] | None:
    """Fields named by the comment line of block 3, or ``None`` when it is not understood.

    Every ``.PAR`` documents its own order on that line (``del_t total_steps salto_impresion
    v.pivnp v.pivlab moister REC PTR``), which is the only reliable way to know it: there are
    dialects carrying the same eight values in a different order.
    """
    text = _PARENTHESES.sub(" ", header)
    _, _, after = text.partition(":")  # drops the leading "BLOCK 3:"
    words = _strip_accents(after or text).lower().split()
    if not words:
        return None

    fields: list[str] = []
    for word in words:
        for key, field_name in _BLOCK3_NAMES:
            if key in word:
                if field_name in fields:  # a repeated field means this is not a name list
                    return None
                fields.append(field_name)
                break
        else:
            return None  # one unrecognized word invalidates the whole line
    return fields


#: Order of the values of block 3 in the single format.
ANALYSIS_ORDER = ("DT", "TOTAL_STEPS", "IMPPAS", "MOISTER", "IVERSION", "IPIVLAB", "ICONTOUR",
                  "IREC", "ITR")

_HOW_TO_CONVERT = ("If it comes from an earlier version of the program, move it to the single "
                   "format with:  pivnp <directory> --convert-par")


def _analysis_block(reader: _ListDirectedReader, header: str) -> dict[str, str]:
    """Block 3 of the single format: nine values in the order of :data:`ANALYSIS_ORDER`."""
    values = reader.line_values()
    if len(values) != len(ANALYSIS_ORDER):
        raise ConfigError(
            f"{reader.source}: block 3 carries {len(values)} values and the single format has "
            f"{len(ANALYSIS_ORDER)} ({' '.join(ANALYSIS_ORDER)}). {_HOW_TO_CONVERT}")
    fields = names_in_header(header)
    if fields is not None and tuple(fields) != ANALYSIS_ORDER:
        raise ConfigError(
            f"{reader.source}: the comment line of block 3 names the fields in another order "
            f"({' '.join(fields)}). {_HOW_TO_CONVERT}")
    return dict(zip(ANALYSIS_ORDER, values, strict=True))


@dataclass(frozen=True)
class RawPar:
    """The values of a ``.PAR``, as they come, already identified by name.

    It keeps the strings and not the numbers, so that converting a file to the single format
    is a matter of reordering what is there, without reformatting or rounding anything.
    """

    title: str
    geometry: list[str]  # NC NN NPC NFIL AXC AYC
    analysis: dict[str, str]  # whatever block 3 carries, by name
    density: str
    porosity: str
    #: The ``.PAR`` is in the current single format. Only the converter reads the earlier ones.
    current_dialect: bool = True

    def value(self, name: str) -> str:
        """Value from block 3, or the one assumed when the file does not carry it."""
        if name in self.analysis:
            return self.analysis[name]
        return _BLOCK3_DEFAULTS[name]


def read_par_blocks(text: str, source: str = "<PAR>") -> RawPar:
    """Read a ``.PAR`` in the single format and return its values by name.

    ``.PAR`` files from earlier versions are not read here: they are converted once with
    ``pivnp <directory> --convert-par``, and from then on there is a single input format.
    """
    reader = _ListDirectedReader(text.splitlines(), source)
    title = reader.text()
    reader.comment()
    geometry = reader.values(6)
    analysis = _analysis_block(reader, reader.comment())
    if reader.at_end():
        raise ConfigError(f"{source}: block 4 is missing, with the soil density and the "
                          f"porosity. {_HOW_TO_CONVERT}")
    reader.comment()
    density, porosity = reader.values(2)
    return RawPar(title, geometry, analysis, density, porosity)


def parse_par(text: str, source: str = "<PAR>") -> CaseConfig:
    """Interpret the contents of a ``.PAR`` file."""
    return config_from_blocks(read_par_blocks(text, source), source)


def config_from_blocks(raw: RawPar, source: str = "<PAR>") -> CaseConfig:
    """Validate the values read and turn them into the configuration of the case."""
    nc, nn, npc, nfil, axc, ayc = raw.geometry
    dt, steps, imppas = (raw.value("DT"), raw.value("TOTAL_STEPS"), raw.value("IMPPAS"))
    moister, iversion = raw.value("MOISTER"), raw.value("IVERSION")
    ipivlab, icontour = raw.value("IPIVLAB"), raw.value("ICONTOUR")
    irec, itr = raw.value("IREC"), raw.value("ITR")
    title, density, porosity = raw.title, raw.density, raw.porosity
    current_dialect = raw.current_dialect

    if _to_int(itr, "ITR") != 0:
        raise ConfigError(
            "ITR: PTV tracking particles are no longer supported. Use ITR=0 and drop block 5 "
            "of the .PAR"
        )

    irec_value = _to_int(irec, "IREC")
    if irec_value not in (0, 1):
        raise ConfigError(f"IREC={irec_value} must be 0 or 1")

    moister_value = _to_int(moister, "MOISTER")
    if not current_dialect and moister_value > 1:
        # In the earlier versions any value other than 0 turned on reading the moisture
        # files; the 2 meaning "compute it from the images" is new.
        log.info("%s: MOISTER=%d in a .PAR from an earlier version; it is taken as 1, read "
                 "the Moist_<n>.TXT", source, moister_value)
        moister_value = 1
    if moister_value not in (0, 1, 2):
        raise ConfigError(f"MOISTER={moister_value} must be 0 (no moisture), 1 (read it from "
                          "the Moist_<n>.TXT) or 2 (compute it from the images)")

    config = CaseConfig(
        title=title,
        n_cells=_to_int(nc, "NC"),
        n_nodes=_to_int(nn, "NN"),
        particles_per_side=_to_int(npc, "NPC"),
        n_rows=_to_int(nfil, "NFIL"),
        cell_width=_to_float(axc, "AXC"),
        cell_height=_to_float(ayc, "AYC"),
        dt=_to_float(dt, "DT"),
        # TOTAL_STEPS is REAL in the original: the DO loop truncates its value.
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
    """Read the case name from ``PIV-NP.TXT`` (first value of the first line)."""
    path = find_file(case_dir, CASE_INDEX_FILE)
    first_line = path.read_text(encoding="latin-1").strip().splitlines()[0].strip()
    if first_line[:1] in ("'", '"'):
        return first_line[1:].split(first_line[0], 1)[0]
    return _TOKEN_SEPARATORS.split(first_line, 1)[0]


def load_case(case_dir: Path, case_name: str | None = None) -> tuple[str, CaseConfig]:
    """Return the case name and its validated configuration."""
    case_dir = Path(case_dir)
    name = case_name or read_case_name(case_dir)
    par_path = find_file(case_dir, f"{name}.PAR")
    return name, parse_par(par_path.read_text(encoding="latin-1"), str(par_path))


def find_file(directory: Path, name: str) -> Path:
    """Look for ``name`` in ``directory``, ignoring case (like Windows)."""
    found = look_up(Path(directory), name)
    if found is None:
        raise FileNotFoundError(Path(directory) / name)
    return found


def look_up(directory: Path, name: str) -> Path | None:
    """``directory / name`` ignoring case, or ``None`` when it is not there."""
    candidate = directory / name
    if candidate.exists():
        return candidate
    if not directory.is_dir():
        return None
    lowered = name.lower()
    for entry in directory.iterdir():
        if entry.name.lower() == lowered:
            return entry
    return None


#: Subfolder of a case where the PIVlab velocity files may live instead of the case root.
VELOCITY_SUBFOLDER = "pivlab"
#: Subfolder of a case where the ``Moist_<n>.TXT`` files may live instead of the case root.
MOISTURE_SUBFOLDER = "moisture"


def find_input_file(directory: Path, name: str, subfolder: str) -> Path:
    """Look for an input file of a case, in the case root or in its own subfolder.

    A case with hundreds of steps holds hundreds of input files, which buries the handful of
    configuration files among them, so the inputs may be grouped in a subfolder:
    ``pivlab/`` for the velocity files and ``moisture/`` for the ``Moist_<n>.TXT``. The root
    is searched first, so a flat case -- every file beside the ``.PAR``, which is how every
    case was laid out before -- keeps working untouched.

    A name present in both places is an error rather than a silent choice: picking one
    without saying so is how someone ends up editing a file that is not the one being read.
    """
    directory = Path(directory)
    in_root = look_up(directory, name)
    in_subfolder = look_up(directory / subfolder, name)
    if in_root and in_subfolder:
        raise ConfigError(
            f"{name} is both in {directory} and in its {subfolder}/ subfolder, and there is "
            f"no way to tell which one you mean. Leave only one of the two."
        )
    found = in_root or in_subfolder
    if found is None:
        raise FileNotFoundError(directory / subfolder / name)
    return found
