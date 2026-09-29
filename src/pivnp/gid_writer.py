"""Escritura de resultados para GiD (antes subrutina ``IMPRES_GiD``).

* ``<caso>.POST.MSH``: malla de puntos (una por partícula) con su material.
* ``<caso>.POST.RES``: resultados por partícula en cada instante impreso.

Los bloques de resultados se describen en la tabla :data:`RESULTS`; añadir o quitar un
resultado es añadir o quitar una línea. El formateo se hace en paralelo (Numba) y la
escritura a disco en un hilo aparte, para solaparla con el cálculo del paso siguiente.
"""

from __future__ import annotations

import os
from collections.abc import Callable
from concurrent.futures import Future, ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .fortran_format import format_e, format_float_rows, format_int_rows
from .solver import ACTIVE
from .state import Nodes, Particles

RES_HEADER = "GiD Post Results File 1.0"


@dataclass(frozen=True)
class ResultSpec:
    """Un bloque ``Result`` del archivo ``.POST.RES``."""

    name: str
    kind: str  # "Vector" o "Scalar" (texto de la cabecera GiD)
    values: Callable[[Particles, Nodes], np.ndarray]
    only_active: bool = True  # excluye partículas con nan_initial != 0
    integer: bool = False  # escritura list-directed de enteros (bloque NaNs)
    needs_moisture: bool = False


def _nodal_count_by_particle(p: Particles, nodes: Nodes) -> np.ndarray:
    # Reproduce el original: imprime un contador nodal indexado por número de partícula.
    counts = np.zeros(p.position.shape[0], dtype=np.int64)
    n = min(counts.size, nodes.active_count.size)
    counts[:n] = nodes.active_count[:n]
    return counts


#: Resultado "NaNs": datos que faltan alrededor de la partícula. Con IVERSION=1, cuántos de
#: los 4 nodos de su elemento no tienen medida; con IVERSION=2, si el punto PIVlab del centro
#: de su elemento no la tiene.
MISSING_DATA = ResultSpec("NaNs", "Scalar", lambda p, n: p.missing_data,
                          only_active=False, integer=True)
LEGACY_MISSING_DATA = ResultSpec("NaNs", "Scalar", _nodal_count_by_particle,
                                 only_active=False, integer=True)

#: Energía cinética: un escalar, la suma de las dos componentes. El original las escribía
#: por separado bajo una cabecera "Scalar", así que GiD solo leía la componente x.
KINETIC_ENERGY = ResultSpec("E_kinetic", "Scalar", lambda p, n: p.kinetic_energy.sum(axis=1))
LEGACY_KINETIC_ENERGY = ResultSpec("E_kinetic", "Scalar", lambda p, n: p.kinetic_energy)

RESULTS: tuple[ResultSpec, ...] = (
    ResultSpec("Displacement", "Vector", lambda p, n: p.displacement),
    ResultSpec("Inst_displacement", "Vector", lambda p, n: p.step_displacement),
    MISSING_DATA,
    ResultSpec("Velocity", "Vector", lambda p, n: p.velocity),
    ResultSpec("Acceleration", "Vector", lambda p, n: p.acceleration),
    ResultSpec("Total_strain", "Vector", lambda p, n: p.strain[:, :3]),
    ResultSpec("Inc_strain", "Vector", lambda p, n: p.strain_increment),
    ResultSpec("Equi_strain", "Scalar", lambda p, n: p.eq_strain),
    ResultSpec("Vol_strain", "Scalar", lambda p, n: p.vol_strain),
    ResultSpec("Ins_vol_strain", "Scalar", lambda p, n: p.vol_strain_increment),
    ResultSpec("In_E_strain", "Scalar", lambda p, n: p.eq_strain_increment),
    ResultSpec("E_potential", "Scalar", lambda p, n: p.potential_energy),
    KINETIC_ENERGY,
    ResultSpec("E_total", "Scalar", lambda p, n: p.total_energy),
    ResultSpec("Moisture", "Scalar", lambda p, n: p.moisture, needs_moisture=True),
    ResultSpec("Saturation", "Scalar", lambda p, n: p.saturation, needs_moisture=True),
)


def result_specs(legacy_compat: bool = False) -> tuple[ResultSpec, ...]:
    """Bloques del ``.POST.RES``.

    En modo compatibilidad, "NaNs" vuelve a ser el contador nodal del original y
    "E_kinetic" sus dos componentes.
    """
    if not legacy_compat:
        return RESULTS
    replacements = {MISSING_DATA: LEGACY_MISSING_DATA, KINETIC_ENERGY: LEGACY_KINETIC_ENERGY}
    return tuple(replacements.get(spec, spec) for spec in RESULTS)


def result_header(name: str, kind: str, time: float) -> str:
    """Cabecera ``Result`` (formatos 11 y 18 del original: nombres en campos A15)."""
    return (f"Result {name[:15]:>15} Isochrones {format_e(time)} {kind} {'OnNodes':>15}")


class GidWriter:
    """Escribe ``<caso>.POST.MSH`` y ``<caso>.POST.RES`` en ``directory``."""

    def __init__(self, directory: Path, case_name: str, eol: str = os.linesep,
                 legacy_compat: bool = False) -> None:
        self.mesh_path = Path(directory) / f"{case_name}.POST.MSH"
        self.results_path = Path(directory) / f"{case_name}.POST.RES"
        self.eol = eol
        self.specs = result_specs(legacy_compat)
        self._eol_bytes = eol.encode()
        self._file = None
        self._executor = ThreadPoolExecutor(max_workers=1, thread_name_prefix="gid-writer")
        self._pending: list[Future] = []

    # --- malla ---------------------------------------------------------------------------
    def write_mesh(self, positions: np.ndarray, nan_initial: np.ndarray) -> None:
        n = positions.shape[0]
        ids = np.arange(1, n + 1)
        mesh_name = '"Moved mesh"'
        lines = [
            "# MESH OF POINTS",
            f"MESH {mesh_name:>15} dimension {2:5d} ElemType {'Point':>15} Nnode {1:5d}",
            "Coordinates",
            "# PARTICLE_NUMBER   COORDINATE_X COORDINATE_Y",
        ]
        with open(self.mesh_path, "wb") as f:
            f.write(self._lines(lines))
            f.write(format_float_rows(ids, positions, self._eol_bytes))
            f.write(self._lines(["End Coordinates", "Elements", "# ELEMENTS MATERIAL"]))
            elements = np.stack([ids, ids, nan_initial.astype(np.int64) + 1], axis=1)
            f.write(format_int_rows(elements, 9, self._eol_bytes))
            f.write(self._lines(["End Elements"]))

    # --- resultados ----------------------------------------------------------------------
    def start_results(self, append: bool = False) -> None:
        """Abre el archivo de resultados.

        ``append`` continúa uno existente (reinicio) en vez de vaciarlo; si no existe, se
        crea con su cabecera GiD.
        """
        if append and self.results_path.exists():
            self._file = open(self.results_path, "ab")  # noqa: SIM115 (se cierra en close())
            return
        self._file = open(self.results_path, "wb")  # noqa: SIM115 (se cierra en close())
        self._submit(self._lines([RES_HEADER]))

    def write_step(self, time: float, particles: Particles, nodes: Nodes,
                   located: np.ndarray, with_moisture: bool) -> None:
        """Formatea todos los bloques del instante y los encola para escritura."""
        if self._file is None:
            raise RuntimeError("Llama a start_results() antes de write_step()")
        active = located & (particles.nan_initial == ACTIVE)
        chunks = []
        for spec in self.specs:
            if spec.needs_moisture and not with_moisture:
                continue
            mask = active if spec.only_active else located
            ids = np.flatnonzero(mask) + 1
            values = spec.values(particles, nodes)[mask]
            chunks.append(self._lines([result_header(spec.name, spec.kind, time), "Values"]))
            if spec.integer:
                chunks.append(format_int_rows(np.stack([ids, values], axis=1), 12,
                                              self._eol_bytes))
            else:
                chunks.append(format_float_rows(ids, values, self._eol_bytes))
            chunks.append(self._lines(["End values"]))
        self._submit(b"".join(chunks))

    def close(self) -> None:
        for future in self._pending:
            future.result()
        self._pending.clear()
        self._executor.shutdown(wait=True)
        if self._file is not None:
            self._file.close()
            self._file = None

    def __enter__(self) -> GidWriter:
        return self

    def __exit__(self, *exc) -> None:
        self.close()

    # --- auxiliares ----------------------------------------------------------------------
    def _lines(self, lines: list[str]) -> bytes:
        return "".join(line + self.eol for line in lines).encode("ascii")

    def _submit(self, data: bytes) -> None:
        # Como mucho dos bloques en cola: limita la memoria si el disco es lento.
        while len(self._pending) >= 2:
            self._pending.pop(0).result()
        self._pending.append(self._executor.submit(self._file.write, data))
