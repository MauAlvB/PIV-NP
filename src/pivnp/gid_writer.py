"""Writing results for GiD (formerly the ``IMPRES_GiD`` subroutine).

* ``<case>.POST.MSH``: mesh of points (one per particle) with its material.
* ``<case>.POST.RES``: per-particle results at every printed step.

The result blocks are described in the :data:`RESULTS` table; adding or removing a result is
adding or removing one line. Formatting runs in parallel (Numba) and writing to disk happens
on a separate thread, so that it overlaps with the computation of the next step.
"""

from __future__ import annotations

import os
from collections.abc import Callable
from concurrent.futures import Future, ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .finite_strain import (
    area_change,
    deformation_gradient,
    equivalent_shear,
    green_lagrange,
    rotation_degrees,
)
from .fortran_format import format_e, format_float_rows, format_int_rows
from .solver import ACTIVE, rotation_angle, vorticity_number
from .state import Nodes, Particles

RES_HEADER = "GiD Post Results File 1.0"


@dataclass(frozen=True)
class ResultSpec:
    """One ``Result`` block of the ``.POST.RES`` file."""

    name: str
    kind: str  # "Vector" or "Scalar" (the text of the GiD header)
    values: Callable[[Particles, Nodes], np.ndarray]
    only_active: bool = True  # leaves out particles with nan_initial != 0
    integer: bool = False  # list-directed writing of integers (the NaNs block)
    needs_moisture: bool = False
    #: The original Fortran had no such block, so compatibility mode must not write it:
    #: the regression suite compares the whole file byte for byte against that version.
    absent_in_fortran: bool = False


#: Last deformation gradient computed, so the four blocks that come out of it do not each
#: redo the work. ``write_step`` evaluates them one after another on the same state.
_LAST_FINITE: tuple[int, int, tuple[np.ndarray, ...]] | None = None


def _finite(p: Particles) -> tuple[np.ndarray, ...]:
    """Strain, equivalent shear, rotation and area change from the deformation gradient.

    ``lost`` is already up to date here: ``output_mask`` refreshes it just before the step is
    written. Particles that left the grid or never had data are kept out of the differences.
    """
    global _LAST_FINITE
    n = p.position.shape[0]
    stamp = (int(p.lost[:n].sum()), float(p.displacement.sum()))
    if _LAST_FINITE is not None and _LAST_FINITE[0] == id(p.displacement) \
            and _LAST_FINITE[1] == stamp:
        return _LAST_FINITE[2]

    usable = ~p.lost[:n] & (p.nan_initial == ACTIVE)
    gradient = deformation_gradient(p, usable)
    strain = green_lagrange(gradient)
    values = (np.column_stack(strain), equivalent_shear(strain),
              rotation_degrees(gradient), area_change(gradient))
    _LAST_FINITE = (id(p.displacement), stamp, values)
    return values


def _nodal_count_by_particle(p: Particles, nodes: Nodes) -> np.ndarray:
    # Reproduces the original: it prints a nodal counter indexed by particle number.
    counts = np.zeros(p.position.shape[0], dtype=np.int64)
    n = min(counts.size, nodes.active_count.size)
    counts[:n] = nodes.active_count[:n]
    return counts


#: The "NaNs" result: data missing around the particle. With IVERSION=1, how many of the 4
#: nodes of its element have no measurement; with IVERSION=2, whether the PIVlab point at the
#: centre of its element has none.
MISSING_DATA = ResultSpec("NaNs", "Scalar", lambda p, n: p.missing_data,
                          only_active=False, integer=True)
LEGACY_MISSING_DATA = ResultSpec("NaNs", "Scalar", _nodal_count_by_particle,
                                 only_active=False, integer=True)

#: Kinetic energy: a scalar, the sum of the two components. The original wrote them
#: separately under a "Scalar" header, so GiD only read the x component.
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
    #: Rotation, which the strains cannot tell you: a rigid rotation does not deform the
    #: material, yet this incremental formulation reports strain for it. These two separate
    #: the two things.
    ResultSpec("Vorticity", "Scalar", lambda p, n: p.vorticity, absent_in_fortran=True),
    ResultSpec("Rotation", "Scalar", lambda p, n: p.rotation, absent_in_fortran=True),
    #: How the two divide the deformation up: the citable number, and the same thing
    #: bounded between 0 and 90 degrees so that it can be drawn without a singularity.
    ResultSpec("Vorticity_num", "Scalar", lambda p, n: vorticity_number(p),
               absent_in_fortran=True),
    ResultSpec("Rot_angle", "Scalar", lambda p, n: rotation_angle(p),
               absent_in_fortran=True),
    #: The same deformation measured by differentiating once over the whole analysis instead
    #: of adding a linear increment per step. Published beside the incremental results, not
    #: instead of them, so that the two can be compared on real work. See finite_strain.
    ResultSpec("Finite_strain", "Vector", lambda p, n: _finite(p)[0],
               absent_in_fortran=True),
    ResultSpec("Fin_equi_strain", "Scalar", lambda p, n: _finite(p)[1],
               absent_in_fortran=True),
    ResultSpec("Finite_rotation", "Scalar", lambda p, n: _finite(p)[2],
               absent_in_fortran=True),
    ResultSpec("Finite_area", "Scalar", lambda p, n: _finite(p)[3],
               absent_in_fortran=True),
    ResultSpec("E_potential", "Scalar", lambda p, n: p.potential_energy),
    KINETIC_ENERGY,
    ResultSpec("E_total", "Scalar", lambda p, n: p.total_energy),
    ResultSpec("Moisture", "Scalar", lambda p, n: p.moisture, needs_moisture=True),
    ResultSpec("Saturation", "Scalar", lambda p, n: p.saturation, needs_moisture=True),
)


def result_specs(legacy_compat: bool = False) -> tuple[ResultSpec, ...]:
    """Blocks of the ``.POST.RES``.

    In compatibility mode, "NaNs" goes back to being the nodal counter of the original,
    "E_kinetic" to its two components, and the blocks the original never wrote are left out,
    so that the file stays identical to the one that version produced.
    """
    if not legacy_compat:
        return RESULTS
    replacements = {MISSING_DATA: LEGACY_MISSING_DATA, KINETIC_ENERGY: LEGACY_KINETIC_ENERGY}
    return tuple(replacements.get(spec, spec) for spec in RESULTS
                 if not spec.absent_in_fortran)


def result_header(name: str, kind: str, time: float) -> str:
    """``Result`` header (formats 11 and 18 of the original: names in A15 fields)."""
    return (f"Result {name[:15]:>15} Isochrones {format_e(time)} {kind} {'OnNodes':>15}")


class GidWriter:
    """Writes ``<case>.POST.MSH`` and ``<case>.POST.RES`` into ``directory``."""

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

    # --- mesh ----------------------------------------------------------------------------
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

    # --- results -------------------------------------------------------------------------
    def start_results(self, append: bool = False) -> None:
        """Open the results file.

        ``append`` continues an existing one (a restart) instead of emptying it; if there is
        none, it is created with its GiD header.
        """
        if append and self.results_path.exists():
            self._file = open(self.results_path, "ab")  # noqa: SIM115 (closed in close())
            return
        self._file = open(self.results_path, "wb")  # noqa: SIM115 (closed in close())
        self._submit(self._lines([RES_HEADER]))

    def write_step(self, time: float, particles: Particles, nodes: Nodes,
                   located: np.ndarray, with_moisture: bool) -> None:
        """Format every block of the instant and queue them to be written."""
        if self._file is None:
            raise RuntimeError("call start_results() before write_step()")
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

    # --- helpers -------------------------------------------------------------------------
    def _lines(self, lines: list[str]) -> bytes:
        return "".join(line + self.eol for line in lines).encode("ascii")

    def _submit(self, data: bytes) -> None:
        # At most two blocks queued: this caps the memory used when the disk is slow.
        while len(self._pending) >= 2:
            self._pending.pop(0).result()
        self._pending.append(self._executor.submit(self._file.write, data))
