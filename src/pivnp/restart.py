"""Archivo de reinicio ``<caso>.REC`` (antes subrutina ``RECOM`` y lectura en PIVLAB_DATA).

Usa el formato Fortran *unformatted sequential* (cada registro va precedido y seguido de
su longitud en bytes, entero de 4 bytes), así que es compatible con los ``.REC`` generados
por el ejecutable original y viceversa.

Registros: NP0, IVERSION, XP(NP0,2), UP(NP0,2), EPS(NP0,4), EPSEQ(NP0), NaN_P(NP0).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import BinaryIO

import numpy as np

_MARKER = np.dtype("<i4")


@dataclass(frozen=True)
class RestartData:
    mesh_version: int
    position: np.ndarray  # (n, 2)
    displacement: np.ndarray  # (n, 2)
    strain: np.ndarray  # (n, 4)
    eq_strain: np.ndarray  # (n,)
    nan_initial: np.ndarray  # (n,)

    @property
    def n_particles(self) -> int:
        return self.position.shape[0]


def _write_record(f: BinaryIO, array: np.ndarray) -> None:
    data = np.ascontiguousarray(array).tobytes()
    marker = np.array([len(data)], dtype=_MARKER).tobytes()
    f.write(marker + data + marker)


def _read_record(f: BinaryIO, dtype: str) -> np.ndarray:
    head = f.read(4)
    if len(head) < 4:
        raise EOFError(f"{f.name}: faltan registros en el archivo de reinicio")
    size = int(np.frombuffer(head, _MARKER)[0])
    data = f.read(size)
    tail = int(np.frombuffer(f.read(4), _MARKER)[0])
    if tail != size or len(data) != size:
        raise ValueError(f"{f.name}: registro Fortran corrupto")
    return np.frombuffer(data, dtype=dtype)


def write_restart(path: Path, data: RestartData) -> None:
    with open(path, "wb") as f:
        _write_record(f, np.array([data.n_particles], dtype="<i4"))
        _write_record(f, np.array([data.mesh_version], dtype="<i4"))
        _write_record(f, data.position.astype("<f8"))
        _write_record(f, data.displacement.astype("<f8"))
        _write_record(f, data.strain.astype("<f8"))
        _write_record(f, data.eq_strain.astype("<f8"))
        _write_record(f, data.nan_initial.astype("<i4"))


def read_restart(path: Path) -> RestartData:
    with open(path, "rb") as f:
        n = int(_read_record(f, "<i4")[0])
        version = int(_read_record(f, "<i4")[0])
        return RestartData(
            mesh_version=version,
            position=_read_record(f, "<f8").reshape(n, 2).copy(),
            displacement=_read_record(f, "<f8").reshape(n, 2).copy(),
            strain=_read_record(f, "<f8").reshape(n, 4).copy(),
            eq_strain=_read_record(f, "<f8").copy(),
            nan_initial=_read_record(f, "<i4").astype(np.int8),
        )
