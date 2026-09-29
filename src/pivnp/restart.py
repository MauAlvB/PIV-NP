"""Archivo de reinicio ``<caso>.REC`` (antes subrutina ``RECOM`` y lectura en PIVLAB_DATA).

Usa el formato Fortran *unformatted sequential* (cada registro va precedido y seguido de
su longitud en bytes, entero de 4 bytes), así que es compatible con los ``.REC`` generados
por el ejecutable original y viceversa.

Registros del original: NP, IVERSION, XP(NP,2), UP(NP,2), EPS(NP,4), EPSEQ(NP), NaN_P(NP).

A continuación se añaden registros con el instante, el número de paso y el estado nodal,
necesarios para que continuar un análisis dé exactamente el mismo resultado que no haberlo
interrumpido. El Fortran original los ignora, porque solo lee los siete primeros.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import BinaryIO

import numpy as np

from .state import Nodes

_MARKER = np.dtype("<i4")
EXTENSION_VERSION = 2


@dataclass(frozen=True)
class NodalState:
    """Campos nodales que se arrastran de un paso al siguiente."""

    velocity: np.ndarray  # (n_measured, 2)
    is_nan: np.ndarray  # (n_measured,)
    momentum: np.ndarray  # (n_mesh, 2)
    momentum_increment: np.ndarray  # (n_mesh, 2)
    moisture: np.ndarray  # (n_mesh,)
    saturation: np.ndarray  # (n_mesh,)

    @classmethod
    def from_nodes(cls, nodes: Nodes) -> NodalState:
        return cls(nodes.velocity, nodes.is_nan, nodes.momentum, nodes.momentum_increment,
                   nodes.moisture, nodes.saturation)

    def apply_to(self, nodes: Nodes) -> None:
        nodes.velocity[:] = self.velocity
        nodes.is_nan[:] = self.is_nan
        nodes.momentum[:] = self.momentum
        nodes.momentum_increment[:] = self.momentum_increment
        nodes.moisture[:] = self.moisture
        nodes.saturation[:] = self.saturation


@dataclass(frozen=True)
class RestartData:
    mesh_version: int
    position: np.ndarray  # (n, 2)
    displacement: np.ndarray  # (n, 2)
    strain: np.ndarray  # (n, 4)
    eq_strain: np.ndarray  # (n,)
    nan_initial: np.ndarray  # (n,)
    time: float = 0.0  # instante alcanzado
    step: int = 0  # pasos ya calculados
    nodes: NodalState | None = None
    initial_position: np.ndarray | None = None  # posición al empezar el análisis original

    @property
    def n_particles(self) -> int:
        return self.position.shape[0]


def _write_record(f: BinaryIO, array: np.ndarray) -> None:
    data = np.ascontiguousarray(array).tobytes()
    marker = np.array([len(data)], dtype=_MARKER).tobytes()
    f.write(marker + data + marker)


def _read_record(f: BinaryIO, dtype: str, optional: bool = False) -> np.ndarray | None:
    head = f.read(4)
    if len(head) < 4:
        if optional:
            return None
        raise EOFError(f"{f.name}: faltan registros en el archivo de reinicio")
    size = int(np.frombuffer(head, _MARKER)[0])
    data = f.read(size)
    tail = int(np.frombuffer(f.read(4), _MARKER)[0])
    if tail != size or len(data) != size:
        raise ValueError(f"{f.name}: registro Fortran corrupto")
    return np.frombuffer(data, dtype=dtype)


def write_restart(path: Path, data: RestartData, extended: bool = True) -> None:
    """Escribe el ``.REC``; con ``extended=False`` solo los registros del original."""
    with open(path, "wb") as f:
        _write_record(f, np.array([data.n_particles], dtype="<i4"))
        _write_record(f, np.array([data.mesh_version], dtype="<i4"))
        _write_record(f, data.position.astype("<f8"))
        _write_record(f, data.displacement.astype("<f8"))
        _write_record(f, data.strain.astype("<f8"))
        _write_record(f, data.eq_strain.astype("<f8"))
        _write_record(f, data.nan_initial.astype("<i4"))
        if not extended or data.nodes is None:
            return
        nodes = data.nodes
        _write_record(f, np.array(
            [EXTENSION_VERSION, nodes.velocity.shape[0], nodes.momentum.shape[0], data.step],
            dtype="<i4"))
        _write_record(f, np.array([data.time], dtype="<f8"))
        _write_record(f, nodes.velocity.astype("<f8"))
        _write_record(f, nodes.is_nan.astype("<i4"))
        _write_record(f, nodes.momentum.astype("<f8"))
        _write_record(f, nodes.momentum_increment.astype("<f8"))
        _write_record(f, nodes.moisture.astype("<f8"))
        _write_record(f, nodes.saturation.astype("<f8"))
        initial = data.position if data.initial_position is None else data.initial_position
        _write_record(f, initial.astype("<f8"))


def read_restart(path: Path) -> RestartData:
    with open(path, "rb") as f:
        n = int(_read_record(f, "<i4")[0])
        version = int(_read_record(f, "<i4")[0])
        position = _read_record(f, "<f8").reshape(n, 2).copy()
        displacement = _read_record(f, "<f8").reshape(n, 2).copy()
        strain = _read_record(f, "<f8").reshape(n, 4).copy()
        eq_strain = _read_record(f, "<f8").copy()
        nan_initial = _read_record(f, "<i4").astype(np.int8)

        header = _read_record(f, "<i4", optional=True)
        time, step, nodes, initial = 0.0, 0, None, None
        if header is not None:
            if header[0] != EXTENSION_VERSION:
                raise ValueError(f"{path}: versión de reinicio {header[0]} desconocida")
            n_measured, n_mesh, step = int(header[1]), int(header[2]), int(header[3])
            time = float(_read_record(f, "<f8")[0])
            nodes = NodalState(
                velocity=_read_record(f, "<f8").reshape(n_measured, 2).copy(),
                is_nan=_read_record(f, "<i4").astype(np.int8),
                momentum=_read_record(f, "<f8").reshape(n_mesh, 2).copy(),
                momentum_increment=_read_record(f, "<f8").reshape(n_mesh, 2).copy(),
                moisture=_read_record(f, "<f8").copy(),
                saturation=_read_record(f, "<f8").copy(),
            )
            initial = _read_record(f, "<f8").reshape(n, 2).copy()
        return RestartData(version, position, displacement, strain, eq_strain, nan_initial,
                           time, step, nodes, initial)
