"""Lectura de los archivos de PIVlab (``datos (n).txt``) y de humedad (``Moist_n.TXT``).

PIVlab guarda los puntos por columnas: x constante y y creciente *hacia abajo* (eje de la
imagen). PIV-NP numera los nodos por filas desde abajo a la izquierda y usa el eje y hacia
arriba, por eso se reordenan los nodos y se cambia el signo de la velocidad vertical.
"""

from __future__ import annotations

import re
from collections.abc import Iterator
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, replace
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:  # solo para los tipos: la humedad no depende de este módulo
    from .humedad.fuente import FuenteHumedad as MoistureFromImages

VELOCITY_PATTERN = "datos ({step}).TXT"
MOISTURE_PATTERN = "Moist_{step}.TXT"

#: Factores de conversión que PIVlab escribe en la segunda línea de cada archivo.
_XY_FACTOR = re.compile(r"px\s*->\s*m\)\s*:\s*([0-9.eE+-]+)")
_UV_FACTOR = re.compile(r"px/frame\s*->\s*m/s\)\s*:\s*([0-9.eE+-]+)")


@dataclass(frozen=True)
class Frame:
    """Datos de un instante, en el orden de los puntos del archivo PIVlab."""

    step: int
    source: Path
    u: np.ndarray  # velocidad x (NaN si PIVlab no la midió)
    v: np.ndarray  # velocidad y, eje de imagen (hacia abajo)
    moisture: np.ndarray
    saturation: np.ndarray
    #: Gris de la imagen del ensayo, ya normalizado, cuando la humedad se calcula desde las
    #: imágenes. Es el paso intermedio: la humedad sale de aplicarle el modelo, que lleva
    #: memoria de los instantes anteriores y por eso no puede calcularse aquí.
    normalized_gray: np.ndarray | None = None


def pivlab_to_node(n_cols: int, n_rows: int) -> np.ndarray:
    """Nodo PIV-NP de cada punto PIVlab (``ICONECTIVIDAD`` en el original).

    El punto ``p = col * (n_rows + 1) + fila_desde_arriba`` corresponde al nodo
    ``(n_rows - fila_desde_arriba) * (n_cols + 1) + col``.
    """
    col, row_from_top = np.divmod(np.arange((n_cols + 1) * (n_rows + 1)), n_rows + 1)
    return (n_rows - row_from_top) * (n_cols + 1) + col


def xy_factor_in_header(path: Path) -> float:
    """Metros por píxel con los que PIVlab exportó el archivo (segunda línea)."""
    lineas = Path(path).read_text(encoding="latin-1").splitlines()
    encontrado = _XY_FACTOR.search(lineas[1]) if len(lineas) >= 2 else None
    if not encontrado:
        raise ValueError(f"{path}: la cabecera no trae el factor de píxeles a metros, que "
                         "hace falta para situar los nodos en la imagen")
    return float(encontrado.group(1))


def frame_interval_in_header(path: Path) -> float | None:
    """Intervalo entre imágenes con el que se exportó el archivo, o ``None`` si no consta.

    PIVlab escribe dos factores de conversión: uno de píxeles a metros y otro de píxeles por
    imagen a metros por segundo. Su cociente es el tiempo entre imágenes que se le indicó,
    que debería coincidir con el DT del ``.PAR``: si no, los desplazamientos salen escalados.
    """
    lineas = Path(path).read_text(encoding="latin-1").splitlines()
    if len(lineas) < 2:
        return None
    xy, uv = _XY_FACTOR.search(lineas[1]), _UV_FACTOR.search(lineas[1])
    if not xy or not uv:
        return None
    try:
        divisor = float(uv.group(1))
        return float(xy.group(1)) / divisor if divisor else None
    except (ValueError, ZeroDivisionError):
        return None


def _parse_values(lines: list[str], count: int, path: Path) -> np.ndarray:
    tokens = " ".join(lines).replace(",", " ").split()
    if len(tokens) < count:
        raise ValueError(f"{path}: se esperaban {count} valores y hay {len(tokens)}")
    return np.array(tokens[:count], dtype=np.float64)


def read_velocity_file(path: Path, n_nodes: int, pivlab_format: int) -> tuple[np.ndarray, ...]:
    """Lee ``x, y, u, v`` de un archivo PIVlab (3 líneas de cabecera).

    ``pivlab_format == 1``: 4 columnas, lectura libre de 4*NN valores.
    Otro valor: 5 columnas por línea (la quinta es un indicador que no se usa).
    """
    lines = Path(path).read_text(encoding="latin-1").splitlines()[3:]
    if pivlab_format == 1:
        data = _parse_values(lines, 4 * n_nodes, path).reshape(n_nodes, 4)
    else:
        if len(lines) < n_nodes:
            raise ValueError(f"{path}: se esperaban {n_nodes} líneas y hay {len(lines)}")
        data = np.array(
            [line.replace(",", " ").split()[:4] for line in lines[:n_nodes]], dtype=np.float64
        )
    return data[:, 0], data[:, 1], data[:, 2], data[:, 3]


def read_moisture_file(path: Path, n_nodes: int) -> tuple[np.ndarray, np.ndarray]:
    """Lee humedad y saturación de ``Moist_n.TXT`` (1 línea de cabecera, 4 columnas)."""
    lines = Path(path).read_text(encoding="latin-1").splitlines()[1:]
    data = _parse_values(lines, 4 * n_nodes, path).reshape(n_nodes, 4)
    return data[:, 2].copy(), data[:, 3].copy()


class FrameSource:
    """Proveedor de instantes PIVlab con lectura anticipada en hilos."""

    def __init__(
        self,
        directory: Path,
        n_nodes: int,
        pivlab_format: int = 1,
        moisture: bool = False,
        prefetch: int = 4,
        images: MoistureFromImages | None = None,
    ) -> None:
        self.directory = Path(directory)
        self.n_nodes = n_nodes
        self.pivlab_format = pivlab_format
        #: Leer la humedad de los ``Moist_<n>.TXT``; incompatible con calcularla.
        self.moisture = moisture and images is None
        self.prefetch = max(0, prefetch)
        self.images = images
        self._index = {p.name.lower(): p for p in self.directory.iterdir()}

    def _path(self, pattern: str, step: int) -> Path:
        name = pattern.format(step=step)
        try:
            return self._index[name.lower()]
        except KeyError:
            raise FileNotFoundError(self.directory / name) from None

    def velocity_path(self, step: int) -> Path:
        """Archivo PIVlab del instante ``step``."""
        return self._path(VELOCITY_PATTERN, step)

    def mesh_in_metres(self, step: int = 1) -> tuple[np.ndarray, np.ndarray, float, np.ndarray]:
        """Nodos en metros, metros por píxel y qué nodos midió PIVlab, de un archivo."""
        path = self._path(VELOCITY_PATTERN, step)
        x, y, u, _ = read_velocity_file(path, self.n_nodes, self.pivlab_format)
        return x, y, xy_factor_in_header(path), np.isfinite(u)

    def read(self, step: int) -> Frame:
        """Lee un instante. No aplica el modelo de humedad: eso va en orden, en ``frames``."""
        path = self._path(VELOCITY_PATTERN, step)
        _, _, u, v = read_velocity_file(path, self.n_nodes, self.pivlab_format)
        if self.moisture:
            moisture, saturation = read_moisture_file(
                self._path(MOISTURE_PATTERN, step), self.n_nodes
            )
        else:
            moisture = saturation = np.zeros(self.n_nodes)
        gris = self.images.gris_normalizado(step, np.isfinite(u)) if self.images else None
        return Frame(step, path, u, v, moisture, saturation, gris)

    def _with_moisture(self, frame: Frame) -> Frame:
        """Aplica el modelo de humedad, que necesita los instantes en orden."""
        if self.images is None or frame.normalized_gray is None:
            return frame
        estado = self.images.desde_gris(frame.step, frame.normalized_gray)
        return replace(frame, moisture=estado.humedad, saturation=estado.saturacion)

    def frames(self, steps: range) -> Iterator[Frame]:
        """Devuelve los instantes en orden, leyendo los siguientes en segundo plano."""
        if self.prefetch == 0:
            yield from (self._with_moisture(self.read(step)) for step in steps)
            return
        with ThreadPoolExecutor(max_workers=self.prefetch) as pool:
            pending = [pool.submit(self.read, s) for s in steps[: self.prefetch]]
            for k in range(len(steps)):
                frame = pending[k].result()
                ahead = k + self.prefetch
                if ahead < len(steps):
                    pending.append(pool.submit(self.read, steps[ahead]))
                pending[k] = None  # libera memoria
                yield self._with_moisture(frame)
