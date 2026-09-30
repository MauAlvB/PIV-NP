"""Lectura de las imágenes del ensayo, selección de canal y filtro gaussiano.

El filtro reproduce el de MATLAB (``imgaussfilt``), que es lo que permite comparar los
resultados con los análisis anteriores: núcleo de tamaño ``2·ceil(2σ)+1``, separable, con
relleno por repetición del borde, y resultado redondeado al mismo tipo de la imagen.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
from numba import njit, prange

#: Canales admitidos. La numeración es la del código MATLAB del grupo.
CANALES: dict[str, int] = {"rojo": 1, "verde": 2, "azul": 3, "gris": 0}

#: Pesos de ``rgb2gray`` de MATLAB (Rec. ITU-R BT.601).
PESOS_GRIS = (0.298936021293776, 0.587043074451121, 0.114020904255103)


def numero_de_canal(canal: int | str) -> int:
    """Traduce ``1``/``"rojo"``/``"gris"``… al número interno (0 = gris)."""
    if isinstance(canal, str):
        clave = canal.strip().lower()
        if clave not in CANALES:
            raise ValueError(f"canal {canal!r} desconocido; usa uno de {sorted(CANALES)} "
                             "o el número 1, 2 o 3")
        return CANALES[clave]
    if canal not in (0, 1, 2, 3):
        raise ValueError(f"canal {canal!r} desconocido; usa 1 (rojo), 2 (verde), 3 (azul) "
                         "o 0 (gris)")
    return int(canal)


def a_gris_matlab(rgb: np.ndarray) -> np.ndarray:
    """Convierte RGB a gris igual que ``rgb2gray`` de MATLAB, redondeando a entero."""
    pesos = np.array(PESOS_GRIS, dtype=np.float64)
    gris = rgb[:, :, :3].astype(np.float64) @ pesos
    if rgb.dtype == np.uint8:
        return np.clip(np.floor(gris + 0.5), 0, 255).astype(np.uint8)
    return gris


def leer_imagen(ruta: Path, canal: int | str = "gris") -> np.ndarray:
    """Lee una imagen y devuelve un solo canal como matriz 2D.

    ``canal`` puede ser 1 (rojo), 2 (verde), 3 (azul) o 0/"gris" para la conversión a
    escala de grises.
    """
    try:
        from PIL import Image
    except ImportError:  # pragma: no cover - depende de la instalación
        raise ImportError("leer imágenes necesita Pillow: pip install pillow") from None

    numero = numero_de_canal(canal)
    with Image.open(ruta) as imagen:
        matriz = np.asarray(imagen.convert("RGB"))
    return matriz[:, :, numero - 1].copy() if numero else a_gris_matlab(matriz)


def nucleo_gaussiano(sigma: float) -> np.ndarray:
    """Núcleo 1D de MATLAB: longitud ``2·ceil(2σ)+1``, normalizado a suma 1."""
    if sigma <= 0:
        raise ValueError(f"sigma debe ser positivo y vale {sigma}")
    radio = int(math.ceil(2 * sigma))
    posiciones = np.arange(-radio, radio + 1, dtype=np.float64)
    nucleo = np.exp(-(posiciones**2) / (2.0 * sigma * sigma))
    return nucleo / nucleo.sum()


@njit(parallel=True, cache=True)
def _convolucion_horizontal(entrada, nucleo, salida):
    filas, columnas = entrada.shape
    radio = nucleo.size // 2
    for f in prange(filas):
        for c in range(columnas):
            total = 0.0
            for k in range(-radio, radio + 1):
                origen = min(max(c + k, 0), columnas - 1)  # relleno por repetición
                total += entrada[f, origen] * nucleo[k + radio]
            salida[f, c] = total


@njit(parallel=True, cache=True)
def _convolucion_vertical(entrada, nucleo, salida):
    filas, columnas = entrada.shape
    radio = nucleo.size // 2
    for f in prange(filas):
        for c in range(columnas):
            total = 0.0
            for k in range(-radio, radio + 1):
                origen = min(max(f + k, 0), filas - 1)  # relleno por repetición
                total += entrada[origen, c] * nucleo[k + radio]
            salida[f, c] = total


def desenfoque_gaussiano(imagen: np.ndarray, sigma: float) -> np.ndarray:
    """Filtro gaussiano equivalente a ``imgaussfilt(imagen, sigma)`` de MATLAB.

    Si la imagen es de enteros, el resultado se redondea y se devuelve con el mismo tipo,
    como hace MATLAB.
    """
    nucleo = nucleo_gaussiano(sigma)
    entrada = np.ascontiguousarray(imagen, dtype=np.float64)
    intermedio = np.empty_like(entrada)
    salida = np.empty_like(entrada)
    _convolucion_horizontal(entrada, nucleo, intermedio)
    _convolucion_vertical(intermedio, nucleo, salida)

    if np.issubdtype(imagen.dtype, np.integer):
        informacion = np.iinfo(imagen.dtype)
        return np.clip(np.floor(salida + 0.5), informacion.min,
                       informacion.max).astype(imagen.dtype)
    return salida
