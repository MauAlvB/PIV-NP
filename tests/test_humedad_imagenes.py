"""Lectura de imágenes, selección de canal y filtro gaussiano."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pivnp.humedad.imagenes import (
    a_gris_matlab,
    desenfoque_gaussiano,
    leer_imagen,
    nucleo_gaussiano,
    numero_de_canal,
)


# --- canales ----------------------------------------------------------------------------
@pytest.mark.parametrize(("entrada", "esperado"), [
    (1, 1), (2, 2), (3, 3), (0, 0),
    ("rojo", 1), ("Verde", 2), (" AZUL ", 3), ("gris", 0),
])
def test_numero_de_canal(entrada, esperado):
    assert numero_de_canal(entrada) == esperado


@pytest.mark.parametrize("entrada", [4, -1, "infrarrojo"])
def test_canal_desconocido(entrada):
    with pytest.raises(ValueError, match="desconocido"):
        numero_de_canal(entrada)


def test_gris_con_los_pesos_de_matlab():
    """``rgb2gray`` pondera 0.2989 R + 0.5870 G + 0.1140 B y redondea al entero."""
    rgb = np.array([[[10, 200, 30]]], dtype=np.uint8)
    esperado = round(0.298936021293776 * 10 + 0.587043074451121 * 200
                     + 0.114020904255103 * 30)
    assert a_gris_matlab(rgb)[0, 0] == esperado == 124
    for tono in (0, 128, 255):  # los grises puros se conservan
        assert a_gris_matlab(np.full((2, 2, 3), tono, dtype=np.uint8))[0, 0] == tono


def test_lee_cada_canal(workdir: Path):
    pytest.importorskip("PIL")
    from PIL import Image
    datos = np.zeros((4, 6, 3), dtype=np.uint8)
    datos[..., 0], datos[..., 1], datos[..., 2] = 10, 200, 30
    ruta = workdir / "prueba.png"
    Image.fromarray(datos).save(ruta)

    assert (leer_imagen(ruta, 1) == 10).all()
    assert (leer_imagen(ruta, "verde") == 200).all()
    assert (leer_imagen(ruta, 3) == 30).all()
    assert (leer_imagen(ruta, "gris") == 124).all()
    assert leer_imagen(ruta, 1).shape == (4, 6)


# --- núcleo -----------------------------------------------------------------------------
@pytest.mark.parametrize(("sigma", "tamano"), [(0.5, 3), (1.0, 5), (15, 61), (40, 161)])
def test_tamano_del_nucleo(sigma, tamano):
    """MATLAB usa 2·ceil(2σ)+1 puntos."""
    nucleo = nucleo_gaussiano(sigma)
    assert nucleo.size == tamano
    assert nucleo.sum() == pytest.approx(1.0)
    np.testing.assert_allclose(nucleo, nucleo[::-1], atol=1e-15)
    assert np.argmax(nucleo) == tamano // 2


def test_sigma_invalido():
    with pytest.raises(ValueError, match="positivo"):
        nucleo_gaussiano(0)


# --- filtro -----------------------------------------------------------------------------
def test_una_imagen_constante_no_cambia():
    """Comprueba a la vez la normalización del núcleo y el relleno por repetición."""
    imagen = np.full((30, 40), 173, dtype=np.uint8)
    np.testing.assert_array_equal(desenfoque_gaussiano(imagen, 5.0), imagen)


def test_coincide_con_la_convolucion_directa():
    """La versión separable debe dar lo mismo que convolucionar con el núcleo 2D."""
    rng = np.random.default_rng(0)
    imagen = rng.uniform(0, 255, size=(37, 41))
    sigma = 3.0
    nucleo = nucleo_gaussiano(sigma)
    radio = nucleo.size // 2
    nucleo2d = np.outer(nucleo, nucleo)

    filas, columnas = imagen.shape
    esperado = np.zeros_like(imagen)
    for f in range(filas):
        for c in range(columnas):
            total = 0.0
            for i in range(-radio, radio + 1):
                for j in range(-radio, radio + 1):
                    fo = min(max(f + i, 0), filas - 1)
                    co = min(max(c + j, 0), columnas - 1)
                    total += imagen[fo, co] * nucleo2d[i + radio, j + radio]
            esperado[f, c] = total
    np.testing.assert_allclose(desenfoque_gaussiano(imagen, sigma), esperado, atol=1e-9)


def test_respuesta_a_un_impulso():
    imagen = np.zeros((21, 21))
    imagen[10, 10] = 1.0
    filtrada = desenfoque_gaussiano(imagen, 2.0)
    nucleo = nucleo_gaussiano(2.0)
    np.testing.assert_allclose(filtrada[10, 10], nucleo[nucleo.size // 2] ** 2, atol=1e-12)
    assert filtrada.sum() == pytest.approx(1.0, abs=1e-9)
    np.testing.assert_allclose(filtrada, filtrada[::-1, :], atol=1e-12)


def test_conserva_el_tipo_y_redondea():
    rng = np.random.default_rng(1)
    imagen = rng.integers(0, 256, size=(20, 20), dtype=np.uint8)
    filtrada = desenfoque_gaussiano(imagen, 2.0)
    assert filtrada.dtype == np.uint8
    exacta = desenfoque_gaussiano(imagen.astype(np.float64), 2.0)
    np.testing.assert_array_equal(filtrada, np.clip(np.floor(exacta + 0.5), 0, 255))
