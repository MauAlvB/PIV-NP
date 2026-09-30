"""Configuración, muestreo en los nodos y modelo de humedad."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pivnp.humedad.calibracion import Calibracion
from pivnp.humedad.configuracion import ConfiguracionError, analizar, leer
from pivnp.humedad.modelo import (
    ModeloHumedad,
    normalizar,
    referencias_por_desplazamiento,
)
from pivnp.humedad.muestreo import Registro, coordenadas_en_pixeles, muestrear

MINIMO = "CALIBRACION = cal.csv\n"


# --- configuración ----------------------------------------------------------------------
def test_valores_por_defecto_como_el_matlab():
    cfg = analizar(MINIMO)
    assert cfg.canal == 1 and cfg.sigma == 40
    assert (cfg.desplazamiento_seco, cfg.desplazamiento_saturado) == (5, -6)
    assert cfg.umbral_saturacion == 0.8 and cfg.incremental
    assert cfg.referencia_seca == "ref2.jpg"
    assert cfg.registro_es_identidad
    assert cfg.primer_instante == "igual"
    assert not cfg.redondeo_legado


def test_se_puede_pedir_el_redondeo_antiguo():
    assert analizar(MINIMO + "REDONDEO_LEGADO = si\n").redondeo_legado
    assert not analizar(MINIMO + "REDONDEO_LEGADO = 0\n").redondeo_legado


def test_lee_claves_comentarios_y_orden(workdir: Path):
    ruta = workdir / "caso.HUM"
    ruta.write_text(
        "! prueba\nSIGMA = 15   ! radio\nCANAL = verde\n\n"
        "UMBRAL_SATURACION = 0.95\nINCREMENTAL = 0\nCALIBRACION = suelo.csv\n"
        "IMAGENES = swir_{n}.jpg\nESCALA_X = 2\n", encoding="latin-1")
    cfg = leer(ruta)
    assert cfg.sigma == 15 and cfg.canal == 2
    assert cfg.umbral_saturacion == 0.95 and not cfg.incremental
    assert cfg.patron_imagenes == "swir_{n}.jpg"
    assert not cfg.registro_es_identidad
    assert cfg.ruta_imagen(workdir, 7).name == "swir_7.jpg"


def test_avisa_de_claves_desconocidas():
    cfg = analizar(MINIMO + "SIGMAA = 3\n")
    assert cfg.desconocidas == ("SIGMAA",)


@pytest.mark.parametrize(("texto", "mensaje"), [
    ("SIGMA 40\n", "CLAVE = valor"),
    (MINIMO + "SIGMA = -1\n", "SIGMA"),
    (MINIMO + "SIGMA = mucho\n", "número"),
    (MINIMO + "CANAL = 9\n", "desconocido"),
    (MINIMO + "UMBRAL_SATURACION = 1.5\n", "UMBRAL"),
    (MINIMO + "IMAGENES = vis.jpg\n", "{n}"),
    (MINIMO + "DESPLAZAMIENTO_SECO = -9\n", "DESPLAZAMIENTO_SECO"),
    (MINIMO + "PRIMER_INSTANTE = otro\n", "PRIMER_INSTANTE"),
    ("SIGMA = 40\n", "CALIBRACION"),
])
def test_configuraciones_invalidas(texto, mensaje):
    with pytest.raises(ConfiguracionError, match=mensaje):
        analizar(texto)


def test_archivo_que_no_existe(workdir: Path):
    with pytest.raises(ConfiguracionError, match="no existe"):
        leer(workdir / "falta.HUM")


# --- muestreo ---------------------------------------------------------------------------
def test_coordenadas_redondean_al_pixel_y_pasan_a_base_cero():
    x = np.array([0.0, 0.00041, 0.0006])  # con 0.0002 m/px: 0, 2.05, 3
    y = np.array([0.0, 0.0002, 0.0005])
    columna, fila = coordenadas_en_pixeles(x, y, 0.0002)
    assert columna.tolist() == [-1, 1, 2]
    assert fila.tolist() == [-1, 0, 2]


def test_el_registro_escala_y_desplaza():
    x = np.array([0.001])
    y = np.array([0.001])
    columna, fila = coordenadas_en_pixeles(x, y, 0.0001, Registro(2.0, 5.0, 0.5, -1.0))
    assert columna.tolist() == [24]  # (10 * 2 + 5) - 1
    assert fila.tolist() == [3]      # (10 * 0.5 - 1) - 1


def test_muestrea_y_descarta_lo_que_no_toca():
    imagen = np.arange(20, dtype=np.uint8).reshape(4, 5)
    columna = np.array([0, 4, 2, 99])
    fila = np.array([0, 3, 1, 1])
    con_dato = np.array([True, True, False, True])
    gris = muestrear(imagen, columna, fila, con_dato)
    assert gris[0] == 0 and gris[1] == 19
    assert np.isnan(gris[2])  # nodo sin dato de PIVlab
    assert np.isnan(gris[3])  # nodo fuera de la imagen


# --- modelo -----------------------------------------------------------------------------
def test_normalizacion_entre_las_referencias():
    referencias = referencias_por_desplazamiento(np.array([100.0, 100.0, 100.0]), 5, -6)
    # el gris de la referencia queda en 6/11 del recorrido
    valores = normalizar(np.array([105.0, 94.0, 100.0]), referencias)
    np.testing.assert_allclose(valores, [100.0, 0.0, 6 / 11 * 100])


def test_la_normalizacion_recorta_los_negativos():
    referencias = referencias_por_desplazamiento(np.array([100.0]), 5, -6)
    assert normalizar(np.array([50.0]), referencias)[0] == 0.0


def test_referencias_incoherentes():
    with pytest.raises(ValueError, match="mayor"):
        referencias_por_desplazamiento(np.array([1.0]), -5, 5)


@pytest.fixture
def calibracion() -> Calibracion:
    gris = np.array([0.0, 25.0, 50.0, 75.0, 100.0])
    return Calibracion(gris, np.array([1.0, 0.9, 0.5, 0.2, 0.0]),
                       np.array([25.0, 20.0, 10.0, 4.0, 0.0]))


def test_la_saturacion_no_baja(calibracion: Calibracion):
    modelo = ModeloHumedad(calibracion, umbral_saturacion=0.8)
    humedo = modelo.evaluar(np.array([50.0]))      # saturación 0.5
    seco = modelo.evaluar(np.array([75.0]))        # daría 0.2, pero no puede bajar
    assert humedo.saturacion[0] == pytest.approx(0.5)
    assert seco.saturacion[0] == pytest.approx(0.5)


def test_al_superar_el_umbral_se_queda_saturado(calibracion: Calibracion):
    modelo = ModeloHumedad(calibracion, umbral_saturacion=0.8)
    modelo.evaluar(np.array([25.0]))               # saturación 0.9 >= 0.8
    despues = modelo.evaluar(np.array([100.0]))    # aunque ahora mida 0
    assert despues.saturacion[0] == 1.0


def test_sin_politica_incremental_puede_secarse(calibracion: Calibracion):
    modelo = ModeloHumedad(calibracion, incremental=False)
    modelo.evaluar(np.array([25.0]))
    assert modelo.evaluar(np.array([100.0])).saturacion[0] == pytest.approx(0.0)


def test_el_umbral_es_ajustable(calibracion: Calibracion):
    estricto = ModeloHumedad(calibracion, umbral_saturacion=0.95)
    estricto.evaluar(np.array([25.0]))             # 0.9 < 0.95: no se fija en 1
    assert estricto.evaluar(np.array([100.0])).saturacion[0] == pytest.approx(0.9)


def test_los_nodos_sin_dato_se_mantienen(calibracion: Calibracion):
    modelo = ModeloHumedad(calibracion)
    estado = modelo.evaluar(np.array([np.nan, 50.0]))
    assert np.isnan(estado.saturacion[0])
    siguiente = modelo.evaluar(np.array([np.nan, 75.0]))
    assert np.isnan(siguiente.saturacion[0])
    assert siguiente.saturacion[1] == pytest.approx(0.5)  # conserva lo alcanzado


def test_umbral_invalido(calibracion: Calibracion):
    with pytest.raises(ValueError, match="umbral"):
        ModeloHumedad(calibracion, umbral_saturacion=0)
