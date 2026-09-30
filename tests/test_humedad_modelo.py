"""Configuración, muestreo en los nodos y modelo de humedad."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pivnp.humedad.calibracion import Calibracion
from pivnp.humedad.configuracion import ConfiguracionError, analizar, leer
from pivnp.humedad.modelo import (
    MEDIDO,
    SIN_DATO,
    TOPE_HUMEDO,
    TOPE_SECO,
    ModeloHumedad,
    avisar_si_la_banda_es_estrecha,
    normalizar,
    referencias_globales,
    referencias_por_desplazamiento,
)
from pivnp.humedad.muestreo import Registro, coordenadas_en_pixeles, muestrear

MINIMO = "CALIBRACION = cal.csv\n"


# --- configuración ----------------------------------------------------------------------
def test_valores_por_defecto_como_el_matlab():
    cfg = analizar(MINIMO)
    # El canal por defecto es el gris, que es con el que se hizo el análisis de referencia.
    assert cfg.canal == 0 and cfg.sigma == 40
    assert (cfg.desplazamiento_seco, cfg.desplazamiento_saturado) == (5, -6)
    # el trinquete y su umbral vienen apagados: son una hipótesis, no una medida
    assert cfg.umbral_saturacion == 0.95 and not cfg.incremental
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


def test_la_banda_global_sustituye_a_la_referencia_por_nodo():
    """Dos intensidades para toda la imagen, como en el flujo SWIR."""
    cfg = analizar(MINIMO + "BANDA_SATURADA = 92\nBANDA_SECA = 132\n")
    assert cfg.banda_global and (cfg.banda_saturada, cfg.banda_seca) == (92.0, 132.0)
    assert not analizar(MINIMO).banda_global

    gris = np.array([92.0, 102.0, 132.0, 80.0])
    referencias = referencias_globales(gris.shape, cfg.banda_seca, cfg.banda_saturada)
    assert avisar_si_la_banda_es_estrecha(referencias) == 40.0
    assert referencias.seco.tolist() == [132.0] * 4
    assert referencias.saturado.tolist() == [92.0] * 4
    # una banda de 40 niveles: cada nivel de gris son 2.5 puntos de la escala
    np.testing.assert_allclose(normalizar(gris, referencias), [0.0, 25.0, 100.0, 0.0])


@pytest.mark.parametrize(("texto", "mensaje"), [
    (MINIMO + "BANDA_SECA = 132\n", "van juntas"),
    (MINIMO + "BANDA_SATURADA = 92\n", "van juntas"),
    (MINIMO + "BANDA_SECA = 90\nBANDA_SATURADA = 92\n", "BANDA_SECA"),
])
def test_bandas_globales_invalidas(texto, mensaje):
    with pytest.raises(ConfiguracionError, match=mensaje):
        analizar(texto)


def test_el_registro_admite_una_homografia():
    """Con dos cámaras que miran desde ángulos distintos hace falta la perspectiva."""
    registro = Registro(escala_x=1.1, origen_x=-90.0, escala_y=1.1, origen_y=-30.0,
                        inclinacion_xy=0.09, inclinacion_yx=0.02,
                        perspectiva_x=1e-5, perspectiva_y=8e-5)
    assert not registro.es_identidad
    x, y = np.array([0.5]), np.array([0.25])  # con 0.001 m/px: 500 y 250 píxeles
    columna, fila = coordenadas_en_pixeles(x, y, 0.001, registro)

    peso = 1e-5 * 500 + 8e-5 * 250 + 1.0
    esperada_col = round((1.1 * 500 + 0.09 * 250 - 90.0) / peso) - 1
    esperada_fil = round((0.02 * 500 + 1.1 * 250 - 30.0) / peso) - 1
    assert columna.tolist() == [esperada_col]
    assert fila.tolist() == [esperada_fil]


def test_la_homografia_se_lee_del_archivo():
    cfg = analizar(MINIMO + "ESCALA_X = 1.104\nORIGEN_X = -92.3\nINCLINACION_XY = 0.0906\n"
                            "PERSPECTIVA_Y = 8.56e-5\n")
    assert not cfg.registro_es_identidad
    registro = cfg.registro
    assert registro.escala_x == 1.104 and registro.origen_x == -92.3
    assert registro.inclinacion_xy == 0.0906 and registro.perspectiva_y == 8.56e-5
    assert registro.escala_y == 1.0 and registro.perspectiva_x == 0.0


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


def test_la_medida_viene_sin_politica_incremental(calibracion: Calibracion):
    """El trinquete es una hipótesis sobre el ensayo, no una medida: viene apagado."""
    modelo = ModeloHumedad(calibracion)
    assert not modelo.incremental and modelo.umbral == 0.95
    modelo.evaluar(np.array([25.0]))
    assert modelo.evaluar(np.array([100.0])).saturacion[0] == pytest.approx(0.0)


def test_la_saturacion_no_baja(calibracion: Calibracion):
    modelo = ModeloHumedad(calibracion, umbral_saturacion=0.8, incremental=True)
    humedo = modelo.evaluar(np.array([50.0]))      # saturación 0.5
    seco = modelo.evaluar(np.array([75.0]))        # daría 0.2, pero no puede bajar
    assert humedo.saturacion[0] == pytest.approx(0.5)
    assert seco.saturacion[0] == pytest.approx(0.5)


def test_el_trinquete_arrastra_tambien_la_humedad(calibracion: Calibracion):
    """Los dos campos salen del mismo gris, así que no pueden contradecirse.

    El código original aplicaba el trinquete solo a la saturación y recalculaba la humedad
    entera, con lo que un nodo podía quedar dado por saturado con la humedad casi a cero.
    """
    modelo = ModeloHumedad(calibracion, umbral_saturacion=0.8, incremental=True)
    humedo = modelo.evaluar(np.array([50.0]))
    seco = modelo.evaluar(np.array([75.0]))
    assert humedo.humedad[0] == pytest.approx(10.0)
    assert seco.humedad[0] == pytest.approx(10.0)   # se conserva, como la saturación
    # y lo que se publica sigue estando sobre la curva del suelo
    assert seco.humedad[0] == pytest.approx(
        np.interp(seco.saturacion[0], calibracion.saturacion[::-1],
                  calibracion.humedad[::-1]), abs=1e-9)


def test_al_superar_el_umbral_se_queda_saturado(calibracion: Calibracion):
    modelo = ModeloHumedad(calibracion, umbral_saturacion=0.8, incremental=True)
    modelo.evaluar(np.array([25.0]))               # saturación 0.9 >= 0.8
    despues = modelo.evaluar(np.array([100.0]))    # aunque ahora mida 0
    assert despues.saturacion[0] == 1.0
    assert despues.humedad[0] == pytest.approx(25.0)  # la humedad del suelo saturado


def test_sin_politica_incremental_puede_secarse(calibracion: Calibracion):
    modelo = ModeloHumedad(calibracion, incremental=False)
    modelo.evaluar(np.array([25.0]))
    assert modelo.evaluar(np.array([100.0])).saturacion[0] == pytest.approx(0.0)


def test_el_umbral_es_ajustable(calibracion: Calibracion):
    estricto = ModeloHumedad(calibracion, umbral_saturacion=0.95, incremental=True)
    estricto.evaluar(np.array([25.0]))             # 0.9 < 0.95: no se fija en 1
    assert estricto.evaluar(np.array([100.0])).saturacion[0] == pytest.approx(0.9)


def test_los_nodos_sin_dato_se_mantienen(calibracion: Calibracion):
    modelo = ModeloHumedad(calibracion, incremental=True)
    estado = modelo.evaluar(np.array([np.nan, 50.0]))
    assert np.isnan(estado.saturacion[0])
    siguiente = modelo.evaluar(np.array([np.nan, 75.0]))
    assert np.isnan(siguiente.saturacion[0])
    assert siguiente.saturacion[1] == pytest.approx(0.5)  # conserva lo alcanzado


def test_la_marca_de_calidad_distingue_medida_de_cota(calibracion: Calibracion):
    """Un nodo en un tope de la banda no es una medida, es un 'al menos' o un 'como mucho'."""
    modelo = ModeloHumedad(calibracion)
    estado = modelo.evaluar(np.array([np.nan, 0.0, 50.0, 100.0]))
    assert estado.calidad.tolist() == [SIN_DATO, TOPE_HUMEDO, MEDIDO, TOPE_SECO]
    assert estado.recortados == 2 and estado.medidos == 1


def test_umbral_invalido(calibracion: Calibracion):
    with pytest.raises(ValueError, match="umbral"):
        ModeloHumedad(calibracion, umbral_saturacion=0)
