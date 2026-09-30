"""Curva de calibración: interpolación, lectura del archivo y recorte fuera de rango."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pivnp.humedad import Calibracion, interpolar_pchip

CALIBRACION = Path(__file__).parent / "data" / "humedad" / "calibracion_slope_rgb.csv"


# --- interpolación ----------------------------------------------------------------------
def test_pchip_reproduce_el_valor_del_caso():
    """El gris del primer instante del caso Slope_RGB da 0.149151, como en sus datos."""
    cal = Calibracion.desde_csv(CALIBRACION)
    resultado = cal.evaluar(np.array([6 / 11 * 100]))
    assert resultado.saturacion[0] == pytest.approx(0.149151, abs=5e-7)
    assert resultado.recortados == 0


def test_pchip_pasa_por_los_puntos_dados():
    x = np.array([0.0, 1.0, 3.0, 6.0])
    y = np.array([0.0, 2.0, 2.5, 9.0])
    np.testing.assert_allclose(interpolar_pchip(x, y, x), y, atol=1e-12)


def test_pchip_es_exacto_con_datos_lineales():
    x = np.array([0.0, 1.0, 2.0, 5.0])
    y = 3.0 * x - 1.0
    consulta = np.array([0.25, 1.5, 3.75, 4.99])
    np.testing.assert_allclose(interpolar_pchip(x, y, consulta), 3.0 * consulta - 1.0,
                               atol=1e-12)


def test_pchip_con_dos_puntos_es_la_recta():
    """Es lo que hace MATLAB, y es la tabla mínima que admite una calibración."""
    x = np.array([10.0, 30.0])
    y = np.array([4.0, 0.0])
    consulta = np.array([10.0, 15.0, 20.0, 30.0])
    np.testing.assert_allclose(interpolar_pchip(x, y, consulta), [4.0, 3.0, 2.0, 0.0],
                               atol=1e-12)


def test_pchip_no_se_pasa_de_los_datos():
    """A diferencia de un spline normal, no inventa máximos entre los puntos."""
    x = np.arange(6.0)
    y = np.array([0.0, 0.0, 0.0, 1.0, 1.0, 1.0])
    consulta = np.linspace(0, 5, 200)
    valores = interpolar_pchip(x, y, consulta)
    assert valores.min() >= -1e-12 and valores.max() <= 1 + 1e-12
    assert np.all(np.diff(valores) >= -1e-12)  # monótona creciente


def test_pchip_rechaza_entradas_invalidas():
    with pytest.raises(ValueError, match="ordenado"):
        interpolar_pchip(np.array([1.0, 0.0]), np.array([0.0, 1.0]), np.array([0.5]))
    with pytest.raises(ValueError, match="mismo tamaño"):
        interpolar_pchip(np.array([0.0, 1.0]), np.array([0.0]), np.array([0.5]))


# --- lectura del archivo ----------------------------------------------------------------
def test_lee_la_calibracion_del_caso():
    cal = Calibracion.desde_csv(CALIBRACION)
    assert cal.gris.size == 33
    assert cal.rango == (pytest.approx(0.00005), pytest.approx(100.0))
    # a más gris (más seco) le corresponde menos saturación
    assert cal.saturacion[0] > cal.saturacion[-1]
    # la tabla no es estrictamente monótona: en el extremo seco los valores medidos suben
    # un poco (0.0079, 0.0074, 0.0079). Por eso interesa pchip, que no se inventa un pico.
    assert np.all(np.diff(cal.saturacion) <= 1e-3)


def test_admite_comentarios_desorden_y_punto_y_coma(workdir: Path):
    ruta = workdir / "cal.csv"
    ruta.write_text("# una nota\ngris,saturacion,humedad\n50;0.5;10\n\n0;1;25\n100;0;0\n")
    cal = Calibracion.desde_csv(ruta)
    assert cal.gris.tolist() == [0.0, 50.0, 100.0]
    assert cal.saturacion.tolist() == [1.0, 0.5, 0.0]


@pytest.mark.parametrize(("contenido", "mensaje"), [
    ("gris,saturacion,humedad\n", "ninguna fila"),
    ("gris,saturacion\n0,1\n50,0.5\n", "3 columnas"),
    ("gris,saturacion,humedad\n0,1,25\n0,0.5,10\n", "sin repeticiones"),
    ("gris,saturacion,humedad\n0,1,25\n50,x,10\n", "no numéricos"),
    ("gris,saturacion,humedad\n0,1,25\n", "al menos 2 puntos"),
])
def test_rechaza_archivos_mal_formados(workdir: Path, contenido: str, mensaje: str):
    ruta = workdir / "cal.csv"
    ruta.write_text(contenido)
    with pytest.raises(ValueError, match=mensaje):
        Calibracion.desde_csv(ruta)


# --- evaluación -------------------------------------------------------------------------
@pytest.fixture
def simple() -> Calibracion:
    return Calibracion(np.array([0.0, 50.0, 100.0]), np.array([1.0, 0.5, 0.0]),
                       np.array([25.0, 10.0, 0.0]))


def test_recorta_fuera_de_rango_y_lo_cuenta(simple: Calibracion):
    resultado = simple.evaluar(np.array([-30.0, 0.0, 50.0, 100.0, 180.0]))
    assert resultado.recortados == 2
    assert resultado.saturacion[0] == pytest.approx(1.0)   # recortado al extremo húmedo
    assert resultado.saturacion[-1] == pytest.approx(0.0)  # recortado al extremo seco
    assert resultado.humedad.min() >= 0.0 and resultado.humedad.max() <= 25.0


def test_los_nodos_sin_dato_siguen_sin_dato(simple: Calibracion):
    resultado = simple.evaluar(np.array([np.nan, 50.0, np.nan]))
    assert np.isnan(resultado.saturacion[[0, 2]]).all()
    assert resultado.saturacion[1] == pytest.approx(0.5)
    assert resultado.recortados == 0


def test_la_humedad_nunca_sale_del_rango_de_la_tabla():
    """Es lo que fallaba antes: extrapolar daba humedades de −1700 %."""
    cal = Calibracion.desde_csv(CALIBRACION)
    rng = np.random.default_rng(0)
    gris = rng.uniform(-500, 500, size=10000)
    resultado = cal.evaluar(gris)
    assert resultado.humedad.min() >= cal.humedad.min() - 1e-9
    assert resultado.humedad.max() <= cal.humedad.max() + 1e-9
    assert 0.0 <= resultado.saturacion.min() and resultado.saturacion.max() <= 1.0
    assert resultado.recortados > 0
