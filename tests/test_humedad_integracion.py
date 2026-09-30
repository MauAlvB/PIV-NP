"""La humedad calculada desde las imágenes dentro de un análisis completo (MOISTER=2).

Se monta un caso mínimo: malla de una celda, cuatro nodos, tres instantes e imágenes
uniformes que se van oscureciendo, que es lo que el método interpreta como suelo mojándose.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pivnp.config import ConfigError
from pivnp.humedad.fuente import Malla, fuente_de_caso
from pivnp.simulation import Simulation, write_moisture_files

CABECERA = ("PIVlab\nFactor de conversion (px -> m): 0.001, (px/frame -> m/s): 0.001\n"
            "x [m],y [m],u [m/s],v [m/s]\n")
#: Nodos en metros. Con 0.001 m/px caen en los píxeles (1,1), (1,4), (4,1) y (4,4).
NODOS = [(0.002, 0.002), (0.002, 0.005), (0.005, 0.002), (0.005, 0.005)]
GRIS_REFERENCIA = 200
#: Tabla de calibración lineal: gris normalizado 0 (saturado) a 100 (seco).
CALIBRACION = "gris,saturacion,humedad\n0,1,30\n100,0,0\n"
PAR = """\
Caso de prueba de humedad desde imágenes
BLOQUE 2: N_cel N_nod n_p/c N_fil Ancho Alto
1  4  2  1  0.003  0.003
BLOQUE 3: del_t Total_steps Impresion Moister Version PIVlab Contour Rec Track
1.0  3  1  {moister}  1  1  0  0  0
BLOQUE 4: Densidad Porosidad
2650.0  0.4
"""
HUM = """\
! configuración de la medición de humedad
IMAGENES = img_{n}.png
CANAL = gris
SIGMA = 1
REFERENCIA_SECA = ref.png
CALIBRACION = suelo.csv
"""


def montar_caso(directorio: Path, moister: int = 2, grises=(199, 197, 194)) -> Path:
    """Escribe todos los archivos del caso y devuelve su directorio."""
    pytest.importorskip("PIL")
    from PIL import Image

    directorio.mkdir(parents=True, exist_ok=True)
    (directorio / "PIV-NP.TXT").write_text("'prueba'\n")
    (directorio / "prueba.PAR").write_text(PAR.format(moister=moister))
    (directorio / "prueba.HUM").write_text(HUM)
    (directorio / "suelo.csv").write_text(CALIBRACION)

    Image.fromarray(np.full((8, 8), GRIS_REFERENCIA, dtype=np.uint8)).save(directorio / "ref.png")
    for paso, gris in enumerate(grises, 1):
        filas = "".join(f"{x},{y},0.0,0.0\n" for x, y in NODOS)
        (directorio / f"datos ({paso}).txt").write_text(CABECERA + filas)
        Image.fromarray(np.full((8, 8), gris, dtype=np.uint8)).save(
            directorio / f"img_{paso}.png")
    return directorio


def normalizado(gris: int) -> float:
    """Gris normalizado en la banda seca-saturada, que aquí son 11 niveles."""
    return (gris - (GRIS_REFERENCIA - 6)) / 11 * 100


def test_fuente_de_caso_lee_la_configuracion_y_la_calibracion(workdir: Path):
    caso = montar_caso(workdir / "caso")
    malla = Malla(np.array([x for x, _ in NODOS]), np.array([y for _, y in NODOS]), 0.001)
    fuente = fuente_de_caso(caso, "prueba", malla, np.ones(4, dtype=bool))

    assert fuente.nodos_fuera == 0
    assert fuente.configuracion.sigma == 1 and fuente.configuracion.canal == 0
    estado = fuente.instante(1, np.ones(4, dtype=bool))
    np.testing.assert_allclose(estado.saturacion, 1 - normalizado(199) / 100, atol=1e-9)
    np.testing.assert_allclose(estado.humedad, 30 * (1 - normalizado(199) / 100), atol=1e-9)


def test_la_saturacion_no_baja_y_crece_al_oscurecerse(workdir: Path):
    caso = montar_caso(workdir / "caso", grises=(199, 197, 199))
    malla = Malla(np.array([x for x, _ in NODOS]), np.array([y for _, y in NODOS]), 0.001)
    fuente = fuente_de_caso(caso, "prueba", malla, np.ones(4, dtype=bool))
    con_dato = np.ones(4, dtype=bool)

    saturaciones = [fuente.instante(paso, con_dato).saturacion[0] for paso in (1, 2, 3)]
    assert saturaciones[1] > saturaciones[0]
    assert saturaciones[2] == saturaciones[1]  # la tercera imagen aclara, pero no se deshace


def test_analisis_completo_con_moister_2(workdir: Path):
    caso = montar_caso(workdir / "caso")
    simulacion = Simulation.from_directory(caso)

    assert simulacion.config.moisture and simulacion.config.moisture_from_images
    assert simulacion.frames.images is not None
    assert not simulacion.frames.moisture  # no se leen Moist_<n>.TXT: no existen

    simulacion.run()
    assert (simulacion.particles.moisture > 0).all()
    assert (simulacion.nodes.saturation_measured > 0).all()
    assert (caso / "prueba.POST.RES").exists()


def test_moister_2_sin_archivo_hum(workdir: Path):
    caso = montar_caso(workdir / "caso")
    (caso / "prueba.HUM").unlink()
    with pytest.raises(FileNotFoundError, match="HUM"):
        Simulation.from_directory(caso)


def test_moister_2_sin_calibracion(workdir: Path):
    caso = montar_caso(workdir / "caso")
    (caso / "suelo.csv").unlink()
    with pytest.raises(FileNotFoundError):
        Simulation.from_directory(caso)


def test_moister_1_sigue_leyendo_los_archivos(workdir: Path):
    """Con MOISTER=1 no se toca ninguna imagen: la humedad viene de los Moist_<n>.TXT."""
    caso = montar_caso(workdir / "caso", moister=1)
    for paso in (1, 2, 3):
        filas = "".join(f"{x},{y},{paso / 10},0.5\n" for x, y in NODOS)
        (caso / f"Moist_{paso}.TXT").write_text("x,y,humedad,saturacion\n" + filas)
    simulacion = Simulation.from_directory(caso)

    assert simulacion.frames.images is None and simulacion.frames.moisture
    simulacion.run()
    assert (simulacion.nodes.saturation_measured == 0.5).all()


def test_escribir_los_archivos_de_humedad(workdir: Path):
    caso = montar_caso(workdir / "caso")
    assert write_moisture_files(caso) == 3

    lineas = (caso / "Moist_2.TXT").read_text().splitlines()
    assert lineas[0] == "x_m,y_m,moisture,saturation_degree"
    assert len(lineas) == 1 + len(NODOS)
    x, y, humedad, saturacion = (float(v) for v in lineas[1].split(","))
    assert (x, y) == NODOS[0]
    assert saturacion == pytest.approx(1 - normalizado(197) / 100, abs=1e-9)
    assert humedad == pytest.approx(30 * saturacion, abs=1e-9)

    # los archivos escritos sirven luego como entrada de un análisis con MOISTER=1
    (caso / "prueba.PAR").write_text(PAR.format(moister=1))
    simulacion = Simulation.from_directory(caso)
    simulacion.run()
    assert (simulacion.particles.moisture > 0).all()


def test_moister_desconocido(workdir: Path):
    caso = montar_caso(workdir / "caso", moister=5)
    with pytest.raises(ConfigError, match="MOISTER"):
        Simulation.from_directory(caso)
