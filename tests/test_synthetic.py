"""Validación con casos sintéticos de solución conocida.

Son campos de velocidad fabricados a mano sobre mallas pequeñas: desplazamiento sin
deformación, deformación horizontal uniforme y variable, deformación de corte y rotación de
sólido rígido. Cada uno tiene solución analítica, así que sirven para comprobar el cálculo
en sí, no solo que no cambie respecto a versiones anteriores.

Los datos y los resultados históricos vienen de los análisis del grupo (`casos.json`
guarda el nombre original de cada caso) y sus `.PAR` están en el formato antiguo, con lo
que además se comprueba que se siguen leyendo.
"""

from __future__ import annotations

import gzip
import math
import shutil
from pathlib import Path

import numpy as np
import pytest

from pivnp.particles import seed_positions
from pivnp.pivlab_io import FrameSource, pivlab_to_node
from pivnp.simulation import Simulation
from pivnp.solver import output_mask
from pivnp.vtk_export import iter_time_steps

SINTETICOS = Path(__file__).parent / "data" / "sinteticos"
CASOS = sorted(d.name for d in SINTETICOS.iterdir() if d.is_dir())
#: En la rotación todas las partículas cambian de celda, así que la comprobación celda a
#: celda no aplica; ese caso tiene su propia prueba.
CASOS_SIN_GIRO = [c for c in CASOS if not c.startswith("rotacion")]


def ejecutar(nombre: str, destino: Path) -> Simulation:
    shutil.copytree(SINTETICOS / nombre, destino, dirs_exist_ok=True,
                    ignore=shutil.ignore_patterns("expected"))
    sim = Simulation.from_directory(destino)
    sim.run()
    return sim


def activas(sim: Simulation) -> np.ndarray:
    p = sim.particles
    return output_mask(p, sim.grid, sim.config.total_steps) & (p.nan_initial == 0)


def deformacion_esperada(sim: Simulation) -> np.ndarray:
    """Deformación acumulada de cada celda, por diferencias finitas de las velocidades.

    Es un camino de cálculo independiente del que usa el programa (derivadas de las
    funciones de forma), pero equivalente para un campo bilineal.
    """
    cfg, grid = sim.config, sim.grid
    conn = pivlab_to_node(cfg.n_cols, cfg.n_rows)
    fuente = FrameSource(sim.case_dir, cfg.n_nodes, cfg.pivlab_format, prefetch=0)
    nodos = grid.cell_nodes(np.arange(grid.n_cells))  # (celdas, 4)
    total = np.zeros((grid.n_cells, 3))
    for paso in range(1, cfg.total_steps + 1):
        frame = fuente.read(paso)
        u = np.zeros(cfg.n_nodes)
        v = np.zeros(cfg.n_nodes)
        u[conn] = np.nan_to_num(frame.u)
        v[conn] = -np.nan_to_num(frame.v)  # el eje y de la imagen apunta hacia abajo
        u1, u2, u3, u4 = (u[nodos[:, k]] for k in range(4))
        v1, v2, v3, v4 = (v[nodos[:, k]] for k in range(4))
        dudx = ((u2 - u1) + (u4 - u3)) / (2 * grid.dx)
        dvdy = ((v3 - v1) + (v4 - v2)) / (2 * grid.dy)
        dudy = ((u3 - u1) + (u4 - u2)) / (2 * grid.dy)
        dvdx = ((v2 - v1) + (v4 - v3)) / (2 * grid.dx)
        total += np.stack([dudx, dvdy, dudy + dvdx], axis=1) * cfg.dt
    return total


@pytest.mark.parametrize("nombre", CASOS_SIN_GIRO)
def test_deformacion_coincide_con_la_solucion_analitica(nombre: str, workdir: Path):
    """Las partículas que no cambian de celda acumulan la deformación de esa celda."""
    sim = ejecutar(nombre, workdir / nombre)
    p, grid = sim.particles, sim.grid
    inicial = seed_positions(grid, sim.config.particles_per_side)
    celda_inicial = grid.locate(inicial)
    celda_final = grid.locate(p.position)
    quietas = activas(sim) & (celda_inicial == celda_final) & (celda_inicial >= 0)
    assert quietas.sum() >= 4

    esperada = deformacion_esperada(sim)[celda_inicial[quietas]]
    np.testing.assert_allclose(p.strain[quietas, :3], esperada, atol=1e-12)


def test_desplazamiento_sin_deformacion(workdir: Path):
    """20 cm en x y −20 cm en y, iguales en todas las partículas y sin deformar el sólido."""
    for nombre in ("desplazamiento_1P", "desplazamiento_4P"):
        sim = ejecutar(nombre, workdir / nombre)
        p = sim.particles
        vivas = activas(sim)
        np.testing.assert_allclose(p.displacement[vivas, 0], 0.20, atol=1e-12)
        np.testing.assert_allclose(p.displacement[vivas, 1], -0.20, atol=1e-12)
        assert not p.strain[vivas].any()
        assert not p.eq_strain[vivas].any()


def test_deformacion_de_corte_uniforme(workdir: Path):
    """Corte puro: γxy = 0.1 en todas las partículas, sin deformación normal ni volumétrica."""
    for nombre in ("corte_1P", "corte_4P"):
        sim = ejecutar(nombre, workdir / nombre)
        p = sim.particles
        vivas = activas(sim)
        np.testing.assert_allclose(p.strain[vivas, 2], 0.10, atol=1e-12)
        assert not p.strain[vivas, 0].any() and not p.strain[vivas, 1].any()
        assert not p.vol_strain[vivas].any()
        np.testing.assert_allclose(p.eq_strain[vivas], 0.10 / math.sqrt(3.0), rtol=1e-12)


def test_deformacion_horizontal_uniforme(workdir: Path):
    """Estiramiento uniforme: εxx = 0.24 en todas, εyy = γxy = 0, y εvol = εxx."""
    for nombre in ("horizontal_uniforme_1P", "horizontal_uniforme_4P"):
        sim = ejecutar(nombre, workdir / nombre)
        p = sim.particles
        vivas = activas(sim)
        np.testing.assert_allclose(p.strain[vivas, 0], 0.24, atol=1e-12)
        assert not p.strain[vivas, 1].any() and not p.strain[vivas, 2].any()
        np.testing.assert_allclose(p.vol_strain[vivas], 0.24, atol=1e-12)


def test_deformacion_horizontal_variable(workdir: Path):
    """Estiramiento creciente hacia la derecha: cada columna de celdas tiene su deformación."""
    for nombre in ("horizontal_1P", "horizontal_4P"):
        sim = ejecutar(nombre, workdir / nombre)
        p = sim.particles
        vivas = activas(sim)
        valores = np.unique(np.round(p.strain[vivas, 0], 9))
        np.testing.assert_allclose(valores[:3], [0.03, 0.06, 0.18], atol=1e-12)
        assert not p.strain[vivas, 1].any() and not p.strain[vivas, 2].any()


@pytest.mark.parametrize("nombre", ["rotacion_1P", "rotacion_4P"])
def test_rotacion_de_solido_rigido(nombre: str, workdir: Path):
    """Una rotación rígida no debería deformar, pero la formulación incremental sí lo hace.

    Es una limitación conocida de calcular la deformación con incrementos lineales: fija el
    error actual para detectar si algún cambio lo empeora. Corregirlo requiere una medida de
    deformación finita (gradiente de deformación), no un cambio local.
    """
    sim = ejecutar(nombre, workdir / nombre)
    p = sim.particles
    vivas = activas(sim)
    assert np.abs(p.strain[vivas, :3]).max() == pytest.approx(0.0439, abs=5e-4)
    assert np.abs(p.vol_strain[vivas]).max() < 0.045
    # el giro es simétrico respecto al centro: los desplazamientos también
    np.testing.assert_allclose(p.displacement[vivas, 0].min(), -p.displacement[vivas, 0].max(),
                               rtol=1e-9)


def test_la_deformacion_no_depende_de_las_particulas_por_celda(workdir: Path):
    """Con 1 y con 4 partículas por celda, cada celda da la misma deformación."""
    for base in ("corte", "horizontal", "horizontal_uniforme", "rotacion"):
        valores = []
        for sufijo in ("1P", "4P"):
            sim = ejecutar(f"{base}_{sufijo}", workdir / f"{base}_{sufijo}")
            vivas = activas(sim)
            valores.append(np.unique(np.round(sim.particles.strain[vivas, :3], 9), axis=0))
        comunes = min(len(valores[0]), len(valores[1]))
        assert comunes > 0
        for fila in valores[0][:comunes]:
            assert any(np.allclose(fila, otra, atol=1e-9) for otra in valores[1]), base


@pytest.mark.parametrize("nombre", CASOS)
def test_reproduce_los_resultados_historicos(nombre: str, workdir: Path):
    """Los análisis que el grupo hizo en su día se reproducen partícula a partícula."""
    sim = ejecutar(nombre, workdir / nombre)
    historico = workdir / f"{nombre}.POST.RES"
    with gzip.open(SINTETICOS / nombre / "expected" / "zapatak.POST.RES.gz", "rb") as fin:
        historico.write_bytes(fin.read())

    antiguos = {round(t, 3): b for t, b in iter_time_steps(historico)}
    nuevos = {round(t, 3): b for t, b in iter_time_steps(sim.case_dir / "zapatak.POST.RES")}
    assert set(antiguos) == set(nuevos)

    t_final = max(antiguos)
    for campo in ("Displacement", "Total_strain", "Equi_strain"):
        a, b = antiguos[t_final][campo], nuevos[t_final][campo]
        comunes = np.intersect1d(a.ids, b.ids)
        assert comunes.size >= 16
        va = a.values[np.searchsorted(a.ids, comunes)]
        vb = b.values[np.searchsorted(b.ids, comunes)]
        np.testing.assert_allclose(vb, va, atol=1e-13, err_msg=f"{nombre}: {campo}")
