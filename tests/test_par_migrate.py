"""Paso de los .PAR de cualquier versión al formato único.

La garantía que se comprueba es la que hace segura la conversión: leer el archivo viejo y
leer el convertido tienen que dar exactamente la misma configuración, salvo las correcciones
que el propio conversor declara.
"""

from __future__ import annotations

import dataclasses
from pathlib import Path

import pytest

from pivnp.config import ConfigError, config_from_blocks, parse_par
from pivnp.par_migrate import (
    ANALYSIS_NAMES,
    convert_file,
    convert_tree,
    read_any_par_blocks,
    to_canonical,
)

from .test_config import PAR

CABECERA_PIVLAB = ("PIVlab\nFactor de conversion (px -> m): 0.001, "
                   "(px/frame -> m/s): 0.001\nx,y,u,v\n")

#: Dialectos reales del bloque 3. El orden cambió entre versiones del programa y hay dos con
#: los mismos ocho valores en distinto orden, así que se distinguen por su cabecera.
DIALECTOS = {
    "antiguo, tres valores": (
        "BLOQUE 3: del_t  Total_steps   Salto de impresión\t(Datos)\n1.\t11\t1\n",
        {"dt": 1.0, "total_steps": 11, "mesh_version": 1, "pivlab_format": 1,
         "moisture": False, "restart": False, "contour": 0,
         "soil_density": 0.0, "porosity": 0.0},
    ),
    "antiguo, con humedad": (
        "BLOQUE 3: del_t  Total_steps   Salto de impresión\t(Datos)\n1. 11 1 1\n",
        {"dt": 1.0, "total_steps": 11, "moisture": True, "soil_density": 0.0},
    ),
    "centrifuga 2022": (
        "BLOQUE 3: del_t\t total_steps  salto_impresión  v.pivnp  v.pivlab  moister  REC  "
        "PTR\t(ANALYSIS TYPE DATA)\n2\t149\t1\t2\t2\t1\t0\t0\n"
        "BLOQUE 4: s_density(kg/m3)  initial_porosity\n3600\t0.5\n",
        {"dt": 2.0, "mesh_version": 2, "pivlab_format": 2, "moisture": True,
         "restart": False, "soil_density": 3600.0, "porosity": 0.5},
    ),
    "artículo, pruebas": (
        "BLOQUE 3: del_t Total_steps Salto_de_impresión Humedad Version  IPIVLAB  IREC  PTV\n"
        "1\t20\t1\t1\t2\t2\t1\t0\n"
        "BLOQUE 4: S_density (kg/m3) porosity \n1385.46\t0.506\n",
        {"dt": 1.0, "mesh_version": 2, "pivlab_format": 2, "moisture": True,
         "restart": True, "soil_density": 1385.46, "porosity": 0.506},
    ),
    "artículo, etapas": (
        "BLOQUE 3: del_t Total_steps Salto_de_impresión Humedad S_density (kg/m3) porosity "
        "Version PTV  IREC\n0.04\t20\t1\t1\t1385.46\t0.506\t1\t0\t1\n",
        {"dt": 0.04, "mesh_version": 1, "pivlab_format": 1, "moisture": True,
         "restart": True, "soil_density": 1385.46, "porosity": 0.506},
    ),
    "slope, ocho sin bloque 4": (
        "BLOQUE 3: del_t Total_steps Salto_de_impresión Humedad S_density (kg/m3) porosity "
        "Version PTV\n1\t149\t1\t1\t1212\t0.444\t1\t0\n",
        {"dt": 1.0, "total_steps": 149, "moisture": True, "mesh_version": 1,
         "pivlab_format": 1, "contour": 0, "restart": False,
         "soil_density": 1212.0, "porosity": 0.444},
    ),
    "2024, el actual": (
        "BLOQUE 3: del_t Total_steps Impresion Moister Version PIVlab Contour Rec Track\n"
        "0.8\t149\t1\t1\t2\t2\t3\t1\t0\n"
        "BLOQUE 4: Densidad Porosidad\n2650.0  0.4\n",
        {"dt": 0.8, "mesh_version": 2, "pivlab_format": 2, "moisture": True,
         "restart": True, "contour": 3, "soil_density": 2650.0, "porosity": 0.4},
    ),
}


def par_de(nombre: str) -> str:
    """Un .PAR completo en el dialecto pedido."""
    return "\n".join(PAR.splitlines()[:3]) + "\n" + DIALECTOS[nombre][0]


def leer_viejo(texto: str):
    """Lo que entiende el conversor de un .PAR de cualquier versión."""
    return config_from_blocks(read_any_par_blocks(texto))


@pytest.mark.parametrize("nombre", list(DIALECTOS))
def test_cada_dialecto_se_lee_por_su_cabecera(nombre):
    """Cada .PAR documenta en su comentario qué es cada valor; hay que hacerle caso."""
    esperado = DIALECTOS[nombre][1]
    cfg = leer_viejo(par_de(nombre))
    for campo, valor in esperado.items():
        assert getattr(cfg, campo) == valor, campo


def test_moister_dos_en_un_dialecto_anterior_significa_leer_los_archivos():
    texto = par_de("centrifuga 2022").replace("\t1\t0\t0", "\t2\t0\t0")
    cfg = leer_viejo(texto)
    assert cfg.moisture and not cfg.moisture_from_images


def test_un_bloque_3_que_no_es_de_ninguna_version_se_rechaza():
    texto = "\n".join(PAR.splitlines()[:4] + ["1. 11 1 0 1 1"]) + "\n"
    with pytest.raises(ConfigError, match="ninguna versión conocida"):
        leer_viejo(texto)


def escribir(workdir: Path, texto: str, nombre: str = "caso.PAR") -> Path:
    ruta = workdir / nombre
    ruta.write_text(texto, encoding="latin-1")
    return ruta


@pytest.mark.parametrize("dialecto", list(DIALECTOS))
def test_convertir_no_cambia_la_configuracion(workdir: Path, dialecto):
    texto = par_de(dialecto)
    ruta = escribir(workdir, texto)
    antes = leer_viejo(texto)

    resultado = convert_file(ruta)

    despues = parse_par(ruta.read_text(encoding="latin-1"))
    assert despues == antes, dialecto
    assert resultado.backup is not None and resultado.backup.name == "caso.PAR.orig"
    assert resultado.backup.read_text(encoding="latin-1") == texto


@pytest.mark.parametrize("dialecto", list(DIALECTOS))
def test_el_convertido_ya_no_necesita_adivinar_nada(workdir: Path, dialecto):
    """El archivo convertido se lee por su cabecera, que nombra los nueve campos."""
    ruta = escribir(workdir, par_de(dialecto))
    convert_file(ruta)
    lineas = ruta.read_text(encoding="latin-1").splitlines()
    assert lineas[3].startswith("BLOQUE 3:")
    assert lineas[3].split(":", 1)[1].split() == list(ANALYSIS_NAMES)
    assert len(lineas[4].split()) == len(ANALYSIS_NAMES)
    assert len(lineas) == 7  # título y tres bloques de dos líneas


def test_convertir_dos_veces_no_hace_nada_la_segunda(workdir: Path):
    ruta = escribir(workdir, par_de("centrifuga 2022"))
    original = ruta.read_text(encoding="latin-1")
    primero = convert_file(ruta)
    convertido = ruta.read_text(encoding="latin-1")

    segundo = convert_file(ruta)
    assert primero.changed and not segundo.changed
    assert ruta.read_text(encoding="latin-1") == convertido
    # la copia sigue siendo la del archivo de verdad, no la del ya convertido
    assert primero.backup.read_text(encoding="latin-1") == original


def test_moister_dos_de_una_version_anterior_se_escribe_como_uno(workdir: Path):
    texto = par_de("centrifuga 2022").replace("\t1\t0\t0", "\t2\t0\t0")
    ruta = escribir(workdir, texto)
    resultado = convert_file(ruta)

    assert any("MOISTER" in nota for nota in resultado.notes)
    convertido = parse_par(ruta.read_text(encoding="latin-1"))
    assert convertido.moisture and not convertido.moisture_from_images
    assert convertido == leer_viejo(texto)


def test_ipivlab_se_toma_de_los_archivos_del_caso(workdir: Path):
    """Es la única corrección que cambia la lectura: los .PAR antiguos no traen el campo."""
    ruta = escribir(workdir, par_de("artículo, etapas"))
    filas = "".join(f"0.{i},0.2,1.0,2.0,1\n" for i in range(2100))
    (workdir / "datos (1).txt").write_text(CABECERA_PIVLAB + filas)

    antes = leer_viejo(ruta.read_text(encoding="latin-1"))
    resultado = convert_file(ruta)
    despues = parse_par(ruta.read_text(encoding="latin-1"))

    assert antes.pivlab_format == 1 and despues.pivlab_format == 2
    assert any("IPIVLAB" in nota and "5 columnas" in nota for nota in resultado.notes)
    # lo demás no se toca
    assert (dataclasses.replace(despues, pivlab_format=1)
            == dataclasses.replace(antes, pivlab_format=1))


def test_sin_archivos_del_caso_no_se_inventa_el_formato(workdir: Path):
    ruta = escribir(workdir, par_de("artículo, etapas"))
    resultado = convert_file(ruta)
    assert parse_par(ruta.read_text(encoding="latin-1")).pivlab_format == 1
    assert any("no se han encontrado" in nota for nota in resultado.notes)


def test_avisa_del_dt_que_no_cuadra_con_pivlab(workdir: Path):
    ruta = escribir(workdir, par_de("artículo, pruebas"))  # DT = 1
    cabecera = ("PIVlab\nFactor de conversion (px -> m): 0.002, "
                "(px/frame -> m/s): 0.001\nx,y,u,v\n")  # intervalo 2 s
    (workdir / "datos (1).txt").write_text(cabecera + "0.1,0.2,1.0,2.0\n")
    resultado = convert_file(ruta)

    assert any("DT" in nota and "PIVlab" in nota for nota in resultado.notes)
    assert parse_par(ruta.read_text(encoding="latin-1")).dt == 1.0  # no se toca


def test_un_par_ilegible_no_se_toca(workdir: Path):
    ruta = escribir(workdir, "solo una linea\n")
    resultado = convert_file(ruta)
    assert resultado.error and not resultado.changed
    assert ruta.read_text(encoding="latin-1") == "solo una linea\n"
    assert not (workdir / "caso.PAR.orig").exists()


def test_convertir_un_arbol_entero(workdir: Path):
    for i, dialecto in enumerate(DIALECTOS):
        carpeta = workdir / f"caso{i}"
        carpeta.mkdir()
        escribir(carpeta, par_de(dialecto), f"caso{i}.PAR")
    resultados = convert_tree(workdir)
    assert len(resultados) == len(DIALECTOS)
    assert all(r.error is None for r in resultados)
    assert sum(r.changed for r in resultados) == len(DIALECTOS)


def test_los_valores_se_copian_sin_reformatear(workdir: Path):
    """Convertir no puede cambiar ni un decimal: se reordena, no se recalcula."""
    texto = par_de("artículo, pruebas").replace("1385.46", "1.38546D+03")
    ruta = escribir(workdir, texto)
    convert_file(ruta)
    assert "1.38546D+03" in ruta.read_text(encoding="latin-1")
    assert parse_par(ruta.read_text(encoding="latin-1")).soil_density == 1385.46


def test_el_formato_unico_se_puede_escribir_a_mano(workdir: Path):
    """Un caso nuevo se escribe directamente así, sin pasar por el conversor."""
    crudo = read_any_par_blocks(par_de("2024, el actual"))
    texto = to_canonical(crudo, moister=1, pivlab_format=2)
    assert parse_par(texto).pivlab_format == 2
    assert texto.count("\n") == 7
