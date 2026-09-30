from pathlib import Path

import pytest

from pivnp.config import ConfigError, find_file, load_case, parse_par, read_case_name

PAR = """\
BLOQUE 1: CALCULO DEFORMACIONES PARA ENSAYO DE LABORATORIO.
BLOQUE 2: N_cel N_nod n_p/c N_fil Ancho Alto (GEOMETRIA)
2006  2100    3    34\t0.212115\t0.212115
BLOQUE 3: del_t Total_steps Impresion Moister Version PIVlab Contour Rec Track
0.8\t149\t1\t0\t1\t1\t0\t0\t0
BLOQUE 4: Densidad Porosidad
2650.0  0.4
"""


def test_parse_reference_case():
    cfg = parse_par(PAR)
    assert cfg.title.startswith("BLOQUE 1")
    assert (cfg.n_cells, cfg.n_nodes, cfg.particles_per_side, cfg.n_rows) == (2006, 2100, 3, 34)
    assert cfg.cell_width == cfg.cell_height == 0.212115
    assert (cfg.dt, cfg.total_steps, cfg.print_every) == (0.8, 149, 1)
    assert not cfg.moisture and not cfg.restart
    assert cfg.mesh_version == 1 and cfg.pivlab_format == 1 and cfg.contour == 0
    assert (cfg.soil_density, cfg.porosity) == (2650.0, 0.4)
    assert cfg.n_cols == 59
    assert cfg.n_particles == 2006 * 9


def test_mesh_version_2_particle_count():
    cfg = parse_par(PAR.replace("0.8\t149\t1\t0\t1", "0.8\t149\t1\t0\t2"))
    assert cfg.n_particles == (59 + 1) * (34 + 1) * 9


def test_values_may_span_lines_commas_and_fortran_exponents():
    text = PAR.replace("2006  2100    3    34\t0.212115\t0.212115",
                       "2006, 2100\n3 34 0.212115 0.212115 sobra")
    text = text.replace("0.8\t149", "8.0D-1\t149.9")
    cfg = parse_par(text)
    assert (cfg.n_cells, cfg.n_nodes, cfg.particles_per_side) == (2006, 2100, 3)
    assert cfg.dt == 0.8
    assert cfg.total_steps == 149  # REAL truncado, como el DO del original


def test_old_par_format_is_accepted():
    """Los .PAR de versiones anteriores traen 3 o 4 valores en el bloque 3 y ningún bloque 4."""
    lineas = PAR.splitlines()
    antiguo = "\n".join(lineas[:4] + ["1.\t11\t1"]) + "\n"
    cfg = parse_par(antiguo)
    assert (cfg.dt, cfg.total_steps, cfg.print_every) == (1.0, 11, 1)
    assert cfg.mesh_version == 1 and cfg.pivlab_format == 1
    assert not cfg.moisture and not cfg.restart and cfg.contour == 0
    assert (cfg.soil_density, cfg.porosity) == (0.0, 0.0)

    con_humedad = "\n".join(lineas[:4] + ["1. 11 1 1"]) + "\n"
    assert parse_par(con_humedad).moisture


def test_par_with_density_inside_block_three():
    """Variante con 8 valores: DT TOTAL_STEPS IMPPAS MOISTER S_DENSITY POROSITY IVERSION PTV."""
    lineas = PAR.splitlines()
    texto = "\n".join(lineas[:4] + ["1\t149\t1\t1\t1212\t0.444\t1\t0"]) + "\n"
    cfg = parse_par(texto)
    assert (cfg.dt, cfg.total_steps, cfg.print_every) == (1.0, 149, 1)
    assert cfg.moisture and cfg.mesh_version == 1
    assert (cfg.soil_density, cfg.porosity) == (1212.0, 0.444)
    assert cfg.pivlab_format == 1 and cfg.contour == 0 and not cfg.restart


def test_moister_selects_where_the_moisture_comes_from():
    """0 sin humedad, 1 de los Moist_<n>.TXT, 2 calculada desde las imágenes."""
    def con_moister(valor):
        return parse_par(PAR.replace("0.8\t149\t1\t0\t1", f"0.8\t149\t1\t{valor}\t1"))

    assert not con_moister(0).moisture and not con_moister(0).moisture_from_images
    assert con_moister(1).moisture and not con_moister(1).moisture_from_images
    assert con_moister(2).moisture and con_moister(2).moisture_from_images
    with pytest.raises(ConfigError, match="MOISTER"):
        con_moister(3)


#: Dialectos reales del bloque 3. El orden cambió entre versiones del programa y hay dos con
#: los mismos ocho valores en distinto orden, así que se distinguen por su cabecera.
DIALECTOS = {
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
    "2024, el actual": (
        "BLOQUE 3: del_t Total_steps Impresion Moister Version PIVlab Contour Rec Track\n"
        "0.8\t149\t1\t1\t2\t2\t3\t1\t0\n"
        "BLOQUE 4: Densidad Porosidad\n2650.0  0.4\n",
        {"dt": 0.8, "mesh_version": 2, "pivlab_format": 2, "moisture": True,
         "restart": True, "contour": 3, "soil_density": 2650.0, "porosity": 0.4},
    ),
}


@pytest.mark.parametrize("nombre", list(DIALECTOS))
def test_block_three_is_read_from_its_own_header(nombre):
    """Cada .PAR documenta en su comentario qué es cada valor; hay que hacerle caso."""
    bloque, esperado = DIALECTOS[nombre]
    cfg = parse_par("\n".join(PAR.splitlines()[:3]) + "\n" + bloque)
    for campo, valor in esperado.items():
        assert getattr(cfg, campo) == valor, campo
    assert cfg.total_steps in (20, 149) and cfg.print_every == 1


def test_old_moister_two_means_read_the_files():
    """En las versiones anteriores cualquier valor distinto de 0 leía los archivos; el 2 de
    'calcular desde las imágenes' solo existe en el dialecto actual."""
    bloque, _ = DIALECTOS["centrifuga 2022"]
    cfg = parse_par("\n".join(PAR.splitlines()[:3]) + "\n" + bloque.replace("\t1\t0\t0", "\t2\t0\t0"))
    assert cfg.moisture and not cfg.moisture_from_images


def test_unreadable_header_falls_back_to_positions():
    """Con una cabecera que no nombra los campos se recurre al número de valores."""
    lineas = PAR.splitlines()
    texto = "\n".join(lineas[:4] + ["1. 11 1 1"]) + "\n"
    assert parse_par(texto).moisture
    # una cabecera con un campo repetido tampoco es fiable
    texto = texto.replace(lineas[3], "BLOQUE 3: del_t Total_steps Salto de impresión (Datos)")
    assert parse_par(texto).moisture


def test_incomplete_analysis_block_is_rejected():
    lineas = PAR.splitlines()
    with pytest.raises(ConfigError, match="bloque 3"):
        parse_par("\n".join(lineas[:4] + ["1. 11 1 0 1 1"]) + "\n")


def test_particle_count():
    cfg = parse_par(PAR)
    assert cfg.n_particles == 2006 * 3**2


@pytest.mark.parametrize(("old", "new", "message"), [
    ("2006  2100", "2005  2100", "no es múltiplo"),
    ("2006  2100", "2006  2000", "no coincide"),
    ("0.8\t149\t1\t0\t1", "0.8\t149\t1\t0\t3", "IVERSION"),
    ("0.8\t149\t1", "0.8\t149\t0", "IMPPAS"),
    ("0\t0\t0\n", "0\t2\t0\n", "IREC"),
    ("2650.0", "abc", "S_DENSITY"),
])
def test_invalid_inputs_are_rejected(old, new, message):
    with pytest.raises(ConfigError, match=message):
        parse_par(PAR.replace(old, new))


def test_truncated_file():
    with pytest.raises(ConfigError, match="fin de archivo"):
        parse_par("\n".join(PAR.splitlines()[:3]))


@pytest.mark.parametrize(("content", "name"), [
    ("zapatak\n", "zapatak"),
    ("  'caso con espacios' resto\n", "caso con espacios"),
    ("caso1, otra cosa\n", "caso1"),
])
def test_case_name(workdir: Path, content, name):
    (workdir / "piv-np.txt").write_text(content)
    assert read_case_name(workdir) == name


def test_load_case_is_case_insensitive(workdir: Path):
    (workdir / "PIV-NP.TXT").write_text("zapatak\n")
    (workdir / "zapataK.PAR").write_text(PAR)
    name, cfg = load_case(workdir)
    assert name == "zapatak" and cfg.n_cells == 2006
    assert find_file(workdir, "ZAPATAK.par").read_text() == PAR
    with pytest.raises(FileNotFoundError):
        find_file(workdir, "otro.PAR")
