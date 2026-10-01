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
    assert cfg.total_steps == 149  # a truncated REAL, like the DO of the original


@pytest.mark.parametrize(("text", "message"), [
    # three values in block 3 and no block 4, the oldest format
    ("1.\t11\t1", "block 3 carries 3 values"),
    # eight values, with the density and the porosity tucked into block 3
    ("1\t149\t1\t1\t1212\t0.444\t1\t0", "block 3 carries 8 values"),
])
def test_old_par_formats_are_rejected_with_instructions(text, message):
    """There is a single input format; the earlier ones are converted once and that is it."""
    old = "\n".join(PAR.splitlines()[:4] + [text]) + "\n"
    with pytest.raises(ConfigError, match=message):
        parse_par(old)
    with pytest.raises(ConfigError, match="convert-par"):
        parse_par(old)


def test_a_header_naming_the_fields_in_another_order_is_rejected():
    """Nine values but with the header of another dialect: no guessing, it warns."""
    lines = PAR.splitlines()
    text = "\n".join(lines[:3] + [
        "BLOQUE 3: del_t Total_steps Salto_de_impresión Humedad S_density porosity "
        "Version PTV IREC",
        "0.04 20 1 1 1385.46 0.506 1 0 1"] + lines[5:]) + "\n"
    with pytest.raises(ConfigError, match="another order"):
        parse_par(text)


def test_the_fourth_block_is_required():
    text = "\n".join(PAR.splitlines()[:5]) + "\n"
    with pytest.raises(ConfigError, match="block 4"):
        parse_par(text)


def test_moister_selects_where_the_moisture_comes_from():
    """0 no moisture, 1 from the Moist_<n>.TXT, 2 computed from the images."""
    def with_moister(value):
        return parse_par(PAR.replace("0.8\t149\t1\t0\t1", f"0.8\t149\t1\t{value}\t1"))

    assert not with_moister(0).moisture and not with_moister(0).moisture_from_images
    assert with_moister(1).moisture and not with_moister(1).moisture_from_images
    assert with_moister(2).moisture and with_moister(2).moisture_from_images
    with pytest.raises(ConfigError, match="MOISTER"):
        with_moister(3)


def test_particle_count():
    cfg = parse_par(PAR)
    assert cfg.n_particles == 2006 * 3**2


@pytest.mark.parametrize(("old", "new", "message"), [
    ("2006  2100", "2005  2100", "not a multiple"),
    ("2006  2100", "2006  2000", "does not match"),
    ("0.8\t149\t1\t0\t1", "0.8\t149\t1\t0\t3", "IVERSION"),
    ("0.8\t149\t1", "0.8\t149\t0", "IMPPAS"),
    ("0\t0\t0\n", "0\t2\t0\n", "IREC"),
    ("2650.0", "abc", "S_DENSITY"),
])
def test_invalid_inputs_are_rejected(old, new, message):
    with pytest.raises(ConfigError, match=message):
        parse_par(PAR.replace(old, new))


def test_truncated_file():
    with pytest.raises(ConfigError, match="end of file"):
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
