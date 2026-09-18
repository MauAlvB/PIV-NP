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
    assert not cfg.moisture and not cfg.restart and cfg.tracking is None
    assert cfg.mesh_version == 1 and cfg.pivlab_format == 1 and cfg.contour == 0
    assert (cfg.soil_density, cfg.porosity) == (2650.0, 0.4)
    assert cfg.n_cols == 59
    assert cfg.n_particles == cfg.n_base_particles == 2006 * 9


def test_mesh_version_2_particle_count():
    cfg = parse_par(PAR.replace("0.8\t149\t1\t0\t1", "0.8\t149\t1\t0\t2"))
    assert cfg.n_base_particles == (59 + 1) * (34 + 1) * 9


def test_values_may_span_lines_commas_and_fortran_exponents():
    text = PAR.replace("0.8\t149\t1\t0\t1\t1\t0\t0\t0", "8.0D-1, 149.9, 1\n0 1 1 0 0 0 ignorado")
    cfg = parse_par(text)
    assert cfg.dt == 0.8
    assert cfg.total_steps == 149  # REAL truncado, como el DO del original


def test_tracking_block():
    text = PAR.replace("0\t0\t0\n", "0\t0\t1\n") + "BLOQUE 5\n2 3 1 1 2 2 3 3\n"
    cfg = parse_par(text)
    assert cfg.tracking.positions() == [(2.0, 3.0), (4.0, 6.0), (6.0, 9.0)]
    assert cfg.n_particles == cfg.n_base_particles + 3


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
        parse_par("\n".join(PAR.splitlines()[:5]))


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
