"""Moving the .PAR files of any version to the single format.

The guarantee checked here is the one that makes the conversion safe: reading the old file
and reading the converted one have to give exactly the same configuration, except for the
corrections the converter itself reports.
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

PIVLAB_HEADER = ("PIVlab\nConversion factor (px -> m): 0.001, "
                 "(px/frame -> m/s): 0.001\nx,y,u,v\n")

#: Real dialects of block 3. The order changed between versions of the program and two of
#: them carry the same eight values in a different order, so they are told apart by header.
DIALECTS = {
    "old, three values": (
        "BLOQUE 3: del_t  Total_steps   Salto de impresión\t(Datos)\n1.\t11\t1\n",
        {"dt": 1.0, "total_steps": 11, "mesh_version": 1, "pivlab_format": 1,
         "moisture": False, "restart": False, "contour": 0,
         "soil_density": 0.0, "porosity": 0.0},
    ),
    "old, with moisture": (
        "BLOQUE 3: del_t  Total_steps   Salto de impresión\t(Datos)\n1. 11 1 1\n",
        {"dt": 1.0, "total_steps": 11, "moisture": True, "soil_density": 0.0},
    ),
    "centrifuge 2022": (
        "BLOQUE 3: del_t\t total_steps  salto_impresión  v.pivnp  v.pivlab  moister  REC  "
        "PTR\t(ANALYSIS TYPE DATA)\n2\t149\t1\t2\t2\t1\t0\t0\n"
        "BLOQUE 4: s_density(kg/m3)  initial_porosity\n3600\t0.5\n",
        {"dt": 2.0, "mesh_version": 2, "pivlab_format": 2, "moisture": True,
         "restart": False, "soil_density": 3600.0, "porosity": 0.5},
    ),
    "paper, test runs": (
        "BLOQUE 3: del_t Total_steps Salto_de_impresión Humedad Version  IPIVLAB  IREC  PTV\n"
        "1\t20\t1\t1\t2\t2\t1\t0\n"
        "BLOQUE 4: S_density (kg/m3) porosity \n1385.46\t0.506\n",
        {"dt": 1.0, "mesh_version": 2, "pivlab_format": 2, "moisture": True,
         "restart": True, "soil_density": 1385.46, "porosity": 0.506},
    ),
    "paper, stages": (
        "BLOQUE 3: del_t Total_steps Salto_de_impresión Humedad S_density (kg/m3) porosity "
        "Version PTV  IREC\n0.04\t20\t1\t1\t1385.46\t0.506\t1\t0\t1\n",
        {"dt": 0.04, "mesh_version": 1, "pivlab_format": 1, "moisture": True,
         "restart": True, "soil_density": 1385.46, "porosity": 0.506},
    ),
    "slope, eight without block 4": (
        "BLOQUE 3: del_t Total_steps Salto_de_impresión Humedad S_density (kg/m3) porosity "
        "Version PTV\n1\t149\t1\t1\t1212\t0.444\t1\t0\n",
        {"dt": 1.0, "total_steps": 149, "moisture": True, "mesh_version": 1,
         "pivlab_format": 1, "contour": 0, "restart": False,
         "soil_density": 1212.0, "porosity": 0.444},
    ),
    "2024, the current one": (
        "BLOQUE 3: del_t Total_steps Impresion Moister Version PIVlab Contour Rec Track\n"
        "0.8\t149\t1\t1\t2\t2\t3\t1\t0\n"
        "BLOQUE 4: Densidad Porosidad\n2650.0  0.4\n",
        {"dt": 0.8, "mesh_version": 2, "pivlab_format": 2, "moisture": True,
         "restart": True, "contour": 3, "soil_density": 2650.0, "porosity": 0.4},
    ),
}


def par_of(name: str) -> str:
    """A full .PAR in the requested dialect."""
    return "\n".join(PAR.splitlines()[:3]) + "\n" + DIALECTS[name][0]


def read_old(text: str):
    """What the converter makes of a .PAR of any version."""
    return config_from_blocks(read_any_par_blocks(text))


def write(workdir: Path, text: str, name: str = "case.PAR") -> Path:
    path = workdir / name
    path.write_text(text, encoding="latin-1")
    return path


@pytest.mark.parametrize("dialect", list(DIALECTS))
def test_every_dialect_is_read_from_its_header(dialect):
    """Every .PAR documents in its comment what each value is; it has to be believed."""
    expected = DIALECTS[dialect][1]
    cfg = read_old(par_of(dialect))
    for field_name, value in expected.items():
        assert getattr(cfg, field_name) == value, field_name


def test_moister_two_in_an_earlier_dialect_means_read_the_files():
    text = par_of("centrifuge 2022").replace("\t1\t0\t0", "\t2\t0\t0")
    cfg = read_old(text)
    assert cfg.moisture and not cfg.moisture_from_images


def test_a_block_three_of_no_known_version_is_rejected():
    text = "\n".join(PAR.splitlines()[:4] + ["1. 11 1 0 1 1"]) + "\n"
    with pytest.raises(ConfigError, match="not recognized as any known version"):
        read_old(text)


@pytest.mark.parametrize("dialect", list(DIALECTS))
def test_converting_does_not_change_the_configuration(workdir: Path, dialect):
    text = par_of(dialect)
    path = write(workdir, text)
    before = read_old(text)

    result = convert_file(path)

    after = parse_par(path.read_text(encoding="latin-1"))
    assert after == before, dialect
    assert result.backup is not None and result.backup.name == "case.PAR.orig"
    assert result.backup.read_text(encoding="latin-1") == text


@pytest.mark.parametrize("dialect", list(DIALECTS))
def test_the_converted_file_no_longer_needs_any_guessing(workdir: Path, dialect):
    """The converted file is read from its header, which names the nine fields."""
    path = write(workdir, par_of(dialect))
    convert_file(path)
    lines = path.read_text(encoding="latin-1").splitlines()
    assert lines[3].startswith("BLOCK 3:")
    assert lines[3].split(":", 1)[1].split() == list(ANALYSIS_NAMES)
    assert len(lines[4].split()) == len(ANALYSIS_NAMES)
    # title and three two-line blocks; from there on, the legend
    assert not any(line.startswith("!") for line in lines[:7])
    assert all(line.startswith("!") for line in lines[7:]) and len(lines) > 7


def test_converting_twice_does_nothing_the_second_time(workdir: Path):
    path = write(workdir, par_of("centrifuge 2022"))
    original = path.read_text(encoding="latin-1")
    first = convert_file(path)
    converted = path.read_text(encoding="latin-1")

    second = convert_file(path)
    assert first.changed and not second.changed
    assert path.read_text(encoding="latin-1") == converted
    # the backup is still the real file, not the already converted one
    assert first.backup.read_text(encoding="latin-1") == original


def test_moister_two_from_an_earlier_version_is_written_as_one(workdir: Path):
    text = par_of("centrifuge 2022").replace("\t1\t0\t0", "\t2\t0\t0")
    path = write(workdir, text)
    result = convert_file(path)

    assert any("MOISTER" in note for note in result.notes)
    converted = parse_par(path.read_text(encoding="latin-1"))
    assert converted.moisture and not converted.moisture_from_images
    assert converted == read_old(text)


def test_ipivlab_is_taken_from_the_files_of_the_case(workdir: Path):
    """The only correction that changes the reading: old .PAR files lack the field."""
    path = write(workdir, par_of("paper, stages"))
    rows = "".join(f"0.{i},0.2,1.0,2.0,1\n" for i in range(2100))
    (workdir / "datos (1).txt").write_text(PIVLAB_HEADER + rows)

    before = read_old(path.read_text(encoding="latin-1"))
    result = convert_file(path)
    after = parse_par(path.read_text(encoding="latin-1"))

    assert before.pivlab_format == 1 and after.pivlab_format == 2
    assert any("IPIVLAB" in note and "5 columns" in note for note in result.notes)
    # nothing else is touched
    assert (dataclasses.replace(after, pivlab_format=1)
            == dataclasses.replace(before, pivlab_format=1))


def test_without_the_files_of_the_case_the_format_is_not_invented(workdir: Path):
    path = write(workdir, par_of("paper, stages"))
    result = convert_file(path)
    assert parse_par(path.read_text(encoding="latin-1")).pivlab_format == 1
    assert any("were not found" in note for note in result.notes)


def test_it_warns_about_a_dt_that_does_not_match_pivlab(workdir: Path):
    path = write(workdir, par_of("paper, test runs"))  # DT = 1
    header = ("PIVlab\nConversion factor (px -> m): 0.002, "
              "(px/frame -> m/s): 0.001\nx,y,u,v\n")  # 2 s interval
    (workdir / "datos (1).txt").write_text(header + "0.1,0.2,1.0,2.0\n")
    result = convert_file(path)

    assert any("DT" in note and "PIVlab" in note for note in result.notes)
    assert parse_par(path.read_text(encoding="latin-1")).dt == 1.0  # not touched


def test_an_unreadable_par_is_left_alone(workdir: Path):
    path = write(workdir, "just one line\n")
    result = convert_file(path)
    assert result.error and not result.changed
    assert path.read_text(encoding="latin-1") == "just one line\n"
    assert not (workdir / "case.PAR.orig").exists()


def test_converting_a_whole_tree(workdir: Path):
    for i, dialect in enumerate(DIALECTS):
        folder = workdir / f"case{i}"
        folder.mkdir()
        write(folder, par_of(dialect), f"case{i}.PAR")
    results = convert_tree(workdir)
    assert len(results) == len(DIALECTS)
    assert all(r.error is None for r in results)
    assert sum(r.changed for r in results) == len(DIALECTS)


def test_the_values_are_copied_without_reformatting(workdir: Path):
    """Converting cannot change a single decimal: it reorders, it does not recompute."""
    text = par_of("paper, test runs").replace("1385.46", "1.38546D+03")
    path = write(workdir, text)
    convert_file(path)
    assert "1.38546D+03" in path.read_text(encoding="latin-1")
    assert parse_par(path.read_text(encoding="latin-1")).soil_density == 1385.46


def test_the_single_format_can_be_written_by_hand(workdir: Path):
    """A new case is written straight like this, without going through the converter."""
    raw = read_any_par_blocks(par_of("2024, the current one"))
    text = to_canonical(raw, moister=1, pivlab_format=2)
    assert parse_par(text).pivlab_format == 2


def test_the_legend_at_the_end_is_not_read(workdir: Path):
    """It explains what every number means; the reader stops at block 4."""
    path = write(workdir, par_of("2024, the current one"))
    convert_file(path)
    text = path.read_text(encoding="latin-1")
    assert "moisture" in text and "2 = it is computed from the test images" in text
    assert "1 = at the nodes of the grid" in text
    assert "4-column export" in text and "5-column export" in text

    without_legend = "\n".join(text.splitlines()[:7]) + "\n"
    assert parse_par(text) == parse_par(without_legend)
    # and whatever is down there makes no difference
    assert parse_par(text + "anything at all\n1 2 3\n") == parse_par(without_legend)
