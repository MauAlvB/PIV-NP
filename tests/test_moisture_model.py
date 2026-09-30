"""Settings, sampling at the nodes and the moisture model."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pivnp.moisture.calibration import Calibration
from pivnp.moisture.model import (
    DRY_LIMIT,
    MEASURED,
    NO_DATA,
    WET_LIMIT,
    MoistureModel,
    global_references,
    normalize,
    references_from_offsets,
    warn_if_band_is_narrow,
)
from pivnp.moisture.sampling import Registration, pixel_coordinates, sample
from pivnp.moisture.settings import SettingsError, parse, read

MINIMAL = "CALIBRATION = cal.csv\n"


# --- settings -----------------------------------------------------------------------------
def test_default_values():
    cfg = parse(MINIMAL)
    # the default channel is gray, which is the one the reference analysis used
    assert cfg.channel == 0 and cfg.sigma == 40
    assert (cfg.dry_offset, cfg.saturated_offset) == (5, -6)
    # the ratchet and its threshold come off: they are a hypothesis, not a measurement
    assert cfg.saturation_threshold == 0.95 and not cfg.incremental
    assert cfg.dry_reference == "ref2.jpg"
    assert cfg.registration_is_identity
    assert cfg.first_step == "same"
    assert not cfg.legacy_rounding


def test_the_old_rounding_can_be_asked_for():
    assert parse(MINIMAL + "LEGACY_ROUNDING = yes\n").legacy_rounding
    assert not parse(MINIMAL + "LEGACY_ROUNDING = 0\n").legacy_rounding


def test_reads_keys_comments_and_any_order(workdir: Path):
    path = workdir / "case.HUM"
    path.write_text(
        "! a test\nSIGMA = 15   ! radius\nCHANNEL = green\n\n"
        "SATURATION_THRESHOLD = 0.95\nINCREMENTAL = 0\nCALIBRATION = soil.csv\n"
        "IMAGES = swir_{n}.jpg\nSCALE_X = 2\n", encoding="latin-1")
    cfg = read(path)
    assert cfg.sigma == 15 and cfg.channel == 2
    assert cfg.saturation_threshold == 0.95 and not cfg.incremental
    assert cfg.image_pattern == "swir_{n}.jpg"
    assert not cfg.registration_is_identity
    assert cfg.image_path(workdir, 7).name == "swir_7.jpg"


def test_warns_about_unknown_keys():
    cfg = parse(MINIMAL + "SIGMAA = 3\n")
    assert cfg.unknown_keys == ("SIGMAA",)


@pytest.mark.parametrize(("text", "message"), [
    ("SIGMA 40\n", "KEY = value"),
    (MINIMAL + "SIGMA = -1\n", "SIGMA"),
    (MINIMAL + "SIGMA = a lot\n", "number"),
    (MINIMAL + "CHANNEL = 9\n", "unknown"),
    (MINIMAL + "SATURATION_THRESHOLD = 1.5\n", "THRESHOLD"),
    (MINIMAL + "IMAGES = vis.jpg\n", "{n}"),
    (MINIMAL + "DRY_OFFSET = -9\n", "DRY_OFFSET"),
    (MINIMAL + "FIRST_STEP = other\n", "FIRST_STEP"),
    ("SIGMA = 40\n", "CALIBRATION"),
])
def test_invalid_settings(text, message):
    with pytest.raises(SettingsError, match=message):
        parse(text)


def test_a_file_that_does_not_exist(workdir: Path):
    with pytest.raises(SettingsError, match="does not exist"):
        read(workdir / "missing.HUM")


# --- sampling -----------------------------------------------------------------------------
def test_coordinates_round_to_the_pixel_and_become_zero_based():
    x = np.array([0.0, 0.00041, 0.0006])  # with 0.0002 m/px: 0, 2.05, 3
    y = np.array([0.0, 0.0002, 0.0005])
    column, row = pixel_coordinates(x, y, 0.0002)
    assert column.tolist() == [-1, 1, 2]
    assert row.tolist() == [-1, 0, 2]


def test_the_registration_scales_and_shifts():
    x = np.array([0.001])
    y = np.array([0.001])
    column, row = pixel_coordinates(x, y, 0.0001, Registration(2.0, 5.0, 0.5, -1.0))
    assert column.tolist() == [24]  # (10 * 2 + 5) - 1
    assert row.tolist() == [3]      # (10 * 0.5 - 1) - 1


def test_the_global_band_replaces_the_per_node_reference():
    """Two intensities for the whole image, as in the SWIR flow."""
    cfg = parse(MINIMAL + "SATURATED_BAND = 92\nDRY_BAND = 132\n")
    assert cfg.has_global_band and (cfg.saturated_band, cfg.dry_band) == (92.0, 132.0)
    assert not parse(MINIMAL).has_global_band

    gray = np.array([92.0, 102.0, 132.0, 80.0])
    references = global_references(gray.shape, cfg.dry_band, cfg.saturated_band)
    assert warn_if_band_is_narrow(references) == 40.0
    assert references.dry.tolist() == [132.0] * 4
    assert references.saturated.tolist() == [92.0] * 4
    # a 40-level band: every gray level is 2.5 points of the scale
    np.testing.assert_allclose(normalize(gray, references), [0.0, 25.0, 100.0, 0.0])


@pytest.mark.parametrize(("text", "message"), [
    (MINIMAL + "DRY_BAND = 132\n", "go together"),
    (MINIMAL + "SATURATED_BAND = 92\n", "go together"),
    (MINIMAL + "DRY_BAND = 90\nSATURATED_BAND = 92\n", "DRY_BAND"),
])
def test_invalid_global_bands(text, message):
    with pytest.raises(SettingsError, match=message):
        parse(text)


def test_the_registration_admits_a_homography():
    """With two cameras looking from different angles the perspective is needed."""
    registration = Registration(scale_x=1.1, offset_x=-90.0, scale_y=1.1, offset_y=-30.0,
                                shear_xy=0.09, shear_yx=0.02,
                                perspective_x=1e-5, perspective_y=8e-5)
    assert not registration.is_identity
    x, y = np.array([0.5]), np.array([0.25])  # with 0.001 m/px: 500 and 250 pixels
    column, row = pixel_coordinates(x, y, 0.001, registration)

    weight = 1e-5 * 500 + 8e-5 * 250 + 1.0
    expected_column = round((1.1 * 500 + 0.09 * 250 - 90.0) / weight) - 1
    expected_row = round((0.02 * 500 + 1.1 * 250 - 30.0) / weight) - 1
    assert column.tolist() == [expected_column]
    assert row.tolist() == [expected_row]


def test_the_homography_is_read_from_the_file():
    cfg = parse(MINIMAL + "SCALE_X = 1.104\nOFFSET_X = -92.3\nSHEAR_XY = 0.0906\n"
                          "PERSPECTIVE_Y = 8.56e-5\n")
    assert not cfg.registration_is_identity
    registration = cfg.registration
    assert registration.scale_x == 1.104 and registration.offset_x == -92.3
    assert registration.shear_xy == 0.0906 and registration.perspective_y == 8.56e-5
    assert registration.scale_y == 1.0 and registration.perspective_x == 0.0


def test_samples_and_drops_what_it_should_not_touch():
    image = np.arange(20, dtype=np.uint8).reshape(4, 5)
    column = np.array([0, 4, 2, 99])
    row = np.array([0, 3, 1, 1])
    has_data = np.array([True, True, False, True])
    gray = sample(image, column, row, has_data)
    assert gray[0] == 0 and gray[1] == 19
    assert np.isnan(gray[2])  # node without PIVlab data
    assert np.isnan(gray[3])  # node outside the image


# --- model --------------------------------------------------------------------------------
def test_normalization_between_the_references():
    references = references_from_offsets(np.array([100.0, 100.0, 100.0]), 5, -6)
    # the gray of the reference sits at 6/11 of the span
    values = normalize(np.array([105.0, 94.0, 100.0]), references)
    np.testing.assert_allclose(values, [100.0, 0.0, 6 / 11 * 100])


def test_the_normalization_clamps_negatives():
    references = references_from_offsets(np.array([100.0]), 5, -6)
    assert normalize(np.array([50.0]), references)[0] == 0.0


def test_inconsistent_references():
    with pytest.raises(ValueError, match="greater"):
        references_from_offsets(np.array([1.0]), -5, 5)


@pytest.fixture
def calibration() -> Calibration:
    gray = np.array([0.0, 25.0, 50.0, 75.0, 100.0])
    return Calibration(gray, np.array([1.0, 0.9, 0.5, 0.2, 0.0]),
                       np.array([25.0, 20.0, 10.0, 4.0, 0.0]))


def test_the_measurement_comes_without_the_incremental_policy(calibration: Calibration):
    """The ratchet is a hypothesis about the test, not a measurement: it comes off."""
    model = MoistureModel(calibration)
    assert not model.incremental and model.threshold == 0.95
    model.evaluate(np.array([25.0]))
    assert model.evaluate(np.array([100.0])).saturation[0] == pytest.approx(0.0)


def test_the_saturation_does_not_go_down(calibration: Calibration):
    model = MoistureModel(calibration, saturation_threshold=0.8, incremental=True)
    wet = model.evaluate(np.array([50.0]))      # saturation 0.5
    dry = model.evaluate(np.array([75.0]))      # would give 0.2, but it cannot go down
    assert wet.saturation[0] == pytest.approx(0.5)
    assert dry.saturation[0] == pytest.approx(0.5)


def test_the_ratchet_carries_the_water_content_along(calibration: Calibration):
    """Both fields come out of the same gray, so they cannot contradict each other.

    The original code applied the ratchet to the saturation only and recomputed the water
    content from scratch, so a node could end up taken as saturated with its water content
    down to almost zero.
    """
    model = MoistureModel(calibration, saturation_threshold=0.8, incremental=True)
    wet = model.evaluate(np.array([50.0]))
    dry = model.evaluate(np.array([75.0]))
    assert wet.moisture[0] == pytest.approx(10.0)
    assert dry.moisture[0] == pytest.approx(10.0)   # preserved, like the saturation
    # and what is published still sits on the curve of the soil
    assert dry.moisture[0] == pytest.approx(
        np.interp(dry.saturation[0], calibration.saturation[::-1],
                  calibration.moisture[::-1]), abs=1e-9)


def test_past_the_threshold_it_stays_saturated(calibration: Calibration):
    model = MoistureModel(calibration, saturation_threshold=0.8, incremental=True)
    model.evaluate(np.array([25.0]))               # saturation 0.9 >= 0.8
    after = model.evaluate(np.array([100.0]))      # even though it now measures 0
    assert after.saturation[0] == 1.0
    assert after.moisture[0] == pytest.approx(25.0)  # the water content of saturated soil


def test_without_the_incremental_policy_it_can_dry(calibration: Calibration):
    model = MoistureModel(calibration, incremental=False)
    model.evaluate(np.array([25.0]))
    assert model.evaluate(np.array([100.0])).saturation[0] == pytest.approx(0.0)


def test_the_threshold_is_adjustable(calibration: Calibration):
    strict = MoistureModel(calibration, saturation_threshold=0.95, incremental=True)
    strict.evaluate(np.array([25.0]))              # 0.9 < 0.95: it is not pinned to 1
    assert strict.evaluate(np.array([100.0])).saturation[0] == pytest.approx(0.9)


def test_nodes_without_data_are_kept(calibration: Calibration):
    model = MoistureModel(calibration, incremental=True)
    state = model.evaluate(np.array([np.nan, 50.0]))
    assert np.isnan(state.saturation[0])
    following = model.evaluate(np.array([np.nan, 75.0]))
    assert np.isnan(following.saturation[0])
    assert following.saturation[1] == pytest.approx(0.5)  # keeps what it reached


def test_the_quality_mark_tells_a_measurement_from_a_bound(calibration: Calibration):
    """A node at a limit of the band is not a measurement, it is an 'at least' or 'at most'."""
    model = MoistureModel(calibration)
    state = model.evaluate(np.array([np.nan, 0.0, 50.0, 100.0]))
    assert state.quality.tolist() == [NO_DATA, WET_LIMIT, MEASURED, DRY_LIMIT]
    assert state.clamped == 2 and state.measured == 1


def test_invalid_threshold(calibration: Calibration):
    with pytest.raises(ValueError, match="threshold"):
        MoistureModel(calibration, saturation_threshold=0)
