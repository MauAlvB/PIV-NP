"""Calibration curve: interpolation, reading the file and clamping outside the range."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pivnp.moisture import Calibration, pchip_interpolate

CALIBRATION = Path(__file__).parent / "data" / "moisture" / "calibration_slope_rgb.csv"


# --- interpolation ------------------------------------------------------------------------
def test_pchip_reproduces_the_value_of_the_case():
    """The gray of the first step of the Slope_RGB case gives 0.149151, as in its data."""
    cal = Calibration.from_csv(CALIBRATION)
    result = cal.evaluate(np.array([6 / 11 * 100]))
    assert result.saturation[0] == pytest.approx(0.149151, abs=5e-7)
    assert result.clamped == 0


def test_pchip_goes_through_the_given_points():
    x = np.array([0.0, 1.0, 3.0, 6.0])
    y = np.array([0.0, 2.0, 2.5, 9.0])
    np.testing.assert_allclose(pchip_interpolate(x, y, x), y, atol=1e-12)


def test_pchip_is_exact_on_linear_data():
    x = np.array([0.0, 1.0, 2.0, 5.0])
    y = 3.0 * x - 1.0
    query = np.array([0.25, 1.5, 3.75, 4.99])
    np.testing.assert_allclose(pchip_interpolate(x, y, query), 3.0 * query - 1.0, atol=1e-12)


def test_pchip_with_two_points_is_the_straight_line():
    """That is what MATLAB does, and it is the smallest table a calibration admits."""
    x = np.array([10.0, 30.0])
    y = np.array([4.0, 0.0])
    query = np.array([10.0, 15.0, 20.0, 30.0])
    np.testing.assert_allclose(pchip_interpolate(x, y, query), [4.0, 3.0, 2.0, 0.0], atol=1e-12)


def test_pchip_does_not_overshoot_the_data():
    """Unlike a plain spline, it does not invent maxima between the points."""
    x = np.arange(6.0)
    y = np.array([0.0, 0.0, 0.0, 1.0, 1.0, 1.0])
    query = np.linspace(0, 5, 200)
    values = pchip_interpolate(x, y, query)
    assert values.min() >= -1e-12 and values.max() <= 1 + 1e-12
    assert np.all(np.diff(values) >= -1e-12)  # monotonically increasing


def test_pchip_rejects_invalid_input():
    with pytest.raises(ValueError, match="increasing order"):
        pchip_interpolate(np.array([1.0, 0.0]), np.array([0.0, 1.0]), np.array([0.5]))
    with pytest.raises(ValueError, match="same size"):
        pchip_interpolate(np.array([0.0, 1.0]), np.array([0.0]), np.array([0.5]))


# --- reading the file ---------------------------------------------------------------------
def test_reads_the_calibration_of_the_case():
    cal = Calibration.from_csv(CALIBRATION)
    assert cal.gray.size == 33
    assert cal.limits == (pytest.approx(0.00005), pytest.approx(100.0))
    # more gray (drier) means less saturation
    assert cal.saturation[0] > cal.saturation[-1]
    # the table is not strictly monotone: at the dry end the measured values rise a little
    # (0.0079, 0.0074, 0.0079). That is why pchip is useful: it does not invent a peak.
    assert np.all(np.diff(cal.saturation) <= 1e-3)


def test_accepts_comments_any_order_and_semicolons(workdir: Path):
    path = workdir / "cal.csv"
    path.write_text("# a note\ngray,saturation,moisture\n50;0.5;10\n\n0;1;25\n100;0;0\n")
    cal = Calibration.from_csv(path)
    assert cal.gray.tolist() == [0.0, 50.0, 100.0]
    assert cal.saturation.tolist() == [1.0, 0.5, 0.0]


@pytest.mark.parametrize(("contents", "message"), [
    ("gray,saturation,moisture\n", "no data rows"),
    ("gray,saturation\n0,1\n50,0.5\n", "3 columns"),
    ("gray,saturation,moisture\n0,1,25\n0,0.5,10\n", "no repeats"),
    ("gray,saturation,moisture\n0,1,25\n50,x,10\n", "non-numeric"),
    ("gray,saturation,moisture\n0,1,25\n", "at least 2 points"),
])
def test_rejects_malformed_files(workdir: Path, contents: str, message: str):
    path = workdir / "cal.csv"
    path.write_text(contents)
    with pytest.raises(ValueError, match=message):
        Calibration.from_csv(path)


# --- evaluation ---------------------------------------------------------------------------
@pytest.fixture
def simple() -> Calibration:
    return Calibration(np.array([0.0, 50.0, 100.0]), np.array([1.0, 0.5, 0.0]),
                       np.array([25.0, 10.0, 0.0]))


def test_clamps_outside_the_range_and_counts_it(simple: Calibration):
    result = simple.evaluate(np.array([-30.0, 0.0, 50.0, 100.0, 180.0]))
    assert result.clamped == 2
    assert result.saturation[0] == pytest.approx(1.0)   # clamped at the wet end
    assert result.saturation[-1] == pytest.approx(0.0)  # clamped at the dry end
    assert result.moisture.min() >= 0.0 and result.moisture.max() <= 25.0


def test_nodes_without_data_stay_without_data(simple: Calibration):
    result = simple.evaluate(np.array([np.nan, 50.0, np.nan]))
    assert np.isnan(result.saturation[[0, 2]]).all()
    assert result.saturation[1] == pytest.approx(0.5)
    assert result.clamped == 0


def test_moisture_never_leaves_the_range_of_the_table():
    """This is what used to fail: extrapolating gave water contents of −1700 %."""
    cal = Calibration.from_csv(CALIBRATION)
    rng = np.random.default_rng(0)
    gray = rng.uniform(-500, 500, size=10000)
    result = cal.evaluate(gray)
    assert result.moisture.min() >= cal.moisture.min() - 1e-9
    assert result.moisture.max() <= cal.moisture.max() + 1e-9
    assert 0.0 <= result.saturation.min() and result.saturation.max() <= 1.0
    assert result.clamped > 0
