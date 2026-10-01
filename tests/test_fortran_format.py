import math

import numpy as np
import pytest

from pivnp.fortran_format import format_e, format_float_rows, format_int_rows

# Real output of gfortran 5.3 with the E14.6 descriptor.
GFORTRAN_E14_6 = [
    (0.0, "  0.000000E+00"),
    (-0.0, " -0.000000E+00"),
    (1.0, "  0.100000E+01"),
    (0.1, "  0.100000E+00"),
    (-0.118643e-2, " -0.118643E-02"),
    (0.5078125, "  0.507812E+00"),  # an exact tie: round half to even
    (0.9999995, "  0.100000E+01"),
    (1234565.0, "  0.123456E+07"),  # an exact tie: round half to even
    (1e-100, "  0.100000E-99"),
    (-2.5e150, " -0.250000+151"),  # a 3-digit exponent: the letter E is dropped
    (9.999999e-10, "  0.100000E-08"),
    (1.0 / 3.0, "  0.333333E+00"),
    (math.nan, "           NaN"),
    (-math.inf, "     -Infinity"),
]


@pytest.mark.parametrize(("value", "expected"), GFORTRAN_E14_6)
def test_format_e_matches_gfortran(value, expected):
    assert format_e(value) == expected


@pytest.mark.parametrize(("value", "expected"), GFORTRAN_E14_6)
def test_numba_rows_match_gfortran(value, expected):
    row = format_float_rows(np.array([7]), np.array([[value]]), eol=b"\n")
    assert row.decode() == f"{7:14d}{expected}\n"


def test_numba_formatter_agrees_with_reference_on_random_values():
    rng = np.random.default_rng(0)
    values = np.concatenate([
        rng.normal(size=20000) * 10.0 ** rng.integers(-30, 30, size=20000),
        rng.uniform(-1, 1, size=5000),
        np.round(rng.uniform(-1, 1, size=5000), 7),  # plenty of near-ties
        10.0 ** np.arange(-20, 21),
        np.nextafter(10.0 ** np.arange(-20, 21), 0),
        [5e-324, 1.7976931348623157e308],
    ]).reshape(-1, 1)
    ids = np.arange(1, values.shape[0] + 1)
    text = format_float_rows(ids, values, eol=b"\n").decode()
    expected = "".join(f"{i:14d}{format_e(v)}\n" for i, v in zip(ids, values[:, 0], strict=True))
    assert text == expected


def test_float_rows_several_columns_and_crlf():
    out = format_float_rows(np.array([1, 22]), np.array([[1.0, -2.0], [0.5, 0.25]]),
                            eol=b"\r\n")
    assert out == (b"             1  0.100000E+01 -0.200000E+01\r\n"
                   b"            22  0.500000E+00  0.250000E+00\r\n")


def test_int_rows():
    out = format_int_rows(np.array([[1, 1, 2], [123, 123, -1]]), 9, eol=b"\n")
    assert out == b"        1        1        2\n      123      123       -1\n"


def test_empty_rows():
    assert format_float_rows(np.array([], dtype=int), np.zeros((0, 2))) == b""
