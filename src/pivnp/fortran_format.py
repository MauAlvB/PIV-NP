"""Number formatting identical to the Fortran ``I<w>`` and ``E14.6`` descriptors.

Writing the results as text is the most expensive part of the program (hundreds of MB).
Here fixed-width rows are formatted straight into a byte buffer with Numba, in parallel. The
few values whose rounding to 6 digits is ambiguous in floating point (near ties) are rewritten
with :func:`format_e`, which uses Python's exact rounding.
"""

from __future__ import annotations

import math

import numpy as np
from numba import njit, prange

E_WIDTH = 14
E_DIGITS = 6
_SPACE, _MINUS, _PLUS, _DOT, _E, _ZERO = (ord(c) for c in " -+.E0")
_POW10 = np.array([10.0**k for k in range(23)])  # powers that are exact in binary
_AMBIGUITY = 1e-7


def format_e(value: float, width: int = E_WIDTH, digits: int = E_DIGITS) -> str:
    """Exact Python reference for ``Ew.d`` (gfortran's format)."""
    if math.isnan(value):
        return "NaN".rjust(width)
    if math.isinf(value):
        return ("-Infinity" if value < 0 else "Infinity").rjust(width)
    sign = "-" if math.copysign(1.0, value) < 0 else ""
    if value == 0:
        mantissa, exponent = "0" * digits, 0
    else:
        text = f"{abs(value):.{digits - 1}e}"  # d.ddddde±xx, correctly rounded
        mantissa = text[0] + text[2 : digits + 1]
        exponent = int(text[digits + 2 :]) + 1
    exp_text = f"E{exponent:+03d}" if abs(exponent) <= 99 else f"{exponent:+04d}"
    return f"{sign}0.{mantissa}{exp_text}".rjust(width)


@njit(cache=True, nogil=True, inline="always")
def _put_int(out, end, value):
    """Write ``value`` right-aligned, ending at ``end`` (exclusive)."""
    negative = value < 0
    v = -value if negative else value
    pos = end - 1
    while True:
        out[pos] = _ZERO + v % 10
        v //= 10
        pos -= 1
        if v == 0:
            break
    if negative:
        out[pos] = _MINUS


@njit(cache=True, nogil=True)
def _put_e14_6(out, start, value):
    """Write ``value`` in E14.6 format into ``out[start:start+14]``.

    Returns True when the rounding is ambiguous and has to be fixed with :func:`format_e`.
    """
    for k in range(E_WIDTH):
        out[start + k] = _SPACE
    if value != value or math.isinf(value):
        return True
    negative = math.copysign(1.0, value) < 0.0
    a = abs(value)
    digits = 0
    exponent = 0
    if a != 0.0:
        exponent = int(math.floor(math.log10(a))) + 1
        found = False
        for _ in range(3):
            k = E_DIGITS - exponent
            if k > 22 or k < -22:
                return True
            s = a * _POW10[k] if k >= 0 else a / _POW10[-k]  # a single rounded operation
            whole = math.floor(s)
            frac = s - whole
            if abs(frac - 0.5) < _AMBIGUITY:
                return True
            digits = int(whole) + (1 if frac > 0.5 else 0)
            if digits >= 1000000:  # rounding moves up a decade: 0.999999x -> 0.100000
                digits = 100000
                exponent += 1
            elif digits < 100000:  # log10 was off by one right at a power of 10
                exponent -= 1
                continue
            found = True
            break
        if not found or abs(exponent) > 99:
            return True

    pos = start + 1
    if negative:
        out[pos] = _MINUS
    out[pos + 1] = _ZERO
    out[pos + 2] = _DOT
    _put_int(out, pos + 3 + E_DIGITS, digits)
    for k in range(pos + 3, pos + 3 + E_DIGITS):  # leading zeros of the mantissa
        if out[k] == _SPACE:
            out[k] = _ZERO
    out[pos + 9] = _E
    out[pos + 10] = _MINUS if exponent < 0 else _PLUS
    e = -exponent if exponent < 0 else exponent
    out[pos + 11] = _ZERO + e // 10
    out[pos + 12] = _ZERO + e % 10
    return False


@njit(parallel=True, cache=True)
def _format_float_rows(ids, values, id_width, eol, out, ambiguous):
    n, k = values.shape
    row_width = id_width + E_WIDTH * k + eol.size
    for r in prange(n):
        base = r * row_width
        for c in range(id_width):
            out[base + c] = _SPACE
        _put_int(out, base + id_width, ids[r])
        for c in range(k):
            ambiguous[r, c] = _put_e14_6(out, base + id_width + E_WIDTH * c, values[r, c])
        for c in range(eol.size):
            out[base + row_width - eol.size + c] = eol[c]


@njit(parallel=True, cache=True)
def _format_int_rows(values, width, eol, out):
    n, k = values.shape
    row_width = width * k + eol.size
    for r in prange(n):
        base = r * row_width
        for c in range(width * k):
            out[base + c] = _SPACE
        for c in range(k):
            _put_int(out, base + width * (c + 1), values[r, c])
        for c in range(eol.size):
            out[base + row_width - eol.size + c] = eol[c]


def format_float_rows(ids: np.ndarray, values: np.ndarray, eol: bytes = b"\n",
                      id_width: int = 14) -> bytes:
    """Rows ``(I<id_width>, k*E14.6)``: one identifier and ``k`` reals per row."""
    values = np.ascontiguousarray(values, dtype=np.float64)
    if values.ndim == 1:
        values = values[:, None]
    ids = np.ascontiguousarray(ids, dtype=np.int64)
    n, k = values.shape
    row_width = id_width + E_WIDTH * k + len(eol)
    out = np.empty(n * row_width, dtype=np.uint8)
    ambiguous = np.zeros((n, k), dtype=np.bool_)
    _format_float_rows(ids, values, id_width, np.frombuffer(eol, np.uint8), out, ambiguous)
    for r, c in zip(*np.nonzero(ambiguous), strict=True):
        start = r * row_width + id_width + E_WIDTH * c
        out[start : start + E_WIDTH] = np.frombuffer(format_e(values[r, c]).encode(), np.uint8)
    return out.tobytes()


def format_int_rows(values: np.ndarray, width: int, eol: bytes = b"\n") -> bytes:
    """Rows of integers, right-aligned in fields of ``width`` characters."""
    values = np.ascontiguousarray(values, dtype=np.int64)
    n, k = values.shape
    out = np.empty(n * (width * k + len(eol)), dtype=np.uint8)
    _format_int_rows(values, width, np.frombuffer(eol, np.uint8), out)
    return out.tobytes()
