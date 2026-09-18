"""Comparación de archivos de resultados GiD (``.POST.RES`` / ``.POST.MSH``).

Uso: ``python -m pivnp.compare referencia.POST.RES nuevo.POST.RES``

Las líneas se comparan primero como texto. Si difieren, se comparan token a token y los
números se aceptan si difieren como mucho en ``ulps`` unidades de la sexta cifra (el
formato E14.6 tiene 6 cifras significativas).
"""

from __future__ import annotations

import argparse
import itertools
import sys
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class ComparisonReport:
    lines: int = 0
    identical_lines: int = 0
    tolerated_lines: int = 0
    mismatched_lines: int = 0
    max_abs_diff: float = 0.0
    examples: list[tuple[int, str, str]] = field(default_factory=list)

    @property
    def equivalent(self) -> bool:
        return self.mismatched_lines == 0

    def summary(self) -> str:
        return (f"{self.lines} líneas: {self.identical_lines} idénticas, "
                f"{self.tolerated_lines} con diferencias de redondeo, "
                f"{self.mismatched_lines} distintas (máx. dif. abs. {self.max_abs_diff:.3e})")


def _last_digit_unit(token: str) -> float:
    """Valor de una unidad en la última cifra de un número en formato E (0.dddddd E±xx)."""
    mantissa, _, exponent = token.upper().partition("E")
    decimals = len(mantissa.split(".")[1]) if "." in mantissa else 0
    return 10.0 ** (int(exponent or 0) - decimals)


def _tokens_match(a: str, b: str, ulps: int, report: ComparisonReport) -> bool:
    if a == b:
        return True
    try:
        x, y = float(a), float(b)
    except ValueError:
        return False
    diff = abs(x - y)
    report.max_abs_diff = max(report.max_abs_diff, diff)
    if "E" not in a.upper():
        return diff == 0
    return diff <= ulps * _last_digit_unit(a) * (1 + 1e-9)


def compare_files(reference: Path, candidate: Path, ulps: int = 1,
                  max_examples: int = 10) -> ComparisonReport:
    report = ComparisonReport()
    with open(reference, encoding="ascii", errors="replace") as fa, \
            open(candidate, encoding="ascii", errors="replace") as fb:
        for number, (la, lb) in enumerate(itertools.zip_longest(fa, fb), start=1):
            report.lines += 1
            la, lb = (la or "").rstrip("\r\n"), (lb or "").rstrip("\r\n")
            if la == lb:
                report.identical_lines += 1
                continue
            ta, tb = la.split(), lb.split()
            if len(ta) == len(tb) and all(_tokens_match(a, b, ulps, report)
                                          for a, b in zip(ta, tb, strict=True)):
                report.tolerated_lines += 1
                continue
            report.mismatched_lines += 1
            if len(report.examples) < max_examples:
                report.examples.append((number, la, lb))
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("reference", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--ulps", type=int, default=1,
                        help="tolerancia en unidades de la última cifra (por defecto 1)")
    args = parser.parse_args(argv)
    report = compare_files(args.reference, args.candidate, args.ulps)
    print(report.summary())
    for number, a, b in report.examples:
        print(f"  línea {number}:\n    ref: {a}\n    new: {b}")
    return 0 if report.equivalent else 1


if __name__ == "__main__":
    sys.exit(main())
