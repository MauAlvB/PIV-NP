"""Interfaz de línea de comandos: ``pivnp <directorio_del_caso>``."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from .config import ConfigError
from .simulation import RunOptions, run_case


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="pivnp",
        description="Desplazamientos y deformaciones de partículas numéricas a partir de "
                    "campos de velocidad de PIVlab (PIV-NP).",
    )
    parser.add_argument("case_dir", type=Path, nargs="?", default=Path("."),
                        help="directorio con PIV-NP.TXT, <caso>.PAR y los datos PIVlab "
                             "(por defecto, el actual)")
    parser.add_argument("--case", help="nombre del caso (si no, se lee de PIV-NP.TXT)")
    parser.add_argument("--threads", type=int,
                        help="hilos de cálculo (por defecto, todos los núcleos)")
    parser.add_argument("--prefetch", type=int, default=4,
                        help="archivos PIVlab leídos por adelantado (0 = sin lectura anticipada)")
    parser.add_argument("--contour-min-neighbors", type=int, default=3,
                        help="con ICONTOUR=1: vecinos con dato necesarios para reconstruir")
    parser.add_argument("--contour-layers", type=int, default=1,
                        help="con ICONTOUR=1 y 3: capas de puntos a reconstruir hacia fuera")
    parser.add_argument("--contour-min-particles", type=int, default=1,
                        help="con ICONTOUR=2: partículas necesarias alrededor del punto")
    parser.add_argument("--legacy-compat", action="store_true",
                        help="reproducir exactamente el comportamiento del Fortran "
                             "original, para repetir análisis hechos con él")
    parser.add_argument("--vtk", action="store_true",
                        help="al terminar, exportar también a VTK para ParaView (<caso>_vtk/)")
    parser.add_argument("-q", "--quiet", action="store_true", help="solo errores")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.WARNING if args.quiet else logging.INFO,
        format="%(asctime)s, %(message)s",
        datefmt="%H:%M:%S",
    )
    if args.threads:
        import numba

        numba.set_num_threads(args.threads)

    options = RunOptions(
        prefetch=args.prefetch,
        contour_min_neighbors=args.contour_min_neighbors,
        contour_layers=args.contour_layers,
        contour_min_particles=args.contour_min_particles,
        legacy_compat=args.legacy_compat,
    )
    log = logging.getLogger("pivnp")
    try:
        summary = run_case(args.case_dir, args.case, options)
        if args.vtk:
            from .vtk_export import export_vtk

            res = Path(args.case_dir) / f"{summary.case_name}.POST.RES"
            log.info("Resultados para ParaView: %s",
                     export_vtk(res, legacy_mesh=args.legacy_compat))
    except (ConfigError, FileNotFoundError, ValueError) as exc:
        log.error("%s", exc)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
