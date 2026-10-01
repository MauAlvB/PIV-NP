"""Command-line interface: ``pivnp <case_directory>``."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from .config import ConfigError
from .simulation import RunOptions, run_case, write_moisture_files
from .sources import available_sources


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="pivnp",
        description="Displacements and strains of numerical particles from PIVlab velocity "
                    "fields (PIV-NP).",
    )
    parser.add_argument("case_dir", type=Path, nargs="?", default=Path("."),
                        help="directory holding PIV-NP.TXT, <case>.PAR and the PIVlab data "
                             "(the current one by default)")
    parser.add_argument("--case", help="case name (otherwise it is read from PIV-NP.TXT)")
    parser.add_argument("--source", default="pivlab", choices=available_sources(),
                        help="where the displacement data comes from (default: pivlab)")
    parser.add_argument("--threads", type=int,
                        help="computation threads (all cores by default)")
    parser.add_argument("--prefetch", type=int, default=4,
                        help="steps read ahead (0 = no read-ahead)")
    parser.add_argument("--contour-min-neighbors", type=int, default=3,
                        help="with ICONTOUR=1: neighbours with data needed to rebuild a node")
    parser.add_argument("--contour-layers", type=int, default=1,
                        help="with ICONTOUR=1 and 3: layers of points to rebuild outwards")
    parser.add_argument("--contour-min-particles", type=int, default=1,
                        help="with ICONTOUR=2: particles needed around the point")
    parser.add_argument("--legacy-compat", action="store_true",
                        help="reproduce exactly the behaviour of the original Fortran, to "
                             "repeat analyses made with it")
    parser.add_argument("--legacy-2023-average", action="store_true",
                        help="with IVERSION=2, average the staggered-grid nodes the way the "
                             "2023 version did (to repeat those analyses)")
    parser.add_argument("--vtk", action="store_true",
                        help="when done, also export to VTK for ParaView (<case>_vtk/)")
    parser.add_argument("--write-moisture", action="store_true",
                        help="compute the moisture from the test images, write it to the "
                             "Moist_<n>.TXT files and exit without running the analysis")
    parser.add_argument("--convert-par", action="store_true",
                        help="move to the single format every .PAR in the directory and "
                             "below it, keeping each original as <case>.PAR.orig")
    parser.add_argument("-q", "--quiet", action="store_true", help="errors only")
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
        source=args.source,
        prefetch=args.prefetch,
        contour_min_neighbors=args.contour_min_neighbors,
        contour_layers=args.contour_layers,
        contour_min_particles=args.contour_min_particles,
        legacy_compat=args.legacy_compat,
        legacy_2023_average=args.legacy_2023_average,
    )
    log = logging.getLogger("pivnp")
    try:
        if args.convert_par:
            from .par_migrate import convert_tree

            results = convert_tree(args.case_dir)
            for result in results:
                (log.error if result.error else log.info)("%s", result)
            converted = sum(r.changed for r in results)
            failed = sum(r.error is not None for r in results)
            log.info("%d .PAR files: %d converted, %d already there, %d with problems",
                     len(results), converted, len(results) - converted - failed, failed)
            return 1 if failed else 0
        if args.write_moisture:
            write_moisture_files(args.case_dir, args.case, options)
            return 0
        summary = run_case(args.case_dir, args.case, options)
        if args.vtk:
            from .vtk_export import export_vtk

            res = Path(args.case_dir) / f"{summary.case_name}.POST.RES"
            log.info("Results for ParaView: %s",
                     export_vtk(res, legacy_mesh=args.legacy_compat))
    except (ConfigError, FileNotFoundError, ValueError) as exc:
        log.error("%s", exc)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
