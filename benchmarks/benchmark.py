"""Mide el tiempo de un caso real con la versión Python (y opcionalmente con el Fortran).

Copia el caso a un directorio temporal para no sobrescribir resultados. Ejemplo::

    python benchmarks/benchmark.py "../caso" --threads 1 4 0 --legacy legacy/pivnp_legacy.exe

``--threads 0`` significa "todos los núcleos". La primera ejecución de Python incluye la
compilación JIT de Numba si la caché está vacía; se descarta con una ejecución de
calentamiento.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import tempfile
import time
from pathlib import Path

import numba

from pivnp.simulation import RunOptions, run_case


def _copy_inputs(case_dir: Path, dest: Path) -> None:
    for item in case_dir.iterdir():
        if item.is_file() and not item.name.upper().endswith((".RES", ".MSH", ".EXE")):
            shutil.copy(item, dest)


def time_python(case_dir: Path, threads: int) -> float:
    numba.set_num_threads(threads or numba.config.NUMBA_NUM_THREADS)
    with tempfile.TemporaryDirectory() as tmp:
        _copy_inputs(case_dir, Path(tmp))
        start = time.perf_counter()
        run_case(Path(tmp), options=RunOptions())
        return time.perf_counter() - start


def time_legacy(case_dir: Path, exe: Path) -> float:
    with tempfile.TemporaryDirectory() as tmp:
        _copy_inputs(case_dir, Path(tmp))
        start = time.perf_counter()
        subprocess.run([str(exe.resolve())], cwd=tmp, check=True, stdout=subprocess.DEVNULL)
        return time.perf_counter() - start


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("case_dir", type=Path)
    parser.add_argument("--threads", type=int, nargs="+", default=[0])
    parser.add_argument("--legacy", type=Path, help="ejecutable Fortran para comparar")
    parser.add_argument("--repeat", type=int, default=2)
    args = parser.parse_args()

    time_python(args.case_dir, 0)  # calentamiento / compilación JIT
    rows = []
    if args.legacy:
        rows.append(("Fortran original", min(time_legacy(args.case_dir, args.legacy)
                                             for _ in range(args.repeat))))
    for threads in args.threads:
        label = f"Python, {threads or numba.config.NUMBA_NUM_THREADS} hilos"
        rows.append((label, min(time_python(args.case_dir, threads)
                                for _ in range(args.repeat))))

    base = rows[0][1]
    print(f"{'Versión':<28}{'Tiempo (s)':>12}{'Aceleración':>14}")
    for label, seconds in rows:
        print(f"{label:<28}{seconds:>12.2f}{base / seconds:>13.1f}x")


if __name__ == "__main__":
    main()
