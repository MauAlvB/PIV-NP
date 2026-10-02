"""Where the laboratory data the experiments need is found.

None of it is in the repository. The photographs are the raw data of published work, they
belong to the people who ran the tests, and they are large; the cases that go with them are
the same. So every path comes from outside the code, and an experiment that cannot find what
it needs says so in a sentence instead of unwinding a stack.

Set these before running the experiments that need them::

    $env:PIVNP_TEST_IMAGES = "C:\\path\\to\\Digital image-based measurement of"
    $env:PIVNP_DAM_BREAK   = "C:\\path\\to\\dam-break-swir"
    $env:SLOPE_RGB_CASE    = "C:\\path\\to\\Slope_RGB Completo"

``PIVNP_TEST_IMAGES`` is expected to hold ``PIVLAB/`` with the masked images PIVlab ran on,
and ``visual 1 fps (1 a 21)/`` with the camera originals. The other two are complete PIV-NP
cases, each with its ``.PAR`` and its data.

``PIVNP_DAM_BREAK`` may be left unset when the case sits at ``examples/dam-break-swir``,
which is where it lives for whoever has it.
"""
from __future__ import annotations

import os
from pathlib import Path

ENV = "PIVNP_TEST_IMAGES"
#: Cases that are not in the repository, and the variable that points at each.
DAM_BREAK = "PIVNP_DAM_BREAK"
SLOPE_RGB = "SLOPE_RGB_CASE"


def image_folder() -> Path:
    """The folder of test photographs, or a message saying how to point at it."""
    raw = os.environ.get(ENV)
    if not raw:
        raise SystemExit(
            f"These experiments need the test photographs, which are not in the repository.\n"
            f"Point at them with the {ENV} environment variable, for example:\n"
            f'    $env:{ENV} = "D:\\data\\Digital image-based measurement of"\n'
            f"The folder should contain PIVLAB/ and the 'visual 1 fps' images."
        )
    folder = Path(raw)
    if not folder.is_dir():
        raise SystemExit(f"{ENV} points at {folder}, which is not a directory.")
    return folder


def masked(step: int) -> Path:
    """One of the images PIVlab was run on, already masked and contrast-enhanced."""
    return image_folder() / "PIVLAB" / f"Masked_{step:03d}.jpg"


def original(name: str) -> Path:
    """One of the untouched camera images."""
    return image_folder() / "visual 1 fps (1 a 21)" / name


def case_folder(variable: str, what: str, fallback: Path | None = None) -> Path:
    """A complete case that is not in the repository, named by an environment variable.

    ``fallback`` is tried when the variable is unset, for a case that may already sit inside
    the working copy. There is deliberately no path baked in: one that happens to exist on
    the machine where the experiment was written makes it run there and fail everywhere
    else, which is worse than not running at all -- the failure at least says what is wrong.
    """
    raw = os.environ.get(variable)
    if raw:
        folder = Path(raw)
        if not folder.is_dir():
            raise SystemExit(f"{variable} points at {folder}, which is not a directory.")
        return folder
    if fallback is not None and fallback.is_dir():
        return fallback
    where = f"\nIt is not at {fallback} either." if fallback is not None else ""
    raise SystemExit(
        f"This experiment needs {what}, which is not in the repository: it is laboratory "
        f"data belonging to the people who ran the test.\n"
        f"Point at it with the {variable} environment variable, for example:\n"
        f'    $env:{variable} = "D:\\data\\the-case-folder"{where}'
    )
