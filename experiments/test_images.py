"""Where the photographs the image experiments need are found.

Those experiments read the camera images of the dam-break test, which are not in the
repository: they are the raw data of a published paper and they are large. So the path comes
from outside the code.

Set it before running any of the image experiments::

    $env:PIVNP_TEST_IMAGES = "C:\\path\\to\\Digital image-based measurement of"

The folder is expected to hold ``PIVLAB/`` with the masked images PIVlab ran on, and
``visual 1 fps (1 a 21)/`` with the camera originals.
"""
from __future__ import annotations

import os
from pathlib import Path

ENV = "PIVNP_TEST_IMAGES"


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
