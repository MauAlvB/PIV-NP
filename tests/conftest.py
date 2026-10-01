from __future__ import annotations

import tempfile
from pathlib import Path

import pytest


@pytest.fixture
def workdir():
    """Temporary directory removed when the test finishes.

    It is used instead of ``tmp_path`` because pytest creates symbolic links in its own
    temporary directory, and on some Windows installs without permission to create symlinks
    each attempt takes about 30 s.
    """
    with tempfile.TemporaryDirectory(prefix="pivnp-test-") as directory:
        yield Path(directory)
