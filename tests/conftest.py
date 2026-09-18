from __future__ import annotations

import tempfile
from pathlib import Path

import pytest


@pytest.fixture
def workdir():
    """Directorio temporal que se borra al terminar la prueba.

    Se usa en lugar de ``tmp_path`` porque pytest crea enlaces simbólicos en su directorio
    temporal y en algunos Windows sin permiso de symlinks cada intento tarda ~30 s.
    """
    with tempfile.TemporaryDirectory(prefix="pivnp-test-") as directory:
        yield Path(directory)
