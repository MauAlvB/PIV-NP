"""Where the displacement data comes from.

PIV-NP moves its particles with a velocity field measured at the points of a grid. Today
that field is read from the files PIVlab exports, but nothing in the computation depends on
it: the analysis only ever sees a :class:`Frame` per step. This module writes that contract
down so that a different origin -- a PIV analysis built into PIV-NP, another PIV package, a
numerical simulation -- can be plugged in without touching the solver.

To add a source, implement the five members of :class:`DisplacementSource` and register a
factory for it::

    from pivnp.sources import Frame, register_source

    class MySource:
        def __init__(self, ...):
            self.moisture = False      # are there Moist_<n>.TXT files to read?
            self.images = None         # moisture computed from images, set by the analysis

        def frames(self, steps):
            for step in steps:
                u, v = my_velocities(step)          # m/s at every grid point, NaN if missing
                zeros = np.zeros(u.size)
                yield Frame(step, f"my data {step}", u, v, zeros, zeros)

        def mesh_in_metres(self, step=1):
            # x, y of every point in metres, metres per pixel, and which points have data
            return x, y, 0.004, np.isfinite(u)

        def frame_interval(self):
            return None                # seconds between images, or None if not known

    @register_source("mine")
    def _build(case_dir, config, prefetch):
        return MySource(...)

and select it with ``pivnp <case> --source mine``. The points must come in the order of the
PIVlab export (by columns, with y growing downwards): that is the order
``pivlab_to_node`` maps to PIV-NP nodes.
"""

from __future__ import annotations

from collections.abc import Callable, Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Protocol, runtime_checkable

import numpy as np

if TYPE_CHECKING:  # types only: the moisture package does not depend on this module
    from .config import CaseConfig
    from .moisture.source import MoistureSource


@dataclass(frozen=True)
class Frame:
    """Data of one step, in the order of the points of the PIVlab grid."""

    step: int
    #: Where this step came from: a file, or any short label. Only used for the log.
    source: Path | str | None
    u: np.ndarray  # x velocity (NaN where it was not measured)
    v: np.ndarray  # y velocity, image axis (downwards)
    moisture: np.ndarray
    saturation: np.ndarray
    #: Gray of the test image, already normalized, when the moisture is computed from the
    #: images. It is the intermediate step: the moisture comes from applying the model to it,
    #: and the model carries memory of the previous steps, so it cannot be computed here.
    normalized_gray: np.ndarray | None = None

    @property
    def label(self) -> str:
        """Short name of the origin of this step, for the log."""
        if isinstance(self.source, Path):
            return self.source.name
        return str(self.source) if self.source else f"step {self.step}"


@runtime_checkable
class DisplacementSource(Protocol):
    """What the analysis needs from whatever provides the velocity field.

    ``moisture`` and ``images`` are attributes, not methods: the analysis sets them when the
    ``.PAR`` asks for the moisture to be read (``MOISTER=1``) or computed from the test
    images (``MOISTER=2``). A source that cannot offer moisture just leaves them at
    ``False`` and ``None``.
    """

    #: Whether this source also provides moisture read from ``Moist_<n>.TXT``.
    moisture: bool
    #: Moisture computed from the test images, set by the analysis when ``MOISTER=2``.
    images: MoistureSource | None

    def frames(self, steps: range) -> Iterator[Frame]:
        """Yield the steps **in order**. The moisture model carries memory, so order matters."""
        ...

    def mesh_in_metres(self, step: int = 1) -> tuple[np.ndarray, np.ndarray, float, np.ndarray]:
        """Points in metres, metres per pixel, and which points have data.

        Only needed when the moisture is computed from the images, to place the grid on them.
        """
        ...

    def frame_interval(self) -> float | None:
        """Seconds between images this data was produced with, or ``None`` if unknown.

        The analysis compares it with the ``DT`` of the ``.PAR`` and warns when they differ,
        because the displacements would come out at a different scale.
        """
        ...


#: Factory of a source: receives the case and returns something satisfying the protocol.
SourceFactory = Callable[..., DisplacementSource]

#: Sources that can be selected by name with ``--source``.
SOURCES: dict[str, SourceFactory] = {}


def register_source(name: str) -> Callable[[SourceFactory], SourceFactory]:
    """Register a factory under ``name`` so that ``--source <name>`` finds it."""

    def decorator(factory: SourceFactory) -> SourceFactory:
        SOURCES[name] = factory
        return factory

    return decorator


def available_sources() -> list[str]:
    """Names that can be given to ``--source``."""
    _load_builtin_sources()
    return sorted(SOURCES)


def _load_builtin_sources() -> None:
    """Import the modules that register the sources shipped with PIV-NP."""
    from . import pivlab_io  # noqa: F401  registers "pivlab"
    from .piv import source  # noqa: F401  registers "images"


def build_source(name: str, case_dir: Path, config: CaseConfig,
                 prefetch: int = 4) -> DisplacementSource:
    """Build the source called ``name`` for a case."""
    _load_builtin_sources()
    try:
        factory = SOURCES[name]
    except KeyError:
        raise ValueError(f"unknown displacement source {name!r}; "
                         f"available: {', '.join(sorted(SOURCES))}") from None
    return factory(case_dir=Path(case_dir), config=config, prefetch=prefetch)
