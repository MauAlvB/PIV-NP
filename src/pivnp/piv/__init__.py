"""Particle image velocimetry inside PIV-NP: a case can be analysed from the photographs.

PIV-NP has always taken its velocities from PIVlab. This package measures them instead, so
that someone with a sequence of images and no PIV experience has a way in, and can decide
later whether to learn a dedicated package and compare the two.

It is a basic PIV on purpose. ``correlation`` has the algorithm, ``settings`` reads the
``<case>.PIV`` file, and ``source`` is the :class:`pivnp.sources.DisplacementSource` that
``--source images`` selects.
"""

from .correlation import Field, one_pass, window_centres

__all__ = ["Field", "one_pass", "window_centres"]
