"""Water content and degree of saturation measured from the test images.

It replaces the set of MATLAB scripts that used to produce the ``Moist_<n>.TXT`` files, so
that the whole analysis can be run with a single program.
"""

from .calibration import Calibration, pchip_interpolate

__all__ = ["Calibration", "pchip_interpolate"]
