"""Medición de humedad y grado de saturación a partir de las imágenes del ensayo.

Reemplaza al conjunto de scripts MATLAB que hasta ahora generaba los archivos
``Moist_<n>.TXT``, de modo que el análisis completo se pueda hacer con un solo programa.
"""

from .calibracion import Calibracion, interpolar_pchip

__all__ = ["Calibracion", "interpolar_pchip"]
