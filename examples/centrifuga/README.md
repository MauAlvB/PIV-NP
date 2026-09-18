# Ejemplo: talud en centrífuga

Datos del caso del artículo (`caso paper centrifuga`): 149 campos PIVlab de 60 × 35 puntos,
celdas de 0.212115 m y Δt = 0.8 s.

El `zapatak.PAR` original de esa carpeta está en el formato antiguo (el del `PIVNP.exe` de
2022, con 4 valores en el bloque 3). El de esta carpeta es el mismo caso en el formato
actual (9 valores en el bloque 3 y un bloque 4).

Los datos (`datos (1).txt` … `datos (149).txt`, 13 MB) no están en el repositorio. Para
ejecutarlo, cópialos aquí y lanza:

```bash
pivnp examples/centrifuga
```

Genera `zapatak.POST.MSH`, `zapatak.POST.RES` (≈ 770 MB) y `zapatak.REC`.
