# Código Fortran original

* `MainCodePIV-NP.for`, `common_PIV-NP.for`: fuentes originales (v.2024.04.17), **sin
  modificar**. Sirven de referencia para las pruebas de regresión.
* `contour_stub.for`: subrutina `CONTOUR` vacía, necesaria para enlazar (el original la
  llama pero no existe; ver `docs/HALLAZGOS.md`, H-16).

## Compilar

```bash
gfortran -O2 -finit-local-zero -ffixed-line-length-none -Wno-tabs -static -o pivnp_legacy.exe MainCodePIV-NP.for contour_stub.for
```

* `-ffixed-line-length-none` es **obligatoria**: hay líneas que pasan de la columna 72 y,
  sin ella, gfortran las trunca en silencio y cambia los resultados (H-20). Con Intel
  Fortran, el equivalente es `/extend-source`.
* `-finit-local-zero` da un valor determinista a las variables locales sin inicializar
  (`NP1`, `GAUSS` con NPC entre 7 y 10; ver H-03 y H-05).

`tools/make_reference.py` compila el ejecutable automáticamente si no existe.
