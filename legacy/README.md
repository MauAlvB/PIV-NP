# Original Fortran code

* `MainCodePIV-NP.for`, `common_PIV-NP.for`: the original sources (v.2024.04.17),
  **unmodified**. They are the reference for the regression tests.
* `contour_stub.for`: an empty `CONTOUR` subroutine, needed to link: the original calls it
  but it was never written.

## Building

```bash
gfortran -O2 -finit-local-zero -ffixed-line-length-none -Wno-tabs -static -o pivnp_legacy.exe MainCodePIV-NP.for contour_stub.for
```

* `-ffixed-line-length-none` is **mandatory**: some lines run past column 72 and, without
  it, gfortran truncates them silently and the results change. With Intel Fortran the
  equivalent is `/extend-source`.
* `-finit-local-zero` gives a deterministic value to the uninitialised local variables
  (`NP1`, and `GAUSS` when NPC is between 7 and 10).

`tools/make_reference.py` builds the executable automatically when it is not there.
