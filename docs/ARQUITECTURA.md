# Arquitectura: estado actual y propuestas

## 1. Qué se ha hecho

El programa Fortran era un único archivo con 7 subrutinas que se comunicaban por bloques
`COMMON` globales (más de 60 arrays de tamaño fijo). La versión 2.0 separa
responsabilidades en módulos pequeños que se prueban por separado:

```
             ┌────────────┐   ┌─────────────┐
 .PAR ─────► │ config.py  │   │ pivlab_io.py│ ◄───── datos (n).txt, Moist_n.TXT
             └─────┬──────┘   └──────┬──────┘   (lectura anticipada en hilos)
                   │ CaseConfig      │ Frame
                   ▼                 ▼
             ┌──────────────────────────────────┐
             │ simulation.py  (bucle de pasos)  │
             └──┬──────────┬──────────┬──────┬──┘
                │          │          │      │
      mesh.py   │  nodal.py│ contour.py  solver.py        restart.py
      particles │  (nodos) │ (CONTOUR)   (partículas,     (.REC binario
      state.py  │          │             Numba paralelo)   compatible)
                ▼
          gid_writer.py + fortran_format.py  ──► .POST.MSH / .POST.RES
          (tabla de resultados, formateo paralelo, escritura en otro hilo)
```

| Fortran | Python |
|---|---|
| `COMMON /PARTICULAS/`, `/NODOS/` | `state.Particles`, `state.Nodes` (arrays con nombre) |
| `COMMON /GEOMETRIA/` + `UCELDA` | `mesh.Grid`, `mesh.locate_point` (O(1) en lugar de búsqueda lineal) |
| `PIVLAB_DATA` | `config.parse_par`, `particles.create_particles`, `simulation._load_restart` |
| `VELOCIDADES` | `pivlab_io.FrameSource` + `nodal.*` |
| `CONTOUR` | `contour.ContourCorrection` (interfaz) |
| `SOLMOV`, `INVAR2` | `solver.advance_particles`, `solver.update_strains`, `solver.deviatoric_q` |
| `IMPRES_GiD` (16 bucles copiados) | `gid_writer.RESULTS` (tabla) + `GidWriter` |
| `RECOM` | `restart.write_restart` |

De dónde sale la aceleración (38 s → 2 s en el caso de la centrífuga):

1. **Escritura de resultados** (lo que más pesaba): el original formatea número a número
   con `WRITE`. Ahora cada bloque se formatea en paralelo directamente en un buffer de bytes
   y el disco se escribe en un hilo aparte mientras se calcula el paso siguiente.
2. **Localización de partículas**: `UCELDA` recorría filas y columnas (O(filas+columnas))
   y se llamaba unas 20 veces por partícula y paso (una por cada bloque de resultados). Ahora
   es O(1) y se hace una vez por fase.
3. **Bucles por partícula** en paralelo con `numba.prange` (cada partícula es independiente).
4. **Lectura anticipada** de los archivos PIVlab en hilos.

Los bucles que acumulan sobre nodos (masa nodal, reparto en la malla desplazada) se dejan
secuenciales a propósito: sumar en otro orden cambiaría el redondeo y dejaría de coincidir
byte a byte con el original.

## 2. Propuestas para seguir mejorando

Ordenadas por relación beneficio/esfuerzo.

### P1. Modo de compatibilidad con el código original — **hecho**

`RunOptions(legacy_compat=True)` (`pivnp --legacy-compat`) reproduce exactamente el
comportamiento del Fortran, de modo que las pruebas de regresión siguen comparando byte a
byte con él y los análisis antiguos se pueden repetir con esta versión.

### P2. Archivo de configuración autodescriptivo (esfuerzo bajo)

El `.PAR` es posicional: un valor fuera de sitio cambia el significado de todo lo que
sigue, y ya hay tres variantes (la del código, la del tutorial y la del caso de ejemplo).
Propuesta: un `caso.toml` con nombres explícitos, unidades y valores por defecto,
manteniendo el lector de `.PAR` para los casos antiguos (`pivnp convert-par`).

```toml
[malla]
celdas = 2006
nodos = 2100
filas = 34
ancho_celda = 0.212115   # m
alto_celda = 0.212115    # m
particulas_por_lado = 3
version = 1

[analisis]
dt = 0.8                  # s
pasos = 149
imprimir_cada = 1

[contorno]
metodo = "media_vecinos"  # "ninguno" | "media_vecinos"
min_vecinos = 3
capas = 1
```

### P3. Estrategias intercambiables para las decisiones físicas (esfuerzo medio)

En el original, las alternativas se gestionaban comentando líneas (por ejemplo, la
deformación en el centro del elemento frente a la deformación en la posición de la
partícula, o masa `1.0` frente a `S_DENSITY*(1-POROSITY)*VVP`). Propuesta: convertir cada
decisión en una interfaz con implementaciones con nombre, elegidas desde la configuración,
igual que ya se ha hecho con `ContourCorrection`:

| Decisión | Opciones | Estado |
|---|---|---|
| Corrección de contorno | ninguna · media de vecinos · media de partículas · extrapolación | **hecho** (`ICONTOUR`) |
| Masa de partícula | unitaria (actual) · densidad seca × volumen (`S_DENSITY`, `POROSITY`) | pendiente |
| Cálculo de la deformación | centro del elemento (actual) · en la posición de la partícula | pendiente |
| Reparto inicial de partículas | Gauss (NPC 4-6) · uniforme (NPC 2-3) | fijado por NPC |

Así cada variante queda documentada, probada y seleccionable sin tocar el código, como ya
ocurre con `ContourCorrection`.

### P4. Formatos de salida modernos (esfuerzo medio, gran ganancia)

El `.POST.RES` de texto del caso de ejemplo ocupa **768 MB**. Propuesta: una interfaz
`ResultWriter` con varias implementaciones:

* `GidAsciiWriter` (la actual, compatibilidad).
* `Hdf5Writer` / `NetCdfWriter`: binario comprimido, típicamente de 10 a 20 veces más
  pequeño; se abre directamente desde Python, MATLAB o ParaView.
* `VtkWriter` (`.vtu` por paso + `.pvd`): visualización en ParaView sin licencia de GiD.
* `TrajectoryCsvWriter` con la trayectoria de puntos concretos, para comparar con
  marcadores de seguimiento PTV de laboratorio.

### P5. Entrada de datos desacoplada de PIVlab (esfuerzo bajo)

`FrameSource` ya devuelve un objeto `Frame` independiente del formato. Añadir lectores
para los `.mat` de PIVlab (evita exportar miles de `.txt`), OpenPIV o DaVis sería añadir
una clase, sin tocar el cálculo.

### P6. API para notebooks (esfuerzo bajo)

`Simulation` ya expone `particles` y `nodes` en memoria. Añadir un método que devuelva
los resultados como `xarray.Dataset` (dimensiones: tiempo × partícula) permitiría analizar
y graficar en Jupyter sin pasar por archivos.

### P7. Integración continua y publicación (esfuerzo bajo)

* GitHub/GitLab con CI que ejecute `pytest` y `ruff` en cada cambio (Windows y Linux).
* `pre-commit` con `ruff format`.
* Versionado semántico, `CHANGELOG.md` y publicación en PyPI o conda-forge.
* DOI de Zenodo para cada versión, citable junto al artículo de 2017.

### P8. Documentación (esfuerzo bajo)

Los módulos ya tienen docstrings. Con MkDocs + mkdocstrings se genera una web con la
referencia de la API, el tutorial de laboratorio (el `.docx` actual) y un notebook de
ejemplo con el caso de la centrífuga.

## 3. Hoja de ruta sugerida

1. **P2** + **P3**: configuración con nombres y opciones de cálculo seleccionables.
2. **P4** (HDF5/VTK) y **P7** (CI) para compartir el código con otros grupos.
3. **P5**, **P6** y **P8** según las necesidades del grupo.
