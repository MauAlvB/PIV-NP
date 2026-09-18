# PIV-NP: Particle Image Velocimetry con partículas numéricas

PIV-NP calcula **desplazamientos, deformaciones, velocidades, aceleraciones, energías y
humedad** de partículas numéricas que se mueven con el campo de velocidades que mide
[PIVlab](https://pivlab.blogspot.com/) sobre una secuencia de imágenes de un ensayo (por
ejemplo, un talud en centrífuga). A diferencia del PIV clásico, que da velocidades en una
malla fija (euleriana), PIV-NP sigue cada punto del material (lagrangiano), así que permite
**grandes desplazamientos** y deformaciones acumuladas.

> Pinyol, N.M. & Alvarado, M. (2017). *Novel PIV-based analysis for large displacement*.
> Canadian Geotechnical Journal 54(7): 933-944.

La versión 2.0 es una reingeniería en Python del código Fortran original
(`legacy/`). **Produce resultados idénticos byte a byte** a los del original y es unas
**19 veces más rápida** en el caso de la centrífuga (38 s → 2 s).

---

## Índice

1. [Instalación](#instalación)
2. [Uso](#uso)
3. [Archivos de entrada](#archivos-de-entrada)
4. [Resultados](#resultados)
5. [Cómo funciona](#cómo-funciona)
6. [Visualización](#visualización)
7. [Corrección de contorno (CONTOUR)](#corrección-de-contorno-contour)
8. [Rendimiento](#rendimiento)
9. [Pruebas](#pruebas)
10. [Estructura del proyecto](#estructura-del-proyecto)
11. [Hallazgos y propuestas](#hallazgos-y-propuestas)

## Instalación

Con [conda](https://docs.conda.io/) (recomendado):

```bash
conda env create -f environment.yml
```

```bash
conda activate pivnp
```

```bash
pip install -e ".[dev]"
```

Solo con pip (Python ≥ 3.10):

```bash
pip install -e ".[dev]"
```

La primera ejecución tarda unos segundos más porque Numba compila los núcleos de cálculo.
La compilación se guarda en caché y las ejecuciones siguientes arrancan al instante.

## Uso

### Línea de comandos

El directorio del caso debe contener `PIV-NP.TXT` (con el nombre del caso), `<caso>.PAR`
y los archivos de PIVlab, igual que con el ejecutable original:

```bash
pivnp ruta/al/caso
```

Opciones:

| Opción | Descripción |
|---|---|
| `--case NOMBRE` | nombre del caso (si no, se lee de `PIV-NP.TXT`) |
| `--threads N` | hilos de cálculo (por defecto, todos los núcleos) |
| `--prefetch N` | archivos PIVlab leídos por adelantado (por defecto 4) |
| `--contour-min-neighbors N` | con `ICONTOUR=1`: vecinos con datos necesarios (por defecto 3) |
| `--contour-layers N` | con `ICONTOUR=1`: capas de nodos a rellenar (por defecto 1) |
| `--vtk` | exportar también a VTK para ParaView (ver [Visualización](#visualización)) |
| `-q` | solo mostrar errores |

También funciona `python -m pivnp ruta/al/caso`.

### Desde Python

```python
from pivnp import Simulation

sim = Simulation.from_directory("ruta/al/caso")
summary = sim.run()

# Estado final en memoria, sin leer archivos:
sim.particles.displacement   # (n, 2) desplazamiento acumulado
sim.particles.eq_strain      # (n,)   deformación de corte equivalente
```

### Comparar resultados

```bash
python -m pivnp.compare referencia.POST.RES nuevo.POST.RES
```

Compara línea a línea y tolera diferencias de una unidad en la sexta cifra.

## Archivos de entrada

### `PIV-NP.TXT`

Una línea con el nombre del caso (sin extensión), por ejemplo `zapatak`.

### `<caso>.PAR`

Formato de lectura libre de Fortran: los valores pueden separarse por espacios,
tabuladores o comas. Las líneas de comentario son obligatorias.

```
Título del análisis
BLOQUE 2: NC NN NPC NFIL AXC AYC
2006  2100  3  34  0.212115  0.212115
BLOQUE 3: DT TOTAL_STEPS IMPPAS MOISTER IVERSION IPIVLAB ICONTOUR IREC ITR
0.8  149  1  0  1  1  0  0  0
BLOQUE 4: S_DENSITY POROSITY
2650.0  0.4
BLOQUE 5 (solo si ITR=1): DXT DYT PTVX1 PTVY1 PTVX2 PTVY2 PTVX3 PTVY3
```

| Parámetro | Significado |
|---|---|
| `NC`, `NN` | celdas y nodos (puntos) de la malla PIVlab; debe cumplirse NN = (NC/NFIL + 1)(NFIL + 1) |
| `NPC` | partículas por lado de celda (NPC² por celda) |
| `NFIL` | filas de celdas |
| `AXC`, `AYC` | ancho y alto de celda [m] |
| `DT` | tiempo entre imágenes [s] |
| `TOTAL_STEPS` | número de archivos PIVlab a procesar |
| `IMPPAS` | se escriben resultados en el paso 1 y cada `IMPPAS` pasos |
| `MOISTER` | 1 = leer también `Moist_<n>.TXT` (humedad y saturación) |
| `IVERSION` | 1 = malla PIV-NP igual a la de PIVlab; 2 = malla desplazada media celda |
| `IPIVLAB` | 1 = archivos de 4 columnas (PIVlab antiguo); otro = 5 columnas |
| `ICONTOUR` | 0 = sin corrección de contorno; 1 = corrección por media de vecinos (nuevo) |
| `IREC` | 0 = análisis nuevo; 1 = continuar desde `<caso>.REC` |
| `ITR` | 1 = añadir 3 partículas de seguimiento (PTV) definidas en el bloque 5 |
| `S_DENSITY`, `POROSITY` | se leen pero no intervienen en el cálculo |

A diferencia del original, se validan los datos (por ejemplo, que NC sea múltiplo de NFIL)
y los errores se explican con un mensaje en lugar de producir resultados sin sentido.

### Datos PIVlab

* `datos (1).txt`, `datos (2).txt`, …: exportación ASCII de PIVlab con 3 líneas de
  cabecera y columnas `x, y, u, v` (o 5 columnas con `IPIVLAB ≠ 1`). Los puntos sin
  medida llevan `NaN`.
* `Moist_1.TXT`, … (solo con `MOISTER=1`): 1 línea de cabecera y columnas
  `x, y, humedad, saturación`.

Los nombres no distinguen mayúsculas.

## Resultados

| Archivo | Contenido |
|---|---|
| `<caso>.POST.MSH` | malla de puntos para GiD (una por partícula). Material 1 = activa, 2 = sin datos en el paso 1, 3 = partícula PTV |
| `<caso>.POST.RES` | resultados por partícula en cada instante impreso (GiD) |
| `<caso>.REC` | estado final para continuar el análisis con `IREC=1` (binario compatible con el original) |

Resultados del `.POST.RES`:

| Nombre | Tipo | Descripción |
|---|---|---|
| `Displacement` | vector | desplazamiento acumulado |
| `Inst_displacement` | vector | desplazamiento del último paso |
| `NaNs` | escalar | contador heredado del original (ver H-11) |
| `Velocity`, `Acceleration` | vector | velocidad y aceleración de la partícula |
| `Total_strain`, `Inc_strain` | 3 comp. | εxx, εyy, γxy acumuladas / del paso |
| `Equi_strain`, `In_E_strain` | escalar | deformación de corte equivalente acumulada / del paso |
| `Vol_strain`, `Ins_vol_strain` | escalar | deformación volumétrica acumulada / del paso |
| `E_potential`, `E_kinetic`, `E_total` | escalar | energías por unidad de masa |
| `Moisture`, `Saturation` | escalar | solo con `MOISTER=1` |

## Cómo funciona

```
 imágenes ──PIVlab──► velocidades en una malla fija ──PIV-NP──► trayectorias y deformaciones
                      (datos (n).txt)                           de partículas del material
```

**1. Malla y partículas.** Los puntos de PIVlab forman una malla rectangular de
(NC/NFIL + 1) × (NFIL + 1) nodos. En cada celda se crean NPC × NPC partículas numéricas.
Con `IVERSION=2` se usa una malla desplazada media celda (cada punto PIVlab queda en el
centro de una celda) y la velocidad de cada punto se reparte a partes iguales entre las 4
esquinas de su celda.

**2. Velocidades nodales.** En cada paso *n* se lee `datos (n).txt`, se reordenan los
puntos (PIVlab los numera por columnas y de arriba a abajo; PIV-NP por filas y de abajo a
arriba) y se cambia el signo de *v* (el eje *y* de la imagen apunta hacia abajo).

**3. Interpolación a las partículas.** Cada partícula se localiza en su celda y se
calculan sus coordenadas locales (ξ, η) ∈ [−1, 1]. Con las funciones de forma bilineales

  Nᵢ(ξ, η) = ¼ (1 + ξ ξᵢ)(1 + η ηᵢ)

se obtienen la velocidad **v**ₚ = Σ Nᵢ **v**ᵢ, la aceleración
**a**ₚ = Σ Nᵢ (**v**ᵢⁿ − **v**ᵢⁿ⁻¹) / Δt y el desplazamiento del paso Δ**u**ₚ = **v**ₚ Δt.

**4. Deformaciones.** Con las derivadas de las funciones de forma en el centro del
elemento (la deformación es uniforme en cada celda):

  Δεxx = Σ ∂Nᵢ/∂x · vx,ᵢ Δt  Δεyy = Σ ∂Nᵢ/∂y · vy,ᵢ Δt  Δγxy = Σ (∂Nᵢ/∂y · vx,ᵢ + ∂Nᵢ/∂x · vy,ᵢ) Δt

Se acumulan y se calculan la deformación volumétrica εv = εxx + εyy y la de corte
equivalente εq = ⅔ q, con q = √(3 J₂) del tensor (εxx, εyy, εzz = 0, εxy = γxy/2).

**5. Actualización.** La partícula se mueve a **x**ₚ + Δ**u**ₚ. Si sale de la malla deja
de calcularse. Las partículas de celdas cuyos 4 nodos no tienen datos en el primer paso
se consideran fuera del material (aire, fondo) y se excluyen de los resultados.

**6. Resultados.** En el paso 1 y cada `IMPPAS` pasos se escriben todos los resultados;
al final se guarda el estado en `<caso>.REC`.

## Visualización

Recomendación: [ParaView](https://www.paraview.org/download/), gratuito y de código
abierto. PIV-NP convierte los resultados a VTK (su formato nativo), también los del
ejecutable Fortran original:

```bash
pivnp-vtk ruta/al/caso/zapatak.POST.RES
```

Después, en ParaView: **File → Open → `zapatak_vtk/zapatak.pvd` → Apply**, colorear por
`Equi_strain` y pulsar **Play**. Guía completa en
[`docs/VISUALIZACION.md`](docs/VISUALIZACION.md).

![Deformación de corte equivalente en el ensayo de centrífuga](docs/img/centrifuga_equi_strain.png)

## Corrección de contorno (CONTOUR)

En el borde del material PIVlab no da velocidad (NaN) en los puntos cuya ventana de
interrogación cae parcialmente fuera. Esos nodos entran en la interpolación con velocidad 0
(o la del paso anterior), y las partículas del borde se mueven menos de lo que deberían.

La propuesta implementada (`pivnp/contour.py`), que se activa con `ICONTOUR=1`:

1. A cada nodo sin datos con al menos `min_neighbors` vecinos con datos (de sus 8 vecinos)
   se le asigna la **media de sus velocidades**.
2. Se puede repetir `layers` veces para avanzar hacia el exterior.
3. Las marcas de "sin datos" originales no cambian, así que **no se activan partículas
   en el aire**: solo mejora la velocidad interpolada de las partículas del borde.

Se promedian nodos y no partículas porque la velocidad de las partículas se interpola desde
los nodos: promediar partículas sería circular y más caro.

Efecto en el caso de la centrífuga (149 pasos, 810 partículas de borde de 9378 activas):

| Configuración | Desplazamiento medio de las partículas de borde |
|---|---|
| Sin corrección (original) | 1.533 m |
| 3 vecinos, 1 capa | 1.621 m (+6 %) |
| 2 vecinos, 2 capas | 1.518 m (−1 %) |

Estos números son solo orientativos: conviene validarlos con marcadores PTV o con
fotografías del ensayo. Para añadir otro método basta con escribir una clase con el método
`apply(nodes, n_cols, n_rows)`.

> El original lee `ICONTOUR` pero lo ignora. Con `ICONTOUR=0` esta versión reproduce
> exactamente el original.

## Rendimiento

Caso de la centrífuga (18 054 partículas, 149 pasos, 768 MB de resultados), en un
procesador de 16 núcleos:

| Versión | Tiempo | Aceleración |
|---|---|---|
| Fortran original (gfortran -O2) | 38.0 s | 1× |
| Python, 1 hilo | 3.1 s | 12× |
| Python, 4 hilos | 2.3 s | 17× |
| Python, 16 hilos | 2.0 s | 19× |

Sin contar la primera compilación de Numba (unos 6 s, una sola vez). Para medirlo en tu
equipo:

```bash
python benchmarks/benchmark.py ruta/al/caso --threads 1 4 0 --legacy legacy/pivnp_legacy.exe
```

El detalle de dónde sale la mejora está en `docs/ARQUITECTURA.md`. La comparación es con
el Fortran compilado con gfortran; el ejecutable original de Intel Fortran puede ser algo
más rápido.

## Pruebas

```bash
pytest
```

* **Regresión** (`tests/test_regression.py`): 9 escenarios ejecutados con el Fortran
  original (malla 1 y 2, NPC = 2, 3, 4, 7 y 11, partículas PTV, humedad, formato de 5
  columnas y reinicio). Se exige que `.POST.RES`, `.POST.MSH` y `.REC` sean
  **idénticos byte a byte**.
* **Unitarias**: cada módulo se compara con una traducción literal de los bucles del
  Fortran (`tests/legacy_reference.py`) o con casos de solución conocida (traslación
  uniforme, campo de deformación lineal, formato E14.6 real de gfortran).

Cobertura (los núcleos de Numba solo se pueden medir sin compilar):

```bash
NUMBA_DISABLE_JIT=1 pytest --cov=pivnp
```

En PowerShell: `$env:NUMBA_DISABLE_JIT=1; pytest --cov=pivnp`. Cobertura actual: 96 %.

Para regenerar las referencias con el Fortran (requiere gfortran):

```bash
python tools/make_reference.py --source "ruta/a/caso paper centrifuga"
```

## Estructura del proyecto

```
piv-np/
├── src/pivnp/
│   ├── config.py          lectura y validación de PIV-NP.TXT y .PAR
│   ├── mesh.py            mallas y localización de partículas (UCELDA)
│   ├── particles.py       generación de partículas
│   ├── state.py           arrays de partículas y nodos (antes COMMON)
│   ├── pivlab_io.py       lectura de datos PIVlab y humedad
│   ├── nodal.py           campos nodales (VELOCIDADES)
│   ├── contour.py         corrección de contorno (CONTOUR)
│   ├── solver.py          movimiento y deformaciones (SOLMOV, INVAR2)
│   ├── gid_writer.py      salida GiD (IMPRES_GiD)
│   ├── fortran_format.py  formateo I/E14.6 rápido y exacto
│   ├── restart.py         archivo .REC (RECOM)
│   ├── simulation.py      bucle principal
│   ├── compare.py         comparación de resultados
│   ├── vtk_export.py      conversión a VTK para ParaView
│   └── cli.py             línea de comandos
├── tests/                 pruebas (unitarias y de regresión con datos reales recortados)
├── legacy/                código Fortran original, sin modificar
├── tools/                 generación de referencias con el Fortran
├── benchmarks/            medición de rendimiento
├── examples/              caso de la centrífuga (formato .PAR actual)
└── docs/
    ├── HALLAZGOS.md       errores y dudas encontrados en el original
    ├── ARQUITECTURA.md    diseño actual y propuestas de mejora
    └── VISUALIZACION.md   guía de ParaView
```

## Hallazgos y propuestas

Durante la migración se encontraron **22 puntos** en el código original. Los más
importantes:

* **H-01**: la aceleración usa una "velocidad anterior" equivocada en el 50 % de los nodos.
* **H-04**: con `IVERSION=2` las deformaciones salen divididas por NPC².
* **H-13**: un reinicio con `IVERSION=2` reparte todas las velocidades a una sola celda.
* **H-08**: la malla GiD se escribe con las posiciones del paso 1.

**No se ha corregido ninguno**: esta versión los reproduce para poder validar la migración.
La lista completa, con propuestas de corrección, está en
[`docs/HALLAZGOS.md`](docs/HALLAZGOS.md). Las propuestas de arquitectura (configuración
autodescriptiva, estrategias intercambiables, salida HDF5/VTK, CI, etc.) están en
[`docs/ARQUITECTURA.md`](docs/ARQUITECTURA.md).

## Licencia y cita

BSD de 4 cláusulas (ver `LICENSE`). Todo material que mencione el uso de este software debe
citar:

> Pinyol, N.M. & Alvarado, M. (2017). Novel PIV-based analysis for large displacement.
> Canadian Geotechnical Journal 54(7): 933-944.
