[English](README.md) · **Español**

# PIV-NP: Particle Image Velocimetry con partículas numéricas

PIV-NP calcula **desplazamientos, deformaciones, velocidades, aceleraciones, energías y
humedad** de partículas numéricas que se mueven con el campo de velocidades que mide
[PIVlab](https://pivlab.blogspot.com/) sobre una secuencia de imágenes de un ensayo (por
ejemplo, un talud en centrífuga). A diferencia del PIV clásico, que da velocidades en una
malla fija (euleriana), PIV-NP sigue cada punto del material (lagrangiano), así que permite
**grandes desplazamientos** y deformaciones acumuladas.

> Pinyol, N.M. & Alvarado, M. (2017). *Novel PIV-based analysis for large displacement*.
> Canadian Geotechnical Journal 54(7): 933-944.

La versión 2.0 está escrita en Python y es unas **19 veces más rápida** que la versión
original en Fortran (`legacy/`): 38 s → 2 s en el caso de la centrífuga.

![Deformación de corte equivalente en el ensayo de centrífuga](docs/img/centrifuge_equi_strain.png)

---

## Índice

1. [Instalación](#instalación)
2. [Uso](#uso)
3. [Archivos de entrada](#archivos-de-entrada)
4. [Resultados](#resultados)
5. [Cómo funciona](#cómo-funciona)
6. [Visualización](#visualización)
7. [Corrección de contorno](#corrección-de-contorno)
8. [Rendimiento](#rendimiento)
9. [Pruebas](#pruebas)
10. [Estructura del proyecto](#estructura-del-proyecto)

## Instalación

Con [conda](https://docs.conda.io/) (recomendado):

```bash
conda env create -f environment.yml
```

```bash
conda activate pivnp
```

```bash
pip install -e .
```

Solo con pip (Python ≥ 3.10):

```bash
pip install -e .
```

La primera ejecución tarda unos segundos más porque Numba compila los núcleos de cálculo.
La compilación se guarda en caché y las siguientes arrancan al instante.

## Uso

### Línea de comandos

El directorio del caso debe contener `PIV-NP.TXT` (con el nombre del caso), `<caso>.PAR`
y los archivos exportados por PIVlab:

```bash
pivnp ruta/al/caso
```

Opciones:

| Opción | Descripción |
|---|---|
| `--case NOMBRE` | nombre del caso (si no, se lee de `PIV-NP.TXT`) |
| `--threads N` | hilos de cálculo (por defecto, todos los núcleos) |
| `--prefetch N` | archivos PIVlab leídos por adelantado (por defecto 4) |
| `--contour-min-neighbors N` | con `ICONTOUR=1`: vecinos con dato necesarios (por defecto 3) |
| `--contour-layers N` | con `ICONTOUR=1` y `3`: capas de puntos a reconstruir (por defecto 1) |
| `--contour-min-particles N` | con `ICONTOUR=2`: partículas necesarias alrededor (por defecto 1) |
| `--legacy-compat` | reproducir exactamente el comportamiento de la versión Fortran |
| `--legacy-2023-average` | reproducir la media en la malla desplazada de la versión de 2023 (solo `IVERSION=2`) |
| `--vtk` | exportar también a VTK para ParaView |
| `--write-moisture` | escribir los `Moist_<n>.TXT` desde las imágenes del ensayo y parar, sin analizar |
| `--convert-par` | reescribir todos los `.PAR` bajo el directorio al formato único y parar |
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

Formato de lectura libre de Fortran: los valores se separan por espacios, tabuladores o
comas. Las líneas de comentario son obligatorias y llevan los nombres de los campos encima
de sus valores, de modo que el archivo se explica solo.

```
Título del análisis
BLOCK 2: n_cells n_nodes n_part_cell n_rows width    height
         2006    2100    3           34     0.212115 0.212115
BLOCK 3: dt  total_steps print_every moisture mesh_version pivlab contour restart track
         0.8 149         1           0        1            1      0       0       0
BLOCK 4: s_density porosity
         2650.0    0.4
!------------------------------------------------------------------------------------------
! What every value means. Nothing below this line is read.
! ...
```

Después del bloque 4 va una leyenda (en inglés, como el resto del repositorio) con el
significado de cada número y las opciones de cada variable. No se lee: está para quien abra
el archivo. El conversor la escribe, y al crear un caso nuevo lo más cómodo es copiar el
`.PAR` de otro y cambiar los valores.

| Parámetro | Significado |
|---|---|
| `NC`, `NN` | celdas y nodos (puntos) de la malla PIVlab; debe cumplirse NN = (NC/NFIL + 1)(NFIL + 1) |
| `NPC` | filas y columnas de partículas por celda (NPC² por celda), de 1 a 6 |
| `NFIL` | filas de celdas |
| `AXC`, `AYC` | ancho y alto de celda [m] |
| `DT` | tiempo entre imágenes [s] |
| `TOTAL_STEPS` | número de archivos PIVlab a procesar |
| `IMPPAS` | se escriben resultados en el paso 1 y cada `IMPPAS` pasos |
| `MOISTER` | 0 = sin humedad; 1 = leerla de `Moist_<n>.TXT`; 2 = calcularla desde las imágenes del ensayo (`<caso>.HUM`) |
| `IVERSION` | 1 = las velocidades PIVlab están en los nodos de la malla; 2 = están en el centro de cada elemento (malla desplazada media celda) |
| `IPIVLAB` | 1 = archivos de 4 columnas (PIVlab antiguo); otro = 5 columnas |
| `ICONTOUR` | corrección de contorno: 0 = ninguna, 1 = media de vecinos, 2 = media de partículas, 3 = extrapolación |
| `IREC` | 0 = análisis nuevo; 1 = continuar desde `<caso>.REC` |
| `ITR` | debe ser 0 |
| `S_DENSITY`, `POROSITY` | se leen pero todavía no intervienen en el cálculo (la masa de partícula es 1) |

Los datos se validan al leerlos: si algo no encaja (por ejemplo, NC no múltiplo de NFIL) se
explica con un mensaje en lugar de producir resultados sin sentido.

#### Casos de versiones anteriores

El bloque 3 cambió de orden varias veces a lo largo de los años, y llegó a haber archivos con
los mismos valores colocados de otra manera. Hay un solo formato de entrada, así que un caso
antiguo se pasa a él una vez:

```bash
pivnp <directorio> --convert-par
```

Reescribe todos los `.PAR` que encuentre por debajo del directorio y guarda cada original al
lado, como `<caso>.PAR.orig`. La conversión reordena los valores tal cual, sin volver a
formatearlos, así que no cambia ni un decimal; los dos casos en los que sí corrige algo
—`IPIVLAB` según las columnas que de verdad tienen los archivos del caso, y el `MOISTER` de
las versiones en las que 2 significaba otra cosa— se avisan por pantalla.

### Datos PIVlab

* `datos (1).txt`, `datos (2).txt`, …: exportación ASCII de PIVlab con 3 líneas de
  cabecera y columnas `x, y, u, v` (o 5 columnas con `IPIVLAB ≠ 1`). Los puntos sin
  medida llevan `NaN`.
* `Moist_1.TXT`, … (solo con `MOISTER=1`): 1 línea de cabecera y columnas
  `x, y, humedad, saturación`.

Los nombres no distinguen mayúsculas.

### `<caso>.HUM` (solo con `MOISTER=2`)

Con `MOISTER=2` la humedad se mide a partir de las fotografías del ensayo en lugar de
leerse de los `Moist_<n>.TXT`. La configuración va en un archivo `<caso>.HUM`, donde cada
valor lleva su nombre, el orden es indiferente y `!` empieza un comentario:

```
! moisture measurement
IMAGES        = vis_{n}.jpg               ! {n} se sustituye por el número de paso
CHANNEL       = gray                      ! 1 rojo, 2 verde, 3 azul, 0 o "gray"
SIGMA         = 40                        ! radio de promediado, en píxeles
DRY_REFERENCE = ref2.jpg                  ! foto del suelo seco
CALIBRATION   = calibration_slope_rgb.csv ! gray,saturation,moisture del suelo
```

Las claves están en inglés, como el resto del repositorio.

| Clave | Por defecto | Significado |
|---|---|---|
| `IMAGES` | `vis_{n}.jpg` | nombre de la imagen de cada paso; debe contener `{n}` |
| `CHANNEL` | `gray` | canal usado: `1` rojo, `2` verde, `3` azul, `0`/`gray` escala de grises |
| `SIGMA` | `40` | radio del desenfoque gaussiano aplicado a cada imagen, en píxeles |
| `DRY_REFERENCE` | `ref2.jpg` | fotografía del suelo seco, la referencia de cada píxel |
| `DRY_OFFSET`, `SATURATED_OFFSET` | `5`, `-6` | niveles de gris por encima/debajo de la referencia que marcan suelo seco y saturado |
| `DRY_BAND`, `SATURATED_BAND` | — | las dos intensidades, fijas para toda la imagen, en lugar de la referencia por píxel. Las dos juntas o ninguna |
| `CALIBRATION` | — | **obligatoria**: CSV con la curva `gray,saturation,moisture` del suelo |
| `INCREMENTAL` | `0` | `1` = el suelo solo puede humedecerse: cada píxel conserva el valor más húmedo alcanzado |
| `SATURATION_THRESHOLD` | `0.95` | con `INCREMENTAL=1`, saturación por encima de la cual el píxel se fija en el extremo saturado |
| `FIRST_STEP` | `same` | `legacy` reproduce el código MATLAB, que daba humedad 0 y saturación constante al primer paso |
| `LEGACY_ROUNDING` | `0` | `1` redondea la imagen desenfocada a enteros, como hacía MATLAB |
| `SCALE_X`, `OFFSET_X`, `SCALE_Y`, `OFFSET_Y` | `1, 0, 1, 0` | dónde cae la malla PIV sobre la imagen de humedad |
| `SHEAR_XY`, `SHEAR_YX`, `PERSPECTIVE_X`, `PERSPECTIVE_Y` | `0` | el resto de términos de la transformación proyectiva, para una segunda cámara (SWIR) que mira desde otro ángulo |

El CSV de calibración lleva la cabecera `gray,saturation,moisture` y una fila por medida de
laboratorio; la curva se interpola entre ellas igual que el `pchip` de MATLAB. Los valores
fuera de la tabla se recortan a sus extremos y se informa de ello.

## Resultados

| Archivo | Contenido |
|---|---|
| `<caso>.POST.MSH` | malla de puntos para GiD (una por partícula), en sus posiciones iniciales. Material 1 = activa, 2 = sin datos en el paso 1 |
| `<caso>.POST.RES` | resultados por partícula en cada instante impreso (GiD) |
| `<caso>.REC` | estado final para continuar el análisis con `IREC=1` |

Resultados del `.POST.RES`:

| Nombre | Tipo | Descripción |
|---|---|---|
| `Displacement` | vector | desplazamiento acumulado |
| `Inst_displacement` | vector | desplazamiento del último paso |
| `NaNs` | escalar | datos que faltan alrededor: nodos sin medida del elemento (IVERSION=1) o punto central sin medida (IVERSION=2) |
| `Velocity`, `Acceleration` | vector | velocidad y aceleración de la partícula |
| `Total_strain`, `Inc_strain` | 3 comp. | εxx, εyy, γxy acumuladas / del paso |
| `Equi_strain`, `In_E_strain` | escalar | deformación de corte equivalente acumulada / del paso |
| `Vol_strain`, `Ins_vol_strain` | escalar | deformación volumétrica acumulada / del paso |
| `E_potential`, `E_kinetic`, `E_total` | escalar | energías por unidad de masa (`E_total` = potencial + cinética) |
| `Moisture`, `Saturation` | escalar | solo con `MOISTER=1` |

Como la malla lleva las posiciones iniciales, `malla + Displacement` es siempre la posición
actual de cada partícula.

### Continuar un análisis

Con `IREC=1` el análisis sigue donde lo dejó el `.REC` anterior: el tiempo y la numeración
de pasos continúan, y los resultados nuevos se **añaden** al `.POST.RES` existente. Los
archivos PIVlab de la continuación se numeran otra vez desde 1. Partir un análisis en dos
da exactamente el mismo resultado que hacerlo de una vez.

## Cómo funciona

```
 imágenes ──PIVlab──► velocidades en una malla fija ──PIV-NP──► trayectorias y deformaciones
                      (datos (n).txt)                           de partículas del material
```

**1. Malla y partículas.** Los puntos de PIVlab forman una malla rectangular de
(NC/NFIL + 1) × (NFIL + 1) nodos. En cada celda se crean NPC × NPC partículas numéricas, en
los puntos de Gauss (NPC de 4 a 6) o repartidas uniformemente (NPC 2 y 3). Con `IVERSION=2`
se usa una malla desplazada media celda, de modo que cada punto PIVlab queda en el centro de
un elemento y reparte su velocidad entre los 4 nodos de ese elemento.

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
abierto. PIV-NP convierte los resultados a VTK, su formato nativo, también los de la
versión Fortran:

```bash
pivnp-vtk ruta/al/caso/zapatak.POST.RES
```

Después, en ParaView: **File → Open → `zapatak_vtk/zapatak.pvd` → Apply**, colorear por
`Equi_strain` y pulsar **Play**. Guía completa en
[`docs/VISUALIZATION.md`](docs/VISUALIZATION.md).

## Corrección de contorno

En el borde del material, PIVlab no da velocidad (NaN) en los puntos cuya ventana de
interrogación cae parcialmente fuera. Esos puntos entran en la interpolación sin un valor
medido, de modo que las partículas del borde se mueven menos de lo que deberían y algunas
acaban despegándose del material.

`ICONTOUR` elige cómo reconstruir la velocidad de esos puntos:

| ICONTOUR | Método | Parámetros |
|---|---|---|
| 0 | ninguno | — |
| 1 | media de los nodos vecinos con dato (8 vecinos) | `--contour-min-neighbors`, `--contour-layers` |
| 2 | media de las velocidades de las partículas de los elementos de alrededor | `--contour-min-particles` |
| 3 | extrapolación lineal desde el interior hacia el exterior | `--contour-layers` |

Con cualquier método distinto de 0, la malla desplazada normaliza además el reparto de cada
punto entre los nodos de su celda, de modo que los nodos del borde reciben la media de los
puntos que aportan y no una fracción.

Los puntos reconstruidos se marcan aparte y **no** cambian las marcas de "sin dato", así que
la corrección no activa partículas en el aire: solo mejora la velocidad interpolada de las
partículas del borde.

Efecto en el caso de la centrífuga (149 pasos, 9378 partículas activas):

| Configuración | Partículas aisladas al final | Desplazamiento máximo |
|---|---|---|
| Sin corrección | 14 | 3.61 m |
| ICONTOUR=1, 3 vecinos, 1 capa | 0 | 3.67 m |
| ICONTOUR=2, media de partículas | 0 | 3.10 m |
| ICONTOUR=3, extrapolación | 0 | 3.39 m |

"Partículas aisladas" son las que acaban sin vecinas a menos de una celda, es decir, las que
se han despegado del material.

![Comparación de los métodos de corrección de contorno](docs/img/contour_comparison.png)

Cuál es el mejor es una cuestión física, no de programación: conviene contrastarlo con
marcadores PTV o con fotografías del ensayo. Para añadir otro método basta con escribir una
clase con un método `apply(ctx)` en `pivnp/contour.py`.

## Rendimiento

Caso de la centrífuga (18 054 partículas, 149 pasos, 768 MB de resultados), en un
procesador de 16 núcleos:

| Versión | Tiempo | Aceleración |
|---|---|---|
| Fortran (gfortran -O2) | 38.0 s | 1× |
| Python, 1 hilo | 3.1 s | 12× |
| Python, 4 hilos | 2.3 s | 17× |
| Python, 16 hilos | 2.0 s | 19× |

Sin contar la primera compilación de Numba (unos 6 s, una sola vez). Para medirlo en tu
equipo:

```bash
python benchmarks/benchmark.py ruta/al/caso --threads 1 4 0 --legacy legacy/pivnp_legacy.exe
```

El detalle de dónde sale la mejora está en [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md).
La comparación es con el Fortran compilado con gfortran.

## Pruebas

```bash
pytest
```

* **Casos sintéticos** (`tests/data/synthetic/`): campos de velocidad de solución conocida
  (desplazamiento sin deformación, deformación horizontal uniforme y variable, corte y
  rotación de sólido rígido). Se comprueba que las deformaciones calculadas coinciden con
  las analíticas y que se reproducen los análisis que el grupo hizo en su día.
* **Regresión**: 7 escenarios ejecutados también con la versión Fortran (las dos mallas,
  NPC = 2, 3 y 4, humedad, formato de 5 columnas y reinicio). En modo `--legacy-compat` se
  exige que `.POST.RES`, `.POST.MSH` y `.REC` sean **idénticos byte a byte**.
* **Comportamiento**: propiedades contrastadas con soluciones conocidas (aceleración
  constante, campo de deformación lineal, un análisis partido en dos equivalente a uno
  seguido).
* **Unitarias**: cada módulo se compara con una traducción literal de los bucles del
  Fortran (`tests/legacy_reference.py`) o con casos de solución conocida.

Cobertura (los núcleos de Numba solo se pueden medir sin compilar):

```bash
NUMBA_DISABLE_JIT=1 pytest --cov=pivnp
```

En PowerShell: `$env:NUMBA_DISABLE_JIT=1; pytest --cov=pivnp`. Cobertura actual: 96 %.

Para regenerar los resultados de referencia con la versión Fortran (requiere gfortran):

```bash
python tools/make_reference.py --source "ruta/a/caso paper centrifuga"
```

Además de las pruebas, [`docs/VALIDATION.md`](docs/VALIDATION.md) recoge lo comprobado
frente a los casos reales analizados con versiones anteriores: dónde esta versión los
reproduce exactamente, qué diferencias resultaron ser regresiones del código antiguo, cuánto
mueve cada decisión del cálculo de humedad y qué límites del método siguen abiertos.

## Estructura del proyecto

```
piv-np/
├── src/pivnp/
│   ├── config.py          lectura y validación de PIV-NP.TXT y .PAR
│   ├── par_migrate.py     conversión de los .PAR antiguos al formato único
│   ├── sources.py         de dónde vienen los desplazamientos (intercambiable)
│   ├── mesh.py            mallas y localización de partículas
│   ├── particles.py       generación de partículas
│   ├── state.py           arrays de partículas y nodos
│   ├── pivlab_io.py       lectura de datos PIVlab y humedad
│   ├── nodal.py           campos nodales
│   ├── contour.py         corrección de contorno
│   ├── solver.py          movimiento y deformaciones
│   ├── gid_writer.py      salida GiD
│   ├── fortran_format.py  formateo I/E14.6 rápido y exacto
│   ├── restart.py         archivo .REC
│   ├── simulation.py      bucle principal
│   ├── compare.py         comparación de resultados
│   ├── vtk_export.py      conversión a VTK para ParaView
│   ├── moisture/          humedad a partir de las imágenes del ensayo (MOISTER=2)
│   └── cli.py             línea de comandos
├── tests/                 pruebas (unitarias, de comportamiento y de regresión)
├── legacy/                código Fortran original, sin modificar
├── tools/                 generación de resultados de referencia
├── benchmarks/            medición de rendimiento
├── examples/              caso de la centrífuga
└── docs/
    ├── ARCHITECTURE.md    diseño y propuestas de mejora
    ├── VALIDATION.md      cómo se ha comprobado, y los límites conocidos
    └── VISUALIZATION.md   guía de ParaView
```

El código, sus comentarios y los documentos de `docs/` están en inglés, para que cualquiera
pueda usar y modificar el repositorio. Este archivo es la traducción al español del
[`README.md`](README.md).

## Licencia y cita

BSD de 4 cláusulas (ver [`LICENSE`](LICENSE)). Todo material que mencione el uso de este
software debe citar:

> Pinyol, N.M. & Alvarado, M. (2017). Novel PIV-based analysis for large displacement.
> Canadian Geotechnical Journal 54(7): 933-944.
