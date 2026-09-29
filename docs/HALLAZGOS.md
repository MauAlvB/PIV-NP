# Hallazgos en el código original

Revisión de `legacy/MainCodePIV-NP.for` (v.2024.04.17) hecha durante la migración.

Gravedad: **Alta** = cambia resultados publicables · **Media** = afecta a casos concretos
o a la visualización · **Baja** = robustez, código muerto o estilo.

Estado: **Corregido** = la versión Python hace lo correcto · **Pendiente** = se reproduce el
comportamiento original, a la espera de decisión · **N/A en Python** = no existe en la
reescritura.

> **Modo compatibilidad.** Los hallazgos corregidos que cambian resultados (H-01, H-04,
> H-08 y H-13) se pueden revertir con `pivnp --legacy-compat` (o
> `RunOptions(legacy_compat=True)`), para repetir análisis antiguos. En ese modo las
> salidas siguen siendo idénticas byte a byte a las del Fortran, y así lo comprueban las
> pruebas de regresión.

## Resumen

| ID | Gravedad | Estado | Qué pasa |
|---|---|---|---|
| H-01 | Alta | **Corregido** | La "velocidad anterior" se guarda con el índice equivocado: la aceleración es incorrecta en el 50 % de los nodos |
| H-02 | Media | Pendiente | Un nodo que pasa a NaN conserva la velocidad del último paso válido |
| H-03 | Baja | N/A en Python | Bucle con `NP1` sin inicializar |
| H-04 | Alta | **Corregido** | Con IVERSION=2 las deformaciones salen divididas por NPC² |
| H-05 | Media | **Corregido** | Con NPC entre 7 y 10 todas las partículas se crean en el centro de la celda |
| H-06 | Baja | Pendiente | Reparto de partículas incoherente según NPC |
| H-07 | Baja | Pendiente | Literales en simple precisión (`9.81`, `1.E-10`, puntos de Gauss) |
| H-08 | Media | **Corregido** | La malla `.POST.MSH` se escribe tras el paso 1, no con las posiciones iniciales |
| H-09 | Baja | Pendiente | La marca `IDONDE` no se actualiza en el paso 1 |
| H-10 | Baja | Pendiente | Con IVERSION=2 el bucle de nodos usa la marca `IDONDE` de las partículas |
| H-11 | Media | Pendiente | El resultado "NaNs" imprime un contador nodal indexado por partícula |
| H-12 | Media | Pendiente | "E_kinetic" se declara `Scalar` pero escribe 2 valores |
| H-13 | Alta* | **Corregido** | Reinicio con IVERSION=2: todas las velocidades van a la celda 1 |
| H-14 | Media | Pendiente | Al reiniciar: tiempo a 0, "Inst_displacement" erróneo en el primer paso y se sobrescriben los resultados previos |
| H-15 | Media | **Eliminado** | Las partículas PTV (ITR) nunca aparecen en los resultados |
| H-16 | Baja | **Implementado** | `CONTOUR` no existe; `ICONTOUR` y `NODE_CONECT` sin uso |
| H-17 | Baja | N/A en Python | Escritura en `IELEMENT_ACTIVE(0)` (fuera de límites) |
| H-18 | Baja | N/A en Python | Código y variables muertos |
| H-19 | Media | **Corregido** | Tamaños fijos (200 000 partículas, 500 filas) sin comprobación |
| H-20 | Media | N/A en Python | Líneas de más de 72 columnas en formato fijo |
| H-21 | Baja | **Corregido** | Datos de entrada sin validar |
| H-22 | Baja | Pendiente | Umbral absoluto `J2 > 1e-10` |
| H-23 | Alta | Pendiente | Con IVERSION=2, los nodos del borde reciben una velocidad infravalorada |

\* Solo con reinicio (IREC=1) y malla desplazada (IVERSION=2).

## Detalle

### H-01 · Aceleración con la velocidad anterior equivocada — Alta

```fortran
DO I=1,NN
   VEL_X_NODO_V(I)=VEL_X_NODO(I)          ! I = índice del punto PIVlab
   ...
   VEL_X_NODO(ICONECTIVIDAD(I))=VEL_X(I)  ! se escribe en el nodo PIV-NP
```

`VEL_X_NODO` está indexado por **nodo PIV-NP**, pero la copia se hace con el índice del
**punto PIVlab** `I` dentro del mismo bucle que va sobrescribiendo `VEL_X_NODO`. Si el
nodo `I` ya se actualizó en una iteración anterior, la "velocidad anterior" es la
velocidad actual y `APV = m·(v − v_ant) = 0`.

* En el caso de la centrífuga (60×35 puntos) afecta a **1050 de 2100 nodos** en cada paso.
* Afecta al resultado **Acceleration**. No afecta a desplazamientos ni deformaciones.

**Corregido** en `nodal.load_measurements`: se copia el campo completo antes del bucle.
Los nodos afectados son los de la mitad superior de la malla, que en el ensayo de la
centrífuga es casi todo aire, así que allí el efecto es pequeño (el 0.6 % de las partículas
activas tenía la aceleración anulada y la media de |a| cambia un 1 %). En un ensayo donde
el material ocupe la parte superior de la imagen, el efecto sería mucho mayor.
Prueba: `tests/test_fixes.py::test_h01_acceleration_is_exact_for_every_node` (campo
uniforme con aceleración constante conocida).

### H-02 · Los nodos sin datos conservan la velocidad anterior — Media (¿intencionado?)

Con IVERSION=1, si `Node_NaN(I)=1`, `PV` y `APV` no se actualizan: conservan el valor del
último paso en que el nodo tuvo datos (o 0 si nunca lo tuvo). Las partículas de celdas con
algún nodo NaN siguen moviéndose con esa velocidad "congelada". En cambio, `VEL_X_NODO` sí
se pone a 0, así que el paso siguiente con datos calcula la aceleración respecto a 0.

**Pregunta:** ¿es intencionado (mantener el último valor conocido) o debería ser 0?
Relacionado con la corrección de contorno. Reproducido en `nodal.compute_nodal_momentum_v1`.

### H-03 · `NP1` sin inicializar — Baja

Con IREC=0 se ejecuta `DO I=1,NP1` y `NP1` es una variable local nunca asignada. Según el
compilador vale 0 (el bucle no hace nada) o basura (escritura fuera de límites o fallo).
Debería ser `NP`. Los arrays ya valen 0, así que el resultado no cambia. En Python no existe.

### H-04 · Deformaciones divididas por NPC² con IVERSION=2 — Alta

En `SOLMOV` las deformaciones usan `F10 = DX0*DT/AM(JJ)` y las velocidades usan `F1 = FN`
(la división por `AM` está comentada). Con IVERSION=1, `AM=1` y no pasa nada. Con
IVERSION=2, `AM` es la suma de las funciones de forma de las partículas (≈ NPC² en nodos
interiores), así que la deformación queda dividida por NPC² y es incoherente con el
desplazamiento.

Comprobado con un campo lineal `u = 0.01·x`, `dt = 0.5` (valor exacto 5·10⁻³):

| NPC | IVERSION=1 | IVERSION=2 |
|---|---|---|
| 2 | 5.000·10⁻³ | 1.250·10⁻³ (÷4) |
| 3 | 5.000·10⁻³ | 0.556·10⁻³ (÷9) |

**Corregido** en `solver._strain_kernel`: la deformación usa el mismo peso que la velocidad,
sin dividir por la masa nodal. Con IVERSION=1 (`AM=1`) los resultados no cambian.
Prueba: `tests/test_fixes.py::test_h04_strain_matches_analytic_field` con las dos mallas y
NPC = 2 y 3.

### H-05 · NPC entre 7 y 10 — Media

`GAUSS` solo se rellena para NPC=1…6, pero se usa si `NPC <= 10`. Para 7…10 el array
queda sin inicializar: con gfortran vale 0 y **todas** las partículas de la celda se crean
en el centro (con otro compilador, posiciones aleatorias).

**Corregido**: NPC queda limitado a 1…6 (`config.MAX_PARTICLES_PER_SIDE`) y un valor mayor
se rechaza con un mensaje claro. También se ha quitado el reparto uniforme para NPC > 10,
que ya no es alcanzable. Prueba:
`tests/test_fixes.py::test_h05_particles_per_side_limited_to_six`.

### H-06 · Reparto de partículas incoherente — Baja (diseño)

NPC=2 y 3 usan subceldas uniformes (±0.5; ±2/3, 0), NPC=4…6 usan puntos de Gauss-Legendre
(no uniformes, los valores uniformes están comentados) y NPC>10 vuelve a ser uniforme.
Todas las partículas tienen el mismo volumen `VVP`, lo que solo es coherente con el reparto
uniforme. Conviene decidir un único criterio.

### H-07 · Literales en simple precisión — Baja

`9.81` (SOLMOV) frente a `9.81d0` (inicialización), `1.E-10` (INVAR2) y los puntos de
Gauss se guardan como REAL*4: pierden precisión a partir de la 7.ª-8.ª cifra, y la energía
potencial inicial y la de cada paso usan una gravedad distinta (9.81 frente a
9.8100004196). El efecto es despreciable, pero si compilabas con `/real-size:64` los
resultados cambian ligeramente. En Python se reproducen en `constants.py`.

### H-08 · Malla GiD con las posiciones del paso 1 — Media

`IMPRES_GiD` escribe el `.POST.MSH` en el paso 1, **después** de `SOLMOV`, así que las
coordenadas ya incluyen el desplazamiento del primer paso. GiD dibuja la deformada como
malla + desplazamiento, con lo que ese primer incremento se cuenta dos veces.

**Corregido** en `Simulation._write_output`: la malla se escribe con `posición −
desplazamiento`, es decir, la posición inicial de cada partícula. Así
`malla + Displacement` es siempre la posición actual, también al reiniciar (donde el
desplazamiento sigue acumulándose desde el análisis anterior). El conversor a VTK usa esa
relación; para resultados antiguos hay que pasarle `--legacy-msh`. Prueba:
`tests/test_fixes.py::test_h08_mesh_plus_displacement_is_the_current_position`.

### H-09 · `IDONDE` no se actualiza en el paso 1 — Baja

`UCELDA` solo actualiza `IDONDE` si `IP /= 1`. Consecuencias: en el paso 1 una partícula
marcada como perdida no se recupera aunque esté dentro, y desde el paso 2 las partículas
NaN (`NaN_P=1`) vuelven a localizarse y se les calculan deformaciones (no se imprimen, pero
se guardan en el `.REC`). Reproducido en `particles.update_lost_flags`.

### H-10 · IVERSION=2 mezcla la marca de partícula con la de nodo — Baja

En VELOCIDADES, `CALL UCELDA(I)` con `IVERSIONCASO=2` localiza el centro de celda `XP2(I)`
pero lee y escribe `IDONDE(I)`, que es la marca de la **partícula** `I`. En el paso 1, si
la partícula `I` está fuera de la malla, la velocidad del nodo `I` no se reparte.
Reproducido en `nodal.compute_nodal_momentum_v2`.

### H-11 · Resultado "NaNs" sin sentido — Media

Imprime `ICOUNT_NO_NAN(I)` con `I` = número de partícula, pero `ICOUNT_NO_NAN` es un
contador por nodo de la malla desplazada (y con IVERSION=1 vale siempre 0). Además no se
filtran las partículas NaN. Si la idea era marcar las partículas sin datos, habría que
imprimir `NaN_P`/`NaN_P2`. Reproducido en `gid_writer._nodal_count_by_particle`.

### H-12 · "E_kinetic" declarado escalar con dos componentes — Media

La cabecera dice `Scalar` pero cada línea tiene dos valores (x, y). GiD lee solo el
primero. Corrección: declararlo `Vector` o escribir la suma. Ver `gid_writer.RESULTS`.

### H-13 · Reinicio con IVERSION=2 — Alta (solo en ese caso)

Los centros de celda `XP2` solo se generan si `IREC=0`. En un reinicio quedan a (0, 0) y
todos los puntos PIVlab reparten su velocidad a los nodos de la celda 1: el campo de
velocidades del resto de la malla es 0.

**Corregido** en `Simulation.__init__`: los centros se generan siempre. Pruebas:
`test_h13_restart_continues_the_analysis` comprueba, con las dos mallas, que 4 pasos
seguidos dan exactamente el mismo estado que 2 pasos + reinicio con los 2 siguientes;
`test_h13_legacy_restart_sends_everything_to_the_first_cell` documenta el comportamiento
antiguo.

### H-14 · Reinicio (IREC=1) — Media

* `TIEMPO` e `IP` vuelven a 0 y se leen otra vez `datos (1).TXT`… (si es intencionado,
  hay que usar otra carpeta con los datos de la continuación).
* En el primer paso `UPO = UP + incremento`: "Inst_displacement" muestra el desplazamiento
  acumulado en lugar del incremento.
* `.POST.RES` y `.POST.MSH` se sobrescriben: se pierden los resultados del análisis previo.
* Si el `.REC` es de otra malla, el original lo cargaba igualmente; la versión Python lo
  rechaza con un error.

### H-15 · Partículas PTV invisibles — Media

Las 3 partículas de seguimiento se marcan con `NaN_P=2`, y todos los resultados se filtran
con `NaN_P=0`: nunca aparecen en el `.POST.RES` ni se guardan en el `.REC`. Solo se ven en
la malla del paso 1 (material 3).

**Eliminado**: era un análisis puntual que nunca llegó a los resultados. `ITR` debe ser 0 y
un valor distinto se rechaza con un mensaje que explica el cambio. El material 3 ya no
existe en el `.POST.MSH`. Si en el futuro hiciera falta seguir puntos concretos, lo suyo
sería un archivo aparte con su trayectoria (ver `docs/ARQUITECTURA.md`, P4).

### H-16 · `CONTOUR` sin implementar — Baja

`CALL CONTOUR` no enlaza porque la subrutina no existe; `ICONTOUR` se lee pero no se usa, y
`NODE_CONECT` (vecinos de cada nodo) nunca se rellena, así que el bloque que cuenta NaN
vecinos (`ICOUNT_NAN`, `Node_NaN2`) nunca se ejecuta. Hay una propuesta implementada en
`contour.py` (ver README).

### H-17 · Escritura fuera de límites en `IELEMENT_ACTIVE` — Baja

Tras localizar cada partícula se escribe `IELEMENT_ACTIVE(INDC)` aunque la partícula no se
haya encontrado: con el `INDC` de la partícula anterior o con `INDC=0` (fuera del array,
pisa memoria del COMMON). El array no se lee nunca. En Python no existe.

### H-18 · Código y variables muertos — Baja

Nunca se leen: `SIG`, `VVP`, `MATP`, `AMASSINI`, `AMASSNODE`, `PVNAN`, `APVNAN`,
`SMOISTURE_NAN`, `SATURATION_NAN`, `Node_NaN_old`, `ICOUNT_NO_NAN_OLD`, `PV_OLD`,
`VEL_*_NODO_V2`, `Node_NaN_M`, `ICOUNT_NAN`, `Node_NaN2`, `IELEMENT_*`, `E_TOTAL_FULL`
(además divide por `count_active`, que puede ser 0), `NaN_P2` (solo alimenta variables
muertas), `POSCX/POSCY`, `S_DENSITY` y `POROSITY` (la masa de partícula está fijada a 1).
Con IVERSION=1, la masa nodal acumulada desde las partículas se sobrescribe con 1.

### H-19 · Tamaños fijos sin comprobación — Media

Los arrays tienen tamaño fijo (200 000 partículas o nodos, 500 filas) y no se comprueba si
se superan. Por ejemplo, el caso del tutorial (7670 celdas) con NPC=6 da
276 120 partículas y desbordaría la memoria sin avisar. La versión Python usa arrays
dinámicos.

### H-20 · Líneas de más de 72 columnas — Media

El archivo es de formato fijo, pero hay líneas que pasan de la columna 72 (por ejemplo la
de generación de partículas con `GAUSS(JJ)*AXC/2.`). Con las opciones por defecto,
gfortran las **trunca en silencio** y cambia los resultados. Solo es correcto compilando con
`-ffixed-line-length-none` (gfortran) o `/extend-source` (Intel). Ver `legacy/README.md`.

### H-21 · Entrada sin validar — Baja

No se comprueba que NC sea múltiplo de NFIL, que NN = (NC/NFIL+1)(NFIL+1), que IVERSION
sea 1 o 2, que IREC sea 0 o 1 ni que IMPPAS > 0 (con 0, `MOD(IP,0)` falla). La versión
Python valida todo esto y da un mensaje claro. Además, el `.PAR` de
`caso paper centrifuga` y el tutorial usan formatos distintos del que lee el código actual.

### H-22 · Umbral absoluto en INVAR2 — Baja

`q = 0` si `J2 <= 1e-10`. Es un umbral absoluto: incrementos de deformación del orden de
1e-5 darían `EPSEQ2 = 0`. En el caso de la centrífuga no afecta (0 % de partículas), pero
en ensayos lentos o con muchas imágenes podría. Un umbral relativo o mucho menor sería más
robusto.

### H-23 · IVERSION=2: velocidad infravalorada en los nodos del borde — Alta (pendiente)

Encontrado al corregir H-04. Con la malla desplazada, cada punto PIVlab reparte 1/4 de su
velocidad a los 4 nodos de su celda, **sumando sin normalizar**. Un nodo interior recibe
4 aportaciones (4 × ¼ = 1) y queda con la velocidad media correcta, pero un nodo del borde
exterior de la malla, o vecino de una zona sin datos, recibe solo 1, 2 o 3 aportaciones y
se queda con ¼, ½ o ¾ de la velocidad que le corresponde.

Esto ya pasaba en el original y afecta a **velocidades y desplazamientos**, no solo a las
deformaciones (allí quedaba parcialmente disimulado por la división por `AM` de H-04).
Medido en el primer instante del ensayo de la centrífuga con IVERSION=2:

* 188 de los 1140 nodos con datos (16 %) están infravalorados.
* 1520 de 9405 partículas activas (16 %) tienen una velocidad más de un 5 % baja.
* En esas partículas, la velocidad es de media el **74 %** de la correcta, y baja hasta el
  **19 %** en el peor caso.

Corrección propuesta: normalizar por el peso acumulado, es decir, velocidad nodal =
Σ wᵢ vᵢ / Σ wᵢ, con Σ wᵢ = 1 en el interior (no cambia nada allí) y el reparto correcto en
el borde. Es el mismo problema de fondo que resuelve la corrección de contorno para
IVERSION=1, así que conviene decidir los dos a la vez. **Pendiente de tu decisión**: cambia
desplazamientos y deformaciones de los análisis con IVERSION=2.
