# Plan de validación caso por caso

Objetivo: reproducir con la versión nueva los resultados de cada caso analizado con versiones
anteriores, y corregir lo que vaya apareciendo. Se trabaja en la rama `humedad`.

Cada caso se ejecuta sobre una copia con enlaces duros en el directorio de trabajo temporal,
nunca sobre la carpeta original: los análisis escriben `.POST.MSH`, `.POST.RES` y `.REC`, y
los datos del ensayo no se tocan.

## Lo que hay

| Caso | Nodos | Pasos | `.PAR` | Resultado a reproducir |
|---|---|---|---|---|
| `centrifuga_humedad_caso_red` | 5432 | 149 | 2022, 8 valores + bloque 4 | `centrifuga.POST.RES` (48 MB) y `centrifuga1.POST.RES` (230 MB) |
| `centrifuga_humedad_caso_blue` | 5432 | 149 | igual | — (solo entradas) |
| `centrifuga_humedad_caso_green` | 5432 | 149 | igual | — |
| `centrifuga_humedad_editado` | 5800 | 149 | igual | — |
| `Digital…` etapa 1 (1fps) | 13806 | 20 | 9 valores sin bloque 4 | `DamFailure_prueba_1NP_v1.POST.RES` |
| `Digital…` etapa 1 - prueba | 13806 | 20 | 8 valores + bloque 4 | idem |
| `Digital…` etapa 1 - prueba 2 | 13806 | 20 | igual, IVERSION=2 | idem |
| `Digital…` etapa 2 (25fps) | 13806 | 20 | 9 valores, IREC=1 | idem |
| `Digital…` etapa 2 - prueba | 13806 | 20 | 8 valores + bloque 4, IREC=1 | idem |
| `Digital…` etapa 2 - prueba 2 | 13806 | 20 | igual, IVERSION=2 | idem |

Las dos etapas del caso del artículo son el mismo ensayo, seguidas: la 1 a 1 fps (imágenes 1
a 21) y la 2 a 25 fps (imágenes 21 a 41). Por eso la etapa 2 lleva `IREC=1`, continuando del
`.REC` de la etapa 1.

## Fase 0 — El lector de `.PAR` (bloquea todo lo demás)

**Hallazgo.** El bloque 3 no tiene un orden fijo: cada versión del programa colocó los
valores de una manera, y el mismo número de valores significa cosas distintas según el año.
Hasta ahora se adivinaba por el número de valores, y eso no basta:

| Dialecto | Bloque 3 | Bloque 4 |
|---|---|---|
| 2024 (actual) | `DT STEPS IMPPAS MOISTER IVERSION IPIVLAB ICONTOUR IREC ITR` | sí |
| Artículo, etapas | `DT STEPS IMPPAS MOISTER S_DENSITY POROSITY IVERSION PTV IREC` | no |
| Slope_RGB | `DT STEPS IMPPAS MOISTER S_DENSITY POROSITY IVERSION PTV` | no |
| Artículo, pruebas | `DT STEPS IMPPAS MOISTER IVERSION IPIVLAB IREC PTV` | sí |
| Centrífuga 2022 | `DT STEPS IMPPAS IVERSION IPIVLAB MOISTER IREC ITR` | sí |
| Antiguo | `DT STEPS IMPPAS [MOISTER]` | no |

Los dos últimos tienen los mismos ocho valores y el mismo bloque 4, y significan cosas
distintas: por posición son indistinguibles.

**Solución.** Leer la línea de comentario que precede a los valores, que es justamente donde
cada archivo documenta su propio orden (`del_t total_steps salto_impresión v.pivnp v.pivlab
moister REC PTR`). Se reconocen los nombres, con o sin acentos y en español o inglés, y se
asigna cada valor a su campo. Si la cabecera no se reconoce, se cae en la heurística actual
por número de valores. Así el lector deja de adivinar y pasa a leer lo que el archivo dice.

**Detalle a decidir:** en el `.PAR` de la centrífuga de 2022 aparece `MOISTER = 2`, y su
`.POST.RES` sí trae humedad y saturación, así que en aquella versión 2 seguía significando
"leer los archivos de humedad". En la versión nueva 2 significa "calcularla desde las
imágenes". Se resuelve por dialecto: en los `.PAR` antiguos, cualquier valor distinto de 0 se
entiende como "leer los archivos"; el 2 nuevo solo tiene su significado en el dialecto actual.

## Resultados de la fase 1

**Lo primero, y lo que da sentido a todo lo demás: la versión nueva es idéntica byte a byte
al Fortran de 2024 compilado.** Comprobado ejecutando los dos programas sobre tres casos
distintos —malla PIV-NP y malla desplazada, 20 y 149 pasos, de 5432 a 13806 nodos— y
comparando los quince resultados y la malla: ni una sola diferencia.

Con eso, cuando un resultado guardado no cuadra, la pregunta deja de ser "¿está bien nuestro
código?" y pasa a ser "¿con qué versión se hizo aquel análisis?".

| Caso | Resultado |
|---|---|
| Artículo, etapa 1 - prueba | **Exacto** |
| Artículo, etapa 1 - prueba 2 (malla desplazada) | **Exacto** con `--legacy-2023-average` |
| Artículo, etapa 2 - prueba (reinicio) | **Exacto** |
| Artículo, etapa 2 - prueba 2 (reinicio, malla desplazada) | Exacto en todos los valores; 4 partículas de diferencia en cuáles se imprimen |
| Artículo, etapa 1 (1fps) | Versión anterior: ver abajo |
| Centrífuga, caso rojo | Primer paso exacto con `--legacy-2023-average`; después divergen del 10 al 15 % de las trayectorias |

### Hallazgo: el reparto de la malla desplazada cambió entre versiones

Los resultados guardados con `IVERSION=2` no cuadraban. La relación entre sus valores y los
nuestros salía 4/3 y 3/2, que es lo que se obtiene al dividir por el número de puntos que
aportan a cada nodo. Se probó y **con esa media los quince resultados coinciden exactamente**,
con un detalle: aquella versión no aplicaba la media al incremento de cantidad de movimiento,
de modo que su aceleración no era coherente con su velocidad. Reproduciendo también eso, la
coincidencia es total.

En el Fortran de 2024 se ve el cambio, con la línea original comentada al lado:

```fortran
F1=0.25d0*1.0D0 !AMASSINI(I) ... !AM(JJ)      !Cuidado con la masa
AMASSNODE(JJ)=1.0D0 !AMASSNODE(JJ)+0.25d0*AMASSINI(I)
```

Es decir, la acumulación del peso nodal se sustituyó por la constante 1, y con ella
desapareció la media. **Ese es el origen del hallazgo H-23**, los nodos del borde con
velocidad infravalorada: no es un descuido de siempre, es una regresión introducida entre
versiones. Y confirma que promediar era lo que se hacía antes.

Se añade `--legacy-2023-average` para repetir aquellos análisis.

### Hallazgo: el formato de los archivos PIVlab no se puede fiar al `.PAR`

Los `.PAR` anteriores a 2024 no traen `IPIVLAB`, así que tomaba el valor 1, que significa
leer 4·NN valores seguidos sin mirar los saltos de línea. Los archivos del caso del artículo
tienen cinco columnas, con lo que a partir del primer valor todo se descoloca y el análisis
sale sin sentido. Ahora el formato se deduce del propio archivo. Con archivos de cuatro
columnas los dos caminos dan lo mismo, así que ningún análisis que ya funcionaba cambia.

### Sobre `etapa 1 (1fps)` y `etapa 2 (25fps)`

Son los análisis originales, hechos con una versión bastante anterior a la de las carpetas
`- prueba`. Comparando las partículas que ambos imprimen, **coinciden exactamente** los
desplazamientos, las velocidades, las energías, la humedad y la saturación. Difieren:

* aquella versión **imprimía las 13552 partículas**, también las que están fuera del material;
  la nuestra imprime solo las 5916 localizadas;
* las **deformaciones coinciden en el paso 1 y van divergiendo** después, lo que apunta a que
  el criterio para marcar una partícula como NaN era distinto, no a otra fórmula.

## Fase 1 — Reproducir los resultados mecánicos

Para cada caso con `.POST.RES` de referencia: ejecutar y comparar con `pivnp.compare`, que
tolera una unidad en la sexta cifra del formato `E14.6`. Se empieza por el más pequeño.

1. `Digital…` etapa 1 - prueba (20 pasos, IVERSION=1, IPIVLAB=2).
2. `Digital…` etapa 1 - prueba 2 (lo mismo con IVERSION=2: prueba la malla desplazada).
3. `Digital…` etapa 1 (1fps) (dialecto sin bloque 4).
4. `Digital…` etapa 2, las tres variantes (prueban el reinicio desde `.REC`).
5. `centrifuga_humedad_caso_red`, los dos resultados.
6. Los tres casos de centrífuga sin resultados: solo que corran y den algo coherente.

Cada diferencia que aparezca se analiza antes de tocar nada, se anota aquí y se corrige.

## Resultados de la fase 2: la humedad del caso del artículo

**Reproducido: el 99,8 % de los nodos sale con la humedad exacta**, y la saturación con una
diferencia media de 0,0024. Se llegó así:

1. **Qué tabla de calibración es.** Los pares (humedad, saturación) de los archivos del caso
   no caían sobre la tabla que usamos para `Slope_RGB`. Como los dos campos salen del mismo
   gris a través de la misma tabla, basta con mirar los pares: son los de
   `calibration_curve.m` de la carpeta `SWIR_SR`, no los de `calibration_curve_JC.m`, que
   tiene la saturación retocada a mano. La humedad coincide hasta el último decimal.

2. **Cuánto mide la banda.** En los 40 instantes solo aparecen **19 pares distintos**, y al
   despejar el gris normalizado de cada uno salen 0, 2,5, 5, 7,5 … 45: una escalera perfecta
   de 2,5. Como el gris de la imagen es entero, eso fija la anchura de la banda en
   **exactamente 40 niveles** (100/2,5). De paso confirma que aquel código sí redondeaba el
   filtro, que es lo que hace que los valores estén cuantizados.

3. **El registro entre las dos cámaras.** El MATLAB lo resolvía con una homografía de cuatro
   puntos marcados a mano, que no se guardó. Se recuperó de los propios datos: el gris
   medido en la imagen SWIR tiene que ser `saturado + gv_n/2,5`, así que se busca la
   homografía que menos dispersión deja en esa resta. Sale **0,294 niveles de gris**, por
   debajo de la cuantización, y desde tres puntos de partida distintos se llega a la misma.

4. **La banda, ya en números.** Saturado **92**, seco **132**.

El campo de humedad de este caso está bien calculado en todos los nodos (el fallo del bucle
que vacía los de `Slope_RGB` no está aquí) y **no lleva política incremental**: ningún nodo
tiene saturación exactamente 1 y la mínima baja y vuelve a subir a lo largo del ensayo. Es
decir, es el método del artículo publicado, y es el que conviene tomar como referencia.

### Lo que se añadió al código para poder reproducirlo

* **Registro por homografía.** `Registro` admite ahora los ocho valores de una homografía, no
  solo escala y origen, que es lo que hace falta cuando las dos cámaras miran desde ángulos
  distintos. En el `.HUM`: `INCLINACION_XY`, `INCLINACION_YX`, `PERSPECTIVA_X`,
  `PERSPECTIVA_Y`, además de las que ya había.
* **Banda global.** `BANDA_SECA` y `BANDA_SATURADA` fijan las dos intensidades para toda la
  imagen, en vez de sacarlas de una imagen de referencia nodo a nodo. Es lo que hace el flujo
  SWIR, y no cambia nada de lo anterior: si no se dan, se sigue usando la referencia por nodo.
* **Tabla de calibración** del suelo del artículo, en `tests/data/humedad/calibracion_swir.csv`.

El `.HUM` de ese caso queda así:

```
IMAGENES        = swir_{n}.jpg
CANAL           = gris
SIGMA           = 40
REDONDEO_LEGADO = 1          ! aquel análisis redondeaba el filtro
INCREMENTAL     = 0          ! el método del artículo no lleva trinquete
CALIBRACION     = calibracion_swir.csv
BANDA_SATURADA  = 92
BANDA_SECA      = 132
ESCALA_X        = 1.10390448    INCLINACION_XY = 0.09063166   ORIGEN_X = -92.27781
INCLINACION_YX  = 0.02232128    ESCALA_Y       = 1.09522258   ORIGEN_Y = -29.29620
PERSPECTIVA_X   = 1.3776e-05    PERSPECTIVA_Y  = 8.5586e-05
```

La diferencia que queda en la saturación (0,0024 de media) está en la tabla, no en la cadena:
su extremo húmedo tiene dos puntos de gris casi pegados (0 y 0,0001) con alturas muy
distintas, y ahí la interpolación es delicada. La humedad, que no pasa por ese punto, sale
exacta.

## Fase 2 — La humedad, contra el caso del artículo

Este caso es mejor banco de pruebas que `Slope_RGB`, y conviene decirlo claro:

* La humedad está calculada en **todos** los nodos con dato (7043 de 13806), no en una sola
  columna: el fallo del bucle que afecta a los datos de `Slope_RGB` no está aquí.
* **Ningún nodo tiene saturación exactamente 1**, y la saturación mínima baja y vuelve a subir
  a lo largo del ensayo (0,279 → 0,427 → 0,279). Es decir, **estos datos se calcularon sin la
  política incremental**: sin trinquete y sin umbral.

Es, por tanto, el método del artículo publicado, y es el que hay que reproducir. Pasos:

1. Ojear las 40 imágenes SWIR y las cuatro referencias, y averiguar cuál se usó y con qué σ.
2. Reconstruir la curva de calibración de este suelo.
3. Reproducir los `Moist_n.txt` nodo a nodo, como se hizo con `Slope_RGB`.
4. Con eso ya se puede comparar de verdad el efecto de las decisiones pendientes (trinquete,
   umbral, banda), porque aquí hay una referencia sin ellas.

Los `sc_gv_sat_swir_n.png` son figuras de MATLAB con barra de color, no datos: sirven para
mirar, no para medir. Los `Moist_n.txt` son la referencia numérica.

## Fase 3 — Mejoras del cálculo de humedad

Se prueban sobre el caso del artículo, midiendo cada una:

1. **Banda de calibración en intensidades absolutas** (lo que quedó aplazado). Aquí se puede
   evaluar de verdad, porque hay una referencia calculada sin trinquete.
2. **Medir en secado**: desactivar la política incremental. Los datos del artículo lo hacen,
   así que sirve para comprobar que nuestra versión no incremental es correcta.
3. **Filtro**: medir cuánto cambia el resultado con σ y si conviene un filtro que respete los
   bordes del material en vez de mezclar con el fondo.
4. **Registro entre cámaras**: el caso tiene imágenes visibles (1920×1080) y SWIR; si las dos
   mallas no coinciden, es el momento de probar el registro con escala y origen.

## Estado

| Fase | Estado |
|---|---|
| 0 | hecha: el lector de `.PAR` entiende los seis dialectos y el de PIVlab deduce el formato |
| 1 | hecha en lo esencial: cuatro casos reproducidos exactamente y los demás explicados |
| 2 | hecha: la humedad del caso del artículo se reproduce en el 99,8 % de los nodos |
| 3 | en marcha: ya están la banda global y el registro por homografía |
