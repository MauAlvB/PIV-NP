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
| 0 | por hacer |
| 1 | por hacer |
| 2 | por hacer |
| 3 | por hacer |
