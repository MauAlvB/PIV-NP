# Rama `humedad`: estado del trabajo

Documento de trabajo de esta rama. Recoge lo decidido, lo averiguado y lo que falta, para
poder retomar el hilo en cualquier momento. `main` queda intacta mientras tanto.

## Objetivo

Calcular la humedad y el grado de saturación desde las imágenes del ensayo dentro de PIV-NP,
sin pasar por el código MATLAB ni por archivos intermedios, de modo que un solo programa
haga todo el análisis.

## Base de partida

De las cuatro colecciones de `codigos humedad`, la que generó los datos del caso
`Slope_RGB Completo` es la de `gv_sat = gv − 6` (`Code_Matlab` / `CODE_def_mod_otro`).

Comprobado reproduciendo su saturación inicial: en el primer instante el gris normalizado
vale 6/(5+6)·100 = 54.5455 % y, pasado por la curva de calibración con interpolación pchip,
da **0.149151**, que es exactamente el único valor que aparece en `Moist_1.txt`. Con
`gv_sat = gv − 8` daría 0.102742.

Se toma como base la **lista de archivos de `Code_Matlab`** (dice qué 18 scripts importan),
los **números de `CODE_def_mod_otro`** y los **lectores con canal de `CODE_def_mod``**.

## Lo que se averiguó sobre el caso de referencia

Verificado por huella SHA-256 de los archivos:

* Las imágenes de humedad **son las mismas que usó PIVlab**: `vis_1.jpg` es idéntica byte a
  byte a `1X-V1 001.jpg`. Por tanto **el registro espacial es la identidad** en este caso y
  no hacen falta los cuatro puntos que se pinchaban a mano.
* Las 150 `swir_*.jpg` son **copias exactas** de las 150 `vis_*.jpg`: en esa carpeta no hay
  imágenes infrarrojas reales.
* `ref1.jpg` y `ref2.jpg` son ambas copias de la **primera imagen del ensayo**. Por eso en el
  primer instante la imagen coincide con la referencia seca.

Consecuencia práctica: el caso sirve para validar la cadena completa de forma exacta.

## Algoritmo a portar

Por cada instante:

1. Leer la imagen y quedarse con un canal (1 rojo, 2 verde, 3 azul) o convertirla a gris.
2. Filtro gaussiano de radio σ (40 en el flujo RGB, 15 en el SWIR).
3. Muestrear el gris en cada nodo de la malla PIV.
4. Normalizar: `gv_n = (gv − gv_sat) / (gv_dry − gv_sat) · 100`, recortado a 0 por abajo,
   con `gv_dry = gv_ref + 5` y `gv_sat = gv_ref − 6` nodo a nodo.
5. Interpolar en la curva de calibración (pchip) para obtener saturación y humedad.
6. Política incremental: la saturación no baja y se fija en 1 al superar el umbral.

Detalle a respetar para reproducir: MATLAB usa núcleo de tamaño `2·ceil(2σ)+1` y relleno por
repetición del borde; las coordenadas de los nodos se redondean al píxel.

## Cambios acordados respecto al MATLAB

| # | Cambio | Estado |
|---|---|---|
| 1 | Umbral de saturación de 0.8 a **0.95**, manteniendo el comportamiento incremental | aplazado: por defecto sigue 0.8 hasta cerrar la comparación; se pide con `UMBRAL_SATURACION` |
| 2 | El **primer instante se calcula como los demás** (hoy devuelve humedad 0 y saturación constante) | acordado |
| 3 | Referencias seca y saturada como **parámetros**, con estudio de sensibilidad (−6 frente a −8 mueve la saturación inicial un 45 %) | acordado, valor por decidir con pruebas |
| 4 | **Recortar** en vez de extrapolar fuera del rango de calibración (origen de los −1710 % de `Moist_40.txt`) | hecho |
| 5 | El filtro **deja de redondear** a niveles enteros de gris; `REDONDEO_LEGADO = 1` recupera el comportamiento antiguo | hecho |

La forma de fijar la banda seca-saturada (los desplazamientos +5 y −6 frente a una banda de
laboratorio en intensidades absolutas) queda **como está** de momento, a la espera de tenerla
mejor definida.

Más adelante: permitir medir también en secado, quitando el carácter incremental.

## Arquitectura prevista

```
src/pivnp/humedad/
├── calibracion.py   tabla gris-saturación-humedad, pchip, control de rango
├── imagenes.py      lectura, selección de canal, filtro gaussiano
├── muestreo.py      coordenadas de los nodos en píxeles, muestreo y máscara
├── referencias.py   referencias seca y saturada
├── modelo.py        gris normalizado -> saturación y humedad, política incremental
└── fuente.py        orquestación; entrega los campos de cada instante
```

Encaja en lo existente a través de `Frame(u, v, humedad, saturación)`: una fuente alternativa
calcula esos dos campos desde las imágenes y el resto del programa no cambia.

Datos de entrada nuevos: un archivo de configuración del ensayo (canal, σ, referencias,
registro, umbral) y un archivo de calibración por suelo, que viene de ensayo de laboratorio.
El `.PAR` no se toca.

## Plan por fases

| Fase | Contenido | Validación |
|---|---|---|
| 0 | Comparador de archivos de humedad y datos de referencia congelados | mide diferencias nodo a nodo |
| 1 | Calibración y modelo (pchip, normalización, política incremental) | reproduce las saturaciones del caso |
| 2 | Imágenes, filtro, muestreo | los `Moist_n.txt` generados coinciden con los del caso |
| 3 | Integración en el análisis y comando propio | el análisis desde imágenes = el análisis leyendo los archivos |
| 4 | Los cuatro cambios acordados | pruebas propias y cuantificación del efecto |
| 5 | Opcional: flujo SWIR con dos cámaras | pendiente de disponer de imágenes SWIR reales |

Como en la migración del Fortran, primero se reproduce el comportamiento actual y solo
después se activan los cambios, para poder separar lo que es traducción de lo que es mejora.

## Resultados de la primera validación contra el caso

Con la cadena completa implementada (imagen → canal → filtro → muestreo → normalización →
calibración → política incremental) y comparando con los `Moist_n.txt` del caso:

Reproduciendo el camino de MATLAB (canal gris, `REDONDEO_LEGADO = 1`):

| Instante | Nodos comparables | Saturación idéntica | Difieren | Diferencia media |
|---|---|---|---|---|
| 1 | 1045 | 1045 (error 3·10⁻¹⁶) | 0 | 0 |
| 2 | 1045 | 930 | 11 % | 0,025 |
| 5 | 1045 | 797 | 24 % | 0,054 |
| 10 | 1023 | 704 | 31 % | 0,085 |
| 149 | 784 | 610 | 22 % | 0,090 |

Entre el 69 % y el 100 % de los nodos coinciden **exactamente**. Los que difieren lo hacen
en saltos de un nivel de gris.

### Hallazgos importantes

**1. En los datos del caso la humedad solo se calculó en una columna de nodos.** En
`gv2sat_JC2.m` el cálculo de la humedad (líneas 39-40) está **fuera del bucle interior**:
se ejecuta una vez por fila, con el último valor de la columna. Verificado en los archivos
del caso: solo 34 de 2100 nodos (1,6 %) tienen humedad distinta de cero, y son siempre los
últimos del archivo, es decir la última columna de la malla. El resto se queda con el cero
con el que se inicializa la matriz.

Consecuencia: **el campo de humedad de esos análisis está esencialmente vacío**. La
saturación sí se calcula bien en todos los nodos, porque ese cálculo está dentro de los dos
bucles. Nuestra implementación calcula la humedad en todos los nodos, que es lo que el
código pretendía hacer.

**2. El método es muy sensible al redondeo del filtro.** La banda entre el gris seco y el
saturado son solo **11 niveles de gris** (+5 y −6), así que **un nivel de diferencia equivale
al 9 % de toda la escala** y puede cambiar la saturación de un nodo hasta en 0,85.

Se midió el tamaño del residuo invirtiendo la curva de calibración: de la saturación de los
archivos del caso se deduce qué gris tenía MATLAB y se compara con el nuestro. Sale
**exactamente un nivel de gris** (mediana 0,99, p90 1,02). No hay diferencia de algoritmo.

Se descartó que `imgaussfilt` filtrase en el dominio de la frecuencia: filtrando por FFT con
el mismo relleno se obtiene la misma imagen que por convolución directa, con una diferencia
máxima de 4·10⁻¹³, y los mismos nodos exactos. La causa que queda es una diferencia
sistemática de menos de un nivel entre la descodificación del JPEG de MATLAB y la de Pillow,
que el redondeo del filtro convierte en un nivel entero. Encaja con que el instante 1
coincida al 100 %: allí la imagen y la referencia son la misma, y cualquier sesgo se cancela
al normalizar.

**Decisión tomada:** quitar el redondeo. No empeora el parecido con los datos antiguos —la
diferencia media de saturación pasa de 0,0249 a 0,0220 en el instante 2 y de 0,0850 a 0,0886
en el 10—, simplemente deja de caer en los mismos valores discretos, y elimina de raíz la
fragilidad.

**3. El caso se procesó en escala de grises, no con un canal.** El instante 1 no lo
distingue, porque la imagen es la propia referencia. Comparando los instantes siguientes con
los cuatro canales, el gris es el que reproduce los datos (diferencia media 0,025 en el
instante 2, frente a 0,034 del verde, 0,089 del azul y 0,113 del rojo).

**4. La humedad y la saturación pueden contradecirse.** La política incremental se aplica
solo a la saturación; la humedad se recalcula entera en cada instante. Como las imágenes de
este ensayo se aclaran con el tiempo (la mediana del gris de un nodo sube 3,8 niveles del
instante 1 al 149), los dos campos acaban diciendo cosas opuestas:

| Paso | Nodos con saturación = 1 | Humedad de esos nodos (mín / mediana / máx) | De ellos, bajan de humedad |
|---|---|---|---|
| 10 | 107 | 10,81 / 18,47 / 24,03 | 54 |
| 40 | 569 | 0,20 / 3,13 / 24,03 | 404 |
| 149 | 634 | 0,20 / **1,09** / 24,03 | 190 |

En la tabla de calibración la saturación 1 corresponde al 24,03 % de humedad. En el MATLAB no
se veía porque la humedad estaba casi toda a cero por el fallo del bucle. **Se deja como está
a propósito**, para que la comparación con los análisis anteriores sea limpia; es lo primero a
revisar cuando esa comparación se cierre.

Efecto del umbral, medido en el caso: con 0,95 en vez de 0,8 quedan al final 491 nodos fijados
en saturación 1 en vez de 634, con una diferencia media de 0,022 y máxima de 0,20.

## Integración en el análisis (fase 3)

`MOISTER`, en el bloque 3 del `.PAR`, pasa a admitir tres valores: 0 sin humedad, 1 leerla de
los `Moist_<n>.TXT` (lo de siempre) y 2 calcularla desde las imágenes con la configuración de
`<caso>.HUM`. El `.PAR` no cambia de formato.

La lectura y el filtrado de cada imagen van con la lectura anticipada de los archivos PIVlab,
en varios hilos. El modelo lleva memoria —la saturación no baja—, así que se aplica después,
con los instantes ya en orden: el instante transporta el gris normalizado y la humedad se
calcula al entregarlo. Con `pivnp <caso> --write-moisture` se generan los `Moist_<n>.TXT` sin
ejecutar el análisis, para revisarlos o compararlos con los de análisis anteriores.

### Comprobación de que la integración no cambia nada más

Se analizó el caso completo (149 pasos, 18054 partículas) por los dos caminos, con la misma
configuración salvo el origen de la humedad, y se compararon los `.POST.RES` resultado a
resultado:

| Resultado | Filas | Distintas | Diferencia máxima |
|---|---|---|---|
| Desplazamientos, velocidades, deformaciones, energías, NaNs (14 resultados) | 54507 cada uno | **0** | 0 |
| Moisture | 54507 | 44404 | 24,0 |
| Saturation | 54507 | 21892 | 0,85 |

La malla también sale idéntica. Los dos únicos resultados que cambian son los que tienen que
cambiar: la humedad, porque en los archivos del caso está casi toda a cero por el fallo del
MATLAB, y la saturación, en la medida ya conocida del nivel de gris.

## Estado

* Rama `humedad` con la cadena completa implementada y 267 pruebas en verde.
* Fases 1 (calibración), 2 (imágenes, muestreo, modelo) y 3 (integración y comando)
  terminadas, probadas y validadas contra el caso; el residuo está explicado y acotado.
* **Comparación cerrada** y conclusiones aplicadas (ver `validacion-casos.md`):
  * La política incremental actúa sobre el **gris normalizado**, no sobre la saturación, así
    que los dos campos salen del mismo valor y no pueden contradecirse. Comprobado sobre
    `Slope_RGB`: la saturación sale **exactamente igual** que con la política anterior
    (diferencia máxima 0), y los 579 nodos dados por saturados pasan de tener humedad mediana
    1,14 % a tener la del suelo saturado, 24,03 %.
  * El trinquete viene **apagado** (`INCREMENTAL = 0`) y su umbral en 0,95: es una hipótesis
    sobre el ensayo, no una medida, y con él puesto no se puede medir el secado.
  * Se **avisa cuando la banda es estrecha** (menos de 20 niveles de gris, donde un nivel pasa
    a pesar más del 5 % de la escala). En `Slope_RGB`, con 11 niveles, salta.
  * Cada nodo lleva una **marca de calidad**: medido, sin dato, o en un tope de la banda, donde
    el valor es una cota y no una medida. Al terminar el análisis se dice cuánto es cada cosa.
* Sigue pendiente: cómo se mide la banda seca-saturada de cada suelo, que es el número más
  sensible de todo el método.
* Más adelante: medir también en secado (sin política incremental), el flujo SWIR con dos
  cámaras, y la masa de las partículas desde `S_DENSITY` y `POROSITY`.
* Nada subido a GitHub; `main` tiene un commit local por delante del remoto.
