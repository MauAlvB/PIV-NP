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
| 1 | Umbral de saturación de 0.8 a **0.95**, manteniendo el comportamiento incremental | acordado |
| 2 | El **primer instante se calcula como los demás** (hoy devuelve humedad 0 y saturación constante) | acordado |
| 3 | Referencias seca y saturada como **parámetros**, con estudio de sensibilidad (−6 frente a −8 mueve la saturación inicial un 45 %) | acordado, valor por decidir con pruebas |
| 4 | **Recortar** en vez de extrapolar fuera del rango de calibración (origen de los −1710 % de `Moist_40.txt`) | acordado, forma por confirmar |

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

| Instante | Nodos comparables | Saturación idéntica | Difieren |
|---|---|---|---|
| 1 | 1045 | 1045 (error 3·10⁻¹⁶) | 0 |
| 2 | 1045 | 930 | 11 % |
| 5 | 1045 | 797 | 24 % |
| 10 | 1023 | 704 | 31 % |
| 149 | 784 | 610 | 22 % |

Entre el 69 % y el 100 % de los nodos coinciden **exactamente**. Los que difieren lo hacen
en saltos de un nivel de gris.

### Dos hallazgos importantes

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
al 9 % de toda la escala** y puede cambiar la saturación de un nodo hasta en 0,85. Las
diferencias que quedan tienen exactamente ese tamaño: un nivel de gris.

La causa más probable es que `imgaussfilt` de MATLAB, con un núcleo de 161 puntos, filtre en
el dominio de la frecuencia (lo hace automáticamente con núcleos grandes), y su redondeo a
entero no coincida con el de la convolución directa en algunos píxeles.

Queda por decidir cómo tratarlo, y es una decisión de fondo: o se replica el camino exacto
de MATLAB, o se acepta la diferencia y se documenta, o se evita el problema trabajando en
coma flotante sin redondear a entero tras el filtro (lo que además elimina esta fragilidad
de raíz, pero deja de reproducir los valores antiguos).

## Estado

* Rama `humedad` con la cadena completa implementada y 251 pruebas en verde.
* Fase 1 (calibración) y fase 2 (imágenes, muestreo, modelo) terminadas y probadas.
* Falta: decidir lo del redondeo del filtro, integrar la fuente en el análisis y el comando
  propio, y aplicar los cambios acordados (umbral 0,95 y primer instante).
* Nada subido a GitHub; `main` tiene un commit local por delante del remoto.
