# Casos sintéticos de validación

Campos de velocidad fabricados a mano sobre mallas pequeñas, con solución conocida. Sirven
para comprobar el cálculo en sí, no solo que no cambie respecto a versiones anteriores.
`casos.json` guarda el nombre original de cada caso en la colección del grupo.

Cada carpeta contiene los archivos PIVlab (`datos (n).txt`), el `.PAR` **en el formato
antiguo** (bloque 3 con 3 valores y sin bloque 4, con lo que también se comprueba que se
sigue leyendo) y, en `expected/`, el `.POST.RES` del análisis que el grupo hizo en su día.

| Caso | Malla | Imposición | Resultado esperado |
|---|---|---|---|
| `desplazamiento_1P`, `_4P` | 4×4 y 5×5 celdas de 0.5 m | 10 s a 2 cm/s en +x y 2 cm/s en −y en todos los nodos | u = (+0.20, −0.20) m en todas las partículas, deformación nula |
| `corte_1P`, `_4P` | 5×4 celdas | 10 s con velocidad horizontal de 2, 1.5, 1, 0.5 y 0 cm/s de arriba a abajo | γxy = 0.10 uniforme, εxx = εyy = 0, εq = 0.1/√3 |
| `horizontal_uniforme_1P`, `_4P` | 5×4 celdas | 24 s de estiramiento uniforme | εxx = 0.24 uniforme, εyy = γxy = 0 |
| `horizontal_1P`, `_4P` | 5×4 celdas | 15 s con velocidad de 0, 0.1, 0.3, 0.9 y 2.4 cm/s por columna | εxx distinta en cada columna de celdas: 0.03, 0.06, 0.18… |
| `rotacion_1P`, `_4P` | 6×6 celdas | 50 s girando 1°/s alrededor del nodo central | el sólido no se deforma, así que la deformación debería ser 0 |

En varios casos algunas partículas salen de la malla y dejan de calcularse; por eso el
número de partículas activas al final es menor que el inicial.

## Sobre la rotación

Una rotación de sólido rígido no deforma el material, así que toda la deformación que sale
es error. Midiendo el caso se ve que tiene **dos orígenes distintos**:

| | Deformación de corte equivalente (media) |
|---|---|
| Sin corrección de contorno | 0.0169 |
| Media de los nodos vecinos | 0.0111 |
| Media de las partículas de alrededor | 0.0125 |
| Extrapolación desde el interior | **0.0051, igual en todas las partículas** |

El caso tiene instantes con 24 de sus 49 nodos sin dato, así que **dos tercios del error
venían del contorno**. Al reconstruir esos nodos por extrapolación, el error restante es
uniforme en todo el sólido y coincide con el que predice la teoría: acumular incrementos de
deformación lineales durante un giro finito deja εxx = εyy = n·(cos Δθ − 1), que para 50
pasos de 1° da una deformación de corte equivalente de 0.00508, frente a los 0.0051 medidos.

Ese residuo sí es una limitación de la formulación, no de esta implementación: eliminarlo
requeriría una medida de deformación finita (a partir del gradiente de deformación). Las
pruebas fijan los dos valores para detectar si algún cambio los empeora.
