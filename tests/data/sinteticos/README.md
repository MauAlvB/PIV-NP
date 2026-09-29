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

Una rotación de sólido rígido no deforma el material, pero acumular incrementos de
deformación lineales sí produce una deformación aparente, que crece con el ángulo girado
(unos 0.044 tras 50°) y no depende del número de partículas por celda. Es una limitación de
la formulación, no un problema de esta implementación: eliminarla requeriría una medida de
deformación finita (a partir del gradiente de deformación). La prueba correspondiente fija
el valor actual para detectar si algún cambio lo empeora.
