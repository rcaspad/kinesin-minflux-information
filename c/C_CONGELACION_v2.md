# Etapa C — Iteración v2 de la congelación (2026-08-13)

## Acta de C-v1: INCONCLUSA por su propio gate

El gate de comparabilidad de `C_CONGELACION.md` §2 (50-200 % de los pasos
de los autores) dejó **3 celdas comparables de 9** (K28C ×3; ratios del
resto 2,24-4,35) — menos de las 5 exigidas → C-v1 INCONCLUSA. Se registra
que las 3 comparables coinciden con el constructo que confirmaba
dwell→dwell en las 3 condiciones de B1v2.

**Diagnóstico (con solo conteos y marginales vistos; NINGÚN estadístico
temporal computado):** KV sin umbral de tamaño resuelve subestructura
< 5 nm que el pipeline de los autores fusiona (su `min_step_h = 5`). Las
tablas de B1v2 definen el evento a granularidad ≥ 5 nm; comparar contra
secuencias KV de granularidad más fina es comparar objetos distintos. El
fallo del gate es de definición de evento, no del motor de segmentación.

## Cambio único de C-v2 (declarado y congelado antes de computar nada temporal)

Tras la segmentación KV (idéntica a C-v1), **fusión iterativa de tramos
adyacentes cuyo salto de nivel sea < 5 nm** (se funde el par con menor
salto, nivel = media ponderada por longitud, y se repite hasta que no
queden saltos < 5 nm). El umbral 5 nm se toma de los autores como
DEFINICIÓN COMÚN DE EVENTO; el motor de detección (BIC voraz sin suavizado)
sigue siendo independiente del suyo (`ischange` 12σ² + mediana móvil 9).

Todo lo demás de `C_CONGELACION.md` queda intacto: gate §2 re-aplicado tras
la fusión, contrastes C-1 y C-2 (con la fusión aplicada TAMBIÉN a las
sustitutas del nulo de artefacto — misma cadena completa), semillas 9/10/11,
regla de decisión §4, y limitaciones declaradas.

No habrá C-v3: si el gate vuelve a fallar, la etapa se cierra INCONCLUSA
definitivamente y la discriminación pasa al verificador externo (gate 4).

Congelado por el agente Khora, 2026-08-13.
