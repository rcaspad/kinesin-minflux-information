# B1v2 — Informe de ejecución sobre trayectorias reales (2026-08-13)

**Prereg v2 FIRMADO** (Raúl, 2026-08-13, tras veredicto A2′=PASA): SHA
post-firma `aa461a6e28b4aaa8198b3c44fd480ce0d7d14f6169e95a887460ef6836f22e6f`
(borrador congelado pre-A2′ `e4aebb2a…b756842`). Script:
`b1v2_analisis.py` (SHA `bc6d8bf5e97d95b560b2f961f7bb46fc6c969ba6a1b5a9a7f1dcaa99eea241e0`),
determinista (semillas tipo 7/8, raíz 20260812), estimadores de A2 sin
cambios. **Primera lectura de estructura temporal sobre los datos reales.**
Datos: Wolff et al. 2023, *Science* 379, 1004-1010 (Zenodo 7565676, CC-BY).

## 1. VEREDICTO FORMAL (criterio congelado): CONFIRMA en las 9 celdas

Canal evento, PT k=4, nulo de orden: **las 9 celdas (3 constructos × 3
[ATP]) rechazan** con p de permutación en el suelo (0,000999; p_Holm
0,002997). I_obs entre 0,061 y 0,513 bits, de 2,6 a 34 veces el umbral p95
del nulo. **Robusto a las dos sensibilidades obligatorias** (k=8 y clases de
paso por terciles: 9/9 en ambas). **El guard de limping se supera**: bajo el
nulo de paridad (que preserva la alternancia par/impar) las 9 celdas siguen
rechazando → la estructura NO es (solo) limping.

H1 del prereg v2 queda CONFIRMADA: la secuencia de eventos del caminar real
contiene información predictiva por encima del nulo de renovación.

## 2. Descomposición por canales — dónde vive la señal (registrada con el veredicto)

| canal | resultado | lectura |
|---|---|---|
| paso→paso | 9/9 celdas, fuerte | **persistencia de clase**: P(8→8)≈0,55-0,59, P(16→16)≈0,64-0,70 vs marginales 0,35-0,41/0,52-0,57 — consistente con la mecánica de subpasos CONOCIDA (un paso de 16 nm resuelto como dos subpasos de ~8 nm produce rachas de la clase 8) |
| dwell→paso | 9/9 celdas, fuerte | estructura cruzada del ciclo de paso (el dwell corto del subpaso intermedio precede al subpaso que completa) — también del ciclo conocido |
| dwell→dwell | **7/9 celdas** (falla E215C 1mM p=0,92; T324C 100µM p=0,26), modesto (0,006-0,044 bits) | candidato a memoria más allá del ciclo; sobrevive al nulo de paridad; correlación lag-1 mediana **negativa** (−0,07 a −0,19) |

## 3. Estatuto honesto del hallazgo (artefacto-vs-encuadre, alarma-8)

**Lo que este resultado ES:** la primera caracterización en términos de
información (estilo Still, lag 1) de trayectorias MINFLUX de un motor
molecular — el claim máximo declarado en el prereg v1 §2. La estructura
dominante que la información detecta es el **ciclo de paso/subpaso ya
descrito por los autores del dataset**: B1v2 lo REPLICA y lo CUANTIFICA
informacionalmente (0,06-0,51 bits/evento según celda y régimen), pero no lo
descubre.

**Lo que este resultado NO es (todavía):** evidencia de memoria mecanística
nueva. El componente candidato (dwell→dwell, 7/9 celdas) tiene una firma
—anticorrelación lag-1— que es **indistinguible, con este diseño, de un
artefacto de segmentación**: un error en la frontera de un paso alarga un
dwell y acorta el siguiente, produciendo exactamente anticorrelación
adyacente que sobrevive al nulo de paridad. Confounders no excluibles
declarados: (a) error de fronteras correlacionado del pipeline de pasos de
los autores; (b) no-estacionariedad dentro de traza (tendencias lentas
cuentan como "memoria" para el nulo de permutación); (c) sesgos de selección
de las tablas de pasos (heredados del análisis primario, prereg v1 §6).

**Diagnósticos post-hoc, no confirmatorios:** matrices de transición y
correlaciones lag-1 de este informe se computaron DESPUÉS del veredicto,
como caza de artefactos mandada por la regla de integridad; no forman parte
del contraste preregistrado.

## 4. Qué haría falta para ascender el componente dwell→dwell

1. Sensibilidad sobre trayectorias crudas con detección de pasos
   INDEPENDIENTE (prereg v1 §6 la contempla) — si la anticorrelación es del
   segmentador, cambiarlo la mueve; si es del motor, persiste.
2. El **verificador externo humano (gate 4)** — obligatorio antes de
   cualquier canonización; la carta a Ariga/escuela Sivak sigue pendiente.
3. B2 (nostalgia_op vs [ATP]) exigiría su adenda pública de validación de
   I_mem (prereg v2 §4); no se ejecuta sin ella ni sin decisión de Raúl.

## 5. Archivos

- `b1v2_analisis.py` — SHA `bc6d8bf5…ea241e0`
- `b1v2_resultados.csv` — todas las combinaciones (2 variantes × 2 nulos ×
  4 canales × 2 k × 4 estimadores × 9 celdas), SHA en commit.
- Tablas completas del primario y las sensibilidades: en el propio stdout
  del script (determinista, re-ejecutable).

Ejecutado y firmado: Agente Khora, 2026-08-13. Veredicto formal: CONFIRMA.
Encuadre: replicación informacional del ciclo conocido + candidato
dwell→dwell pendiente de discriminación artefacto/motor.
