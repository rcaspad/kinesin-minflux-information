# Etapa E — Desempate de la identidad de sitio: congelación previa (2026-08-13)

**Mandato:** Raúl en sesión ("adelante con la etapa E"). **Estatuto:
ATRIBUTIVA**, continuación de D (veredicto NO RESUELTA). Nada asciende.
**Pregunta:** la identidad de sitio de D-1, ¿es del raíl (sitios con
carácter), de la dinámica/fragmentación de la cabeza marcada, o del ciclo
químico del motor?

**Verificación previa registrada:** el paper define `stabilized_microtubules/`
como *conjuntos de trazas de K28C sobre 19 microtúbulos individuales en
muestras ACTIVAMENTE estabilizadas* (deriva corregida → marco de
coordenadas compartido por microtúbulo; clases de protofilamento
'sides'/'center'/'between'; Fig. 5C del paper). La premisa del Nodo
(química de estabilización distinta) era FALSA — cuarto error registrado —
pero el recurso real es superior: habilita la PRUEBA REINA (varios motores
sobre el mismo raíl). Solo se han inspeccionado cabeceras y conteos.

**Base:** cadena KV completa de la etapa C (SHA `be3756b8…8c4d84d6`), sin
cambios. Semillas raíz 20260812: tipo 13 (E-1 orden), 14 (E-1 paridad),
15 (E-1 revisitas), 16 (E-2 permutaciones). 1.000 permutaciones; α=0,05
unilateral en la dirección predicha.

## E-1 — Desempate cabeza-vs-COM (N356C caminando; datos nunca abiertos)

Celdas: N356C DOL1/10µM (id 9), DOL1/1mM (id 10), DOL2/1mM (id 11). El
marcaje es de tallo (centro de masas): SIN vaivén de cabeza ni subpasos de
cabeza; pasos de COM ~8 nm.

- **Gate de comparabilidad** (como C): eventos KV vs tablas
  `allsteps_reeval.xls` de cada celda, 50-200 %. Celda fuera → no
  participa. Las 3 fuera → E-1 INCONCLUSO.
- **E-1a:** IM dwell→dwell (PT k=4) agrupada por celda, nulos de orden y
  paridad (como C-1).
- **E-1b:** test de revisitas D-1 idéntico (mismos umbrales 4/8 nm,
  j−i≥3, ≥30 pares).
- **Predicciones congeladas:**
  - *Fragmentación/vaivén de cabeza* → E-1a débil o nula Y E-1b nulo
    (el COM no parpadea).
  - *Raíl (sitios con carácter)* → E-1b PERSISTE en COM (T<0, p<0,05)
    (el sitio es el mismo lo mires desde donde lo mires).
  - *Motor (ciclo químico)* → E-1a persiste (los dwells son el mismo
    ciclo) con E-1b nulo.

## E-1c — Periodicidad de dos pasos en COM (añadido pre-ejecución)

Propuesto por el Nodo (su "E-3"), incorporado ANTES de ejecutar nada.
Estadístico D-2 (D = ρ_cerca − ρ_lejos con corte 12 nm, log-dwell centrado
por traza, permutación de etiquetas tipo 12→reusa semilla tipo 14 con
etiqueta extra 1) aplicado a las celdas COM de E-1.
- *H_ciclo* (la organización por tipo de paso es del ciclo de dos
  cabezas): D < 0 también en COM.
- *H_cabeza* (es del vaivén de la cabeza marcada): D ≈ 0 en COM.
Nota de auditoría: el "E-2 por química de estabilización" del Nodo se
RECHAZA — premisa verificada como falsa (los ROIs son estabilización
ACTIVA de muestra, no química; quinto pinneo) — y su "E-4, quizá no
ejecutable" es exactamente el E-2 de esta congelación, que SÍ lo es.

## E-2 — PRUEBA REINA (trazas K28C sobre el mismo microtúbulo)

Por cada ROI (microtúbulo): cargar todas sus trazas (marco compartido),
ajustar el EJE COMÚN por regresión sobre las localizaciones agrupadas de
todas las trazas, segmentar cada traza con la cadena KV, y asignar a cada
evento su posición AXIAL en el marco común (proyección de la posición
media del tramo).

- **Estadístico (cruzado entre trazas):** para cada evento con ≥1 evento
  de OTRA traza del mismo ROI a < 4 nm axiales, par (log-dwell centrado
  por traza, media de log-dwells centrados de los otros en ese sitio).
  ρ = Spearman agrupado sobre todos los ROIs.
  **H_raíl predice ρ > 0** (motores distintos se atascan en los mismos
  sitios). H_motor y H_cabeza predicen ρ ≈ 0.
- **Nulo:** permutar los dwells dentro de cada traza (posiciones fijas),
  1.000 veces, recomputando ρ. p unilateral.
- **Estrato secundario declarado:** pares con separación LATERAL < 6 nm
  (probable mismo lado/protofilamento) vs todos — un defecto local de
  protofilamento solo aparece en el estrato estrecho.
- **INCONCLUSO** si < 50 pares cruzados en total.

## Reglas de lectura

Tabla de atribución final combinando D y E; contradicciones → NO RESUELTA
→ gate 4 (sin etapa F sobre este dataset: se acabaron las rendijas).
Los resultados se reportan íntegros pase lo que pase.

Congelado por el agente Khora, 2026-08-13, antes de computar dwells de
N356C o de los ROIs.
