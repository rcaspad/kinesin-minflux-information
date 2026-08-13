# Etapa D — Atribución del residuo: congelación previa (2026-08-13)

**Mandato:** Raúl en sesión ("congela la etapa D y ejecutala").
**Estatuto de la etapa: ATRIBUTIVA.** No confirma ni refuta H1 (eso quedó
en B1v2/C); asigna el residuo dwell→dwell de las 3 celdas supervivientes
(E215C 10µM, K28C 10µM, T324C 1mM) entre las cuatro fuentes de la
taxonomía (motor / raíl / instrumento / ambiente). Nada asciende a canon.
**Predicciones pre-dato:** registradas en `../../DEBATE_RAIL_20260813.md`
§4 (commit dcc1ad78, ANTERIOR a esta congelación).

**Base de datos:** secuencias KV de la etapa C (segmentador independiente,
fusión <5nm, pasos CON SIGNO), 8 celdas comparables; foco atributivo en
las 3 supervivientes, resto como contexto. Posición relativa del evento m:
pos_m = suma de pasos firmados anteriores (marco local de traza).
**Semillas:** `semilla_de(12, test_id, id_celda_real)`, raíz 20260812;
1.000 permutaciones; α = 0,05 unilateral en la dirección predicha.

## D-1 — Revisitas (test_id=1): ¿el sitio tiene identidad?

- Par revisita: eventos i<j de la misma traza con j−i ≥ 3 y
  |pos_i − pos_j| < 4 nm (mismo sitio de red, con ≥2 eventos por medio —
  guard anti-parpadeo).
- Δ_same = |log(dwell_j/dwell_i)|. Controles emparejados: pares (i′,j′) de
  la MISMA traza con j′−i′ = j−i y |pos_i′ − pos_j′| ≥ 8 nm (sitios
  distintos, misma separación temporal-de-eventos).
- Estadístico de celda: T = mediana(Δ_same) − mediana(Δ_diff).
  **H_raíl predice T < 0** (pausas más parecidas en el mismo sitio).
- Nulo: 1.000 permutaciones reasignando la etiqueta same/diff dentro de
  cada (traza, separación). p unilateral. INCONCLUSO si < 30 pares
  revisita en la celda.

## D-2 — Decaimiento espacial vs temporal (test_id=2)

- Pares lag-1 (log dwell_k, log dwell_{k+1}) con log-dwell centrado por
  traza (resta de la media de la traza — mata heterogeneidad entre trazas).
- Estratos por salto espacial intermedio: |paso_k| ≤ 12 nm («cerca», clase
  4/8) vs > 12 nm («lejos», clase 16).
- Estadístico: D = ρ_cerca − ρ_lejos (Spearman agrupado por celda).
  **H_raíl predice D > 0** (la similitud muere con la distancia).
- Nulo: 1.000 permutaciones de las etiquetas cerca/lejos dentro de traza.
  p unilateral.

## D-3 — Rachas: árbol, barro o caminante (test_id=3)

- Por traza: lento = log dwell > mediana de la traza. Rachas de «lento»:
  longitud en eventos (L_ev) y extensión en nm (L_nm = Σ|paso| interno).
- Estadístico: media de L_ev por celda vs nulo de orden barajado dentro de
  traza (1.000 perms). p unilateral (exceso).
- Guía interpretativa congelada (no veredicto automático): exceso con
  L_nm concentrado en una escala → parche (barro); sin exceso pero con
  outliers aislados de dwell → defecto puntual (árbol); exceso cuya L_ev
  varía con [ATP] dentro de constructo → motor. Se reportan las tres
  lecturas con números.

## D-4 — Control estacionario N356C (test_id=4)

- La MISMA cadena (carga cruda → filtro → proyección → KV → fusión <5nm →
  eventos) aplicada a `N356C/DOL1` y `N356C/DOL2` (motor clavado: todo
  «paso» detectado es instrumento puro).
- Se reporta: nº de trazas, eventos por traza (esperado ≈ 0 si la cadena
  no fabrica pasos), y si ≥ 30 trazas con ≥ 2 eventos: IM dwell→dwell con
  permutación (como C-1). Cota de instrumento sobre ruido REAL,
  complementaria a las sustitutas sintéticas de C-2.

## Reglas de lectura de etapa

Los cuatro tests se reportan íntegros, celda a celda, pase lo que pase.
La atribución final es la combinación coherente de los cuatro (tabla en el
informe); si los tests se contradicen entre sí, la etapa declara
ATRIBUCIÓN NO RESUELTA y el paquete pasa tal cual al gate 4. Sin D-v2:
lo que no salga de aquí, lo decide el verificador externo o Wirth.

Congelado por el agente Khora, 2026-08-13, antes de computar dwells para
estos tests (solo conteos de viabilidad previos: nº de retrocesos y
revisitas brutas, registrados en conversación).
