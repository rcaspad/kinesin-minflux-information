# Preregistro — Estructura informacional del caminar de la kinesina-1 (F1-B sobre datos MINFLUX)

**Estado:** `FIRMADO Y CONGELADO — v1.0 — 2026-08-12`
**Firma:** Raúl, autorización expresa en sesión ("firmo, ejecuta los gates y
arranca A1"). Este documento queda inmutable; cambios exigen v2 pública.
**SHA-256 del borrador pre-firma:** `e1b216f58a85ebacff1d6fe3781dcfbf0a3dbab5e613d052721c77397ba2ef9c`
**Fecha de redacción:** 2026-08-12
**Redactor:** Agente Khora. **Gobernador:** Raúl (firmado 2026-08-12).
**Estudio:** 14 — Termodinámica de la Información, FASE 1, gate 3.

**Declaración de integridad de partida:** en el momento de redactar este
prereg, NINGÚN archivo de trayectorias del dataset ha sido abierto. La
información usada se limita a: el inventario `unzip -l` (nombres y tamaños,
en `../raw/inbox/kinesin_minflux_wolff2023_zenodo7565676/MANIFEST.md`), el
paper compañero (Wolff et al., *Science* 379, 1004-1010, 2023) y el corpus
verificado del estudio. Tras la firma de Raúl este documento se congela
(SHA-256 registrado en `memoria/decisiones.md`) y cualquier cambio exige v2
explícita y pública.

---

## 1. Datos

- **Dataset:** Zenodo `10.5281/zenodo.7565676` (Wolff & Scheiderer, MPI;
  CC-BY-4.0), depositado con MD5 verificado y SHA-256
  `69a335b5…dfb9a1930`. Trayectorias MINFLUX de kinesina-1 caminando sin
  carga sobre microtúbulos, ~1,7 nm / 1 ms.
- **Estructura:** constructos con marcado en tallo (E215C) y en cabeza
  (K28C), tres concentraciones de ATP (10 µM, 100 µM, 1 mM), muestras
  fechadas; tablas de pasos de los autores (`allsteps_reeval.xls`) y scripts
  Matlab propios.
- **Regla raw:** el ZIP es inmutable; la extracción de trabajo se hace en
  `analisis/` de esta fase, nunca dentro de `raw/`.

## 2. Pregunta y estatuto del claim

**Pregunta:** ¿contiene el caminar real de la kinesina estructura
informacional temporal — memoria que predice y memoria que no ("nostalgia"
operacional) — y cómo cambia entre regímenes de alimentación de ATP?

**Estatuto máximo alcanzable declarado de antemano:** primera
caracterización empírica de la estructura informacional (estilo Still) de
trayectorias MINFLUX de un motor molecular. **NO es un test del teorema de
Still** — sus supuestos (señal exógena registrada, sin feedback) no se
cumplen aquí, y así se declarará en cualquier salida. NO valida ninguna
tesis de Khora; es empiria acotada del Estudio 14.

## 3. Diseño por etapas (cada una con criterio de muerte)

### Etapa A1 — Gate de reproducibilidad (obligatoria, sin claim)

Reproducir con los scripts y tablas de los autores los estadísticos básicos
del paper (tamaños de paso, dwell times por condición).
- **PASA** si nuestras lecturas de las tablas/trayectorias reproducen los
  histogramas de los autores (concordancia visual + medianas dentro de ±10 %).
- **MUERE el pipeline** (no el frente) si no: se documenta la discrepancia y
  no se computa nada más hasta resolverla.

### Etapa A2 — Validación del estimador en sintético (antes de tocar dato real)

Los estimadores de información (sec. 4) se validan sobre datos sintéticos de
un proceso de renovación semi-Markov calibrado SOLO con los parámetros
publicados en el paper (velocidades, tamaños de paso, dwells medios por
condición — nunca con los archivos).
- **PASA** si el estimador recupera I_pred = 0 en el nulo de renovación
  (falsos positivos ≤ 5 %) y detecta memoria inyectada conocida con potencia
  ≥ 80 % a los N disponibles.
- **MUERE la formulación** si el estimador no pasa: se reporta "no medible
  con este dataset/estimador" y el frente se replantea sin tocar los datos.

### Etapa B1 — CONFIRMATORIA: memoria más allá de la renovación

**H1:** la secuencia de eventos del caminar — pares (dwell time, tamaño de
paso) — contiene información predictiva por encima del nulo de renovación
(eventos i.i.d. con las mismas marginales), en al menos una condición de ATP.

- Estadístico por trayectoria: I_pred(1 evento) = I[evento_k ; evento_{k+1}]
  estimada según sec. 4; significación contra el nulo por permutación
  (sec. 5).
- **Guard de prior-art (anti-alarma-8), innegociable:** el "limping"
  (alternancia par/impar de dwells en kinesina) es fenómeno CONOCIDO de la
  literatura. Ablación obligatoria: recomputar H1 condicionando a la paridad
  del paso. Si la memoria desaparece al condicionar → el resultado se
  reporta como **replicación de limping con datos MINFLUX**, no como
  hallazgo. Las citas canónicas de limping se verifican con gate ANTES de
  ejecutar B1 [PENDIENTE-GATE en sec. 8].
- **CONFIRMA** si ≥ 60 % de las trayectorias elegibles de una condición
  superan el percentil 95 del nulo (con Holm entre las 3 condiciones).
- **REFUTA** si ninguna condición lo logra. **INCONCLUSO** si N elegible
  < 10 trayectorias en todas las condiciones.

### Etapa B2 — EXPLORATORIA→condicionada: nostalgia operacional vs. régimen de ATP

**H2 (direccional, registrada ahora):** definiendo el estado coarse-grained
s_k (fase del ciclo de paso según el esquema del paper) y
nostalgia_op = I_mem − I_pred respecto al siguiente evento, la nostalgia_op
por evento **aumenta con la concentración de ATP** (10 µM < 100 µM < 1 mM).
*Racional registrado:* a ATP saturante el motor está más lejos del
equilibrio (mayor fuerza motriz) y el paso mecánico domina sobre la espera
markoviana de ATP; si el marco de Still tiene tracción aquí, más driving
debería dejar más información retenida no predictiva. Es una apuesta
falsable, no un teorema.
- Solo asciende a confirmatoria si A2 valida el estimador de I_mem con los
  mismos criterios. Si no, queda descriptiva.
- **CONFIRMA** la ordenación si es monótona en las medianas por condición
  con IC95 bootstrap por trayectoria que separe los extremos.
  **REFUTA** si la ordenación es la inversa. **INCONCLUSO** en otro caso.

## 4. Estimadores (declarados; elección congelada tras A2)

- Eventos discretos: I por tablas de contingencia con corrección de sesgo de
  Panzeri-Treves o, alternativamente, estimador NSB — la elección entre
  ambos se hace EN A2 (sintético) y queda congelada antes de B1.
- Variables continuas (dwells): discretización por cuantiles (k=4 y k=8,
  ambas reportadas) declarada aquí; sensibilidad obligatoria a k.
- Longitud mínima de trayectoria elegible: ≥ 50 eventos (registrado ahora).

## 5. Nulos y aleatorización

- Nulo primario: permutación de eventos dentro de trayectoria (rompe orden,
  preserva marginales) — 1.000 permutaciones por trayectoria.
- Nulo de paridad (guard de limping): permutación dentro de paridades.
- Semilla raíz: `20260812`; derivación determinista por (condición,
  trayectoria).
- **Por-instancia:** todo se reporta por trayectoria y por condición; el
  agregado nunca sustituye al desglose.

## 6. Confounders declarados

Ruido de localización y photobleaching (filtros de los autores heredados en
análisis primario; sensibilidad con detección de pasos propia), deriva
instrumental, heterogeneidad entre muestras/fechas (efecto fecha como
covariable de sensibilidad), sesgo de selección de las tablas de pasos de
los autores (análisis primario sobre SUS tablas para heredar su validación;
sensibilidad sobre crudo), diferencias entre constructos (E215C y K28C se
analizan por separado; nunca se mezclan).

## 7. Qué NO puede concluir este análisis (registrado de antemano)

- Nada sobre el teorema de Still como tal (supuestos no satisfechos).
- Nada sobre eficiencias termodinámicas (sin medición de fuerzas ni
  contabilidad de ATP por paso en estos datos).
- Nada universal sobre "los motores" — una proteína, un ensayo, sin carga.
- Ningún ascenso a canon sin el verificador externo del gate 4.

## 8. Gates previos a la ejecución (además de la firma)

1. [PENDIENTE-GATE] Verificar citas canónicas del limping y de correlaciones
   de dwell en kinesina (búsqueda + CrossRef) e incorporarlas al corpus.
2. [PENDIENTE-GATE] Búsqueda dirigida: ¿alguien ha computado ya información
   predictiva sobre trayectorias MINFLUX de motores (2023-2026)? Si existe →
   regla de cancelación.
3. Gate 4 en marcha: el prereg congelado se enviará al verificador externo
   identificado (candidatos: Ariga/Osaka; escuela Sivak) ANTES de ejecutar
   B2, o como muy tarde antes de canonizar resultado alguno.

## 9. Atribución y salida

Todo resultado cita: Wolff, Scheiderer, Engelhardt, Engelhardt, Matthias &
Hell (2023), *Science* 379, 1004-1010 (datos CC-BY-4.0, Zenodo 7565676).
Resultados —confirmen o refuten— se depositan en
`resultados/` de esta fase con código determinista, semillas y SHA; el
informe final declara este prereg y su SHA post-firma.

## 10. Firma del gobernador

- [x] **FIRMADO** — Raúl, 2026-08-12, en sesión. Se ejecutan A1→A2 y, si
  pasan, B1→B2.

*Tras la firma: recalcular SHA-256, registrar en `memoria/decisiones.md`,
y solo entonces abrir el primer archivo de datos.*
