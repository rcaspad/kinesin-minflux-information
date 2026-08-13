# B1 — Congelación previa a la ejecución (2026-08-13)

**Prereg:** v1.0 FIRMADO, SHA `c73041ab…cc6c71b2` (+ errata v1.1).
**Decisión del gobernador (Raúl, 2026-08-13, en sesión):** "adopto la lectura
binomial, elige tú el estimador y ejecuta B1". Con ello **A2 PASA** — la
lectura binomial se adopta como interpretación del criterio congelado (i)
(documentada como interpretación en `../a2/RECEPCION_A2.md` §4-5, no como
cambio del prereg).

## 1. Estimador elegido (delegado en el agente, congelado aquí)

**Panzeri-Treves con k = 4** para la decisión primaria. Razón: empate exacto
NSB(k=4) ↔ PT(k=4) en los tres criterios de desempate de A2; a igualdad
empírica se aplica el criterio declarado por el propio informe A2 §8 — PT es
computacionalmente más barato, y la ventaja teórica de NSB (menos sesgo
residual con alfabetos grandes) es invisible en la rejilla validada. k = 8 se
reporta SIEMPRE junto a k = 4 (prereg §4); plug-in, Miller-Madow y NSB se
registran como diagnóstico en el por-instancia, sin papel decisorio.

## 2. Operacionalización sobre datos reales (declarada ANTES de computar)

Todo lo siguiente queda fijado antes de calcular ningún estadístico de B1.

- **Fuente primaria:** tablas de pasos de los autores (`allsteps_reeval.xls`
  por constructo × [ATP]), semántica de columnas decodificada de su
  `process_MF_data.m` L71 (ver A1). Análisis por constructo SEPARADO
  (E215C, K28C, T324C; N356C es control estacionario, sin caminar — errata
  v1.1), nunca mezclados (prereg §6).
- **Evento:** par (dwell `tau`, tamaño de paso `stepx`) de cada fila con
  `end_flag == 0` y `tau > 0`, en el orden del fichero dentro de cada traza
  (las filas `end_flag == 1` son terminadores de traza).
- **Clases de paso (3, nativas del pipeline validado):** |stepx| se asigna a
  {4, 8, 16} nm por centro más cercano → bordes en 6 y 12 nm. Elección
  física (picos 16/8/4 nm del paper), determinista y sin ajuste; declarada
  aquí antes de mirar histogramas por traza. Sensibilidad: terciles
  empíricos por celda (secundaria, se reporta si cambia el veredicto).
- **Bins de dwell:** cuantiles agrupados sobre TODAS las trayectorias
  elegibles de la celda (constructo × condición), k = 8 nativo y k = 4 por
  engrosamiento exacto — espejo de `a2_estimadores.simbolizar`.
- **Estadístico por trayectoria (prereg §3-B1):** I[evento_k ; evento_{k+1}]
  sobre los pares consecutivos de la traza; nulo por permutación del orden
  de eventos dentro de la traza, 1.000 permutaciones; rechaza si
  I_obs > percentil 95 del nulo propio.
- **Elegibilidad:** trazas con ≥ 50 eventos (prereg §4, congelado).
  INCONCLUSO si N elegible < 10 en todas las condiciones.
- **Semillas:** `semilla_de(3, constructo_id, cond_id, tray_id)` para el
  nulo primario; `semilla_de(4, constructo_id, cond_id, tray_id)` para el
  nulo de paridad. constructo_id: E215C=0, K28C=1, T324C=2; cond_id:
  10µM=0, 100µM=1, 1mM=2; tray_id = índice de traza en el orden del fichero.
- **Criterio de confirmación (operacionalización del prereg §3-B1):** por
  condición, fracción de trazas elegibles que rechazan. p-valor de condición
  = P[Binomial(n_eleg, 0.05) ≥ observados] (cola superior, tasa nominal del
  test por traza), corregido por **Holm entre las 3 condiciones dentro de
  cada constructo**. CONFIRMA si alguna condición tiene fracción ≥ 60 % Y
  p_Holm < 0.05. REFUTA si ninguna condición de ningún constructo llega.
  INCONCLUSO según elegibilidad (arriba).
- **Guard de limping (innegociable, prereg §3-B1 y §5):** recomputación
  íntegra con el nulo de paridad — permutación del orden SOLO dentro de las
  posiciones pares e impares de cada traza (preserva la estructura de
  alternancia, destruye el resto del orden). Si el rechazo desaparece bajo
  este nulo (la condición ya no cumple el criterio), el resultado se reporta
  como **replicación de limping con datos MINFLUX**, no como hallazgo.
  Cita canónica verificada: Asbury, Fehr & Block 2003 (corpus #23).
- **Advertencia de alcance (registrada antes de ejecutar):** A2 validó la
  calibración con instancias agrupadas (n_tray ∈ {50,100,300}); el test por
  traza individual hereda su exactitud de la construcción por permutación
  (válida para cualquier n), pero las potencias de A2 NO se transfieren al
  caso por-traza — a igual magnitud, la potencia por traza es menor. El
  umbral del 60 % de trazas es el listón del prereg, no una potencia
  estimada.

## 3. Qué NO se hace

Sin decisión posterior de Raúl: nada de B2, nada de sensibilidad sobre crudo
(detección de pasos propia), nada de mezclar constructos, ningún ascenso de
resultado a canon (gate 4 — verificador externo — sigue pendiente y es
previo a B2).

Congelado por el agente Khora, 2026-08-13, antes de abrir las tablas para B1.
