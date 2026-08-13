# Acta de recepción y verificación — Etapa A2 (2026-08-13)

**Prereg:** v1.0 FIRMADO Y CONGELADO, SHA `c73041ab…cc6c71b2` (+ errata v1.1).
**Origen:** delegación a Claude Science (proyecto `proj_c9c0e12c998d`, sesión
"Synthetic Data Validation Pipeline (A2)", Opus 5) — ver
`../DELEGACION_A2_CLAUDE_SCIENCE.md`. Estatuto: los resultados del delegado
son disparadores; el veredicto lo fija la verificación en este repo.

## 1. Qué se recibió

Barrido completo: 594 tareas (54 celdas distintas × 11 condiciones; las 216
celdas del prereg se recuperan por la degeneración exacta del eje
`dwell_medio_ms`, ver `A2_informe.md` §2), 500 réplicas por celda×condición,
1.000 permutaciones por instancia, semilla raíz 20260812. Ejecutado en la
máquina local (9 procesos, ~6,5 h, sin coste de tokens en el cómputo).
16 archivos depositados en este directorio; inventario y SHA-256 en
`A2_informe.md` §11.

## 2. Verificación de integridad (bit a bit)

- Los 16 archivos del workspace del delegado coinciden con la tabla SHA-256
  del propio `A2_informe.md` (§11), incluidos los dos grandes:
  `a2_por_instancia.csv.gz` (548.744.455 bytes, `38d36f49…13f0839`) y
  `a2_por_instancia_54celdas.parquet` (114.491.231 bytes, `0012b4da…12f37c0`).
- Tras la copia a este directorio se recalcularon TODOS los SHA: idénticos.
- Discrepancia detectada y resuelta durante la recepción: `a2_barrido.py`
  fue editado por el delegado (corrección de un comentario, hallazgo de su
  Reviewer) DESPUÉS de generar el informe; el informe se regeneró y la tabla
  §11 quedó consistente antes de la recepción. Las salidas regeneradas
  (CSV, veredicto, tablas) conservaron SHA idéntico — regeneración
  determinista confirmada de paso.

## 3. Re-ejecución independiente ("manos, no juez")

Entornos deliberadamente distintos:

| | delegado (Claude Science) | este repo (.venv) |
|---|---|---|
| Python | 3.13.15 | 3.14.6 |
| numpy | 2.4.6 | 2.4.4 |
| scipy | 1.18.0 | 1.17.1 |

- **Suite `a2_pruebas.py` (13 comprobaciones): TODAS PASAN** en el entorno
  del repo (determinismo, reanudabilidad con hash, marginales del generador,
  pasos i.i.d. en toda condición, degeneración del eje, engrosamiento
  k=4/k=8, estimadores contra la uniforme, convergencia NSB).
- **Re-ejecución de bloques desde la semilla: 4/4 IDÉNTICOS BIT A BIT**
  (2.000 réplicas recomputadas: celdas 0054 y 0063 × {H0, limping δ=0,5,
  AR(1) ρ=0,5}) contra los bloques `.npy` del barrido original.
  Script y log: `verificacion/`. La igualdad bit a bit entre versiones
  distintas de Python/numpy/scipy es el resultado fuerte de esta acta.

**VEREDICTO DE RECEPCIÓN: ÍNTEGRO Y REPRODUCIDO.**

## 4. Dictamen de A2 (del informe, sin editar)

> NO MEDIBLE bajo el criterio congelado en su lectura literal (alguna celda
> supera FPR 0.05); MEDIBLE con reservas bajo la lectura binomial. La
> discrepancia se documenta en §3.1 y la decide el responsable del prereg.

Los tres hechos que el responsable necesita para decidir:

1. **Los estimadores no fabrican señal.** Control interno decisivo: el canal
   «paso» es i.i.d. por construcción en TODA condición, y aun con memoria
   inyectada en los dwells su tasa de rechazo es 0,0501 (2.160
   celdas×condición, los 4 estimadores). FPR bajo H0: media 0,048-0,050.
2. **El criterio literal es insatisfacible por construcción.** Con 500
   réplicas (EE=0,010), un estimador PERFECTAMENTE calibrado supera 0,05 en
   ~la mitad de las celdas por azar. Test binomial exacto + BH sobre las 54
   celdas distintas: **0 celdas** significativamente sobre 0,05 en los 8
   pares (estimador, k). El delegado no relajó el umbral: reportó ambas
   lecturas.
3. **La potencia depende de N, no del estimador.** En régimen realista
   (magnitud ≤ 0,20): 9 % de celdas alcanzan potencia 0,80 con limping,
   63 % con AR(1). Mínima detectable mediana: ρ ≈ 0,05-0,10; δ ≈ 0,20-0,35.
   Empate exacto NSB(k=4) ↔ Panzeri-Treves(k=4) en los tres criterios de
   desempate; la elección entre ellos queda abierta y es del responsable.

## 5. Decisión pendiente (Raúl, gobernador del prereg)

- **(a)** Adoptar la lectura binomial como interpretación correcta del
  criterio congelado (documentándolo como interpretación, no como cambio):
  A2 PASA → congelar estimador+k y ejecutar B1.
  Si se adopta: falta elegir entre NSB(k=4) y Panzeri-Treves(k=4).
- **(b)** Sostener la lectura literal: A2 = NO MEDIBLE → muere la
  formulación, el frente se replantea sin tocar los datos (prereg §3-A2).

Nada de B1 se ejecuta hasta esa decisión. Gate 4 (verificador externo
humano) sigue pendiente y es previo a B2 en todo caso.

## 6. Notas de coste y limpieza

- El cómputo del barrido fue CPU local; el gasto de tokens de la sesión se
  limitó a los chequeos periódicos de progreso nocturnos.
- Los bloques crudos (`bloques_a2/`, 798 MB) permanecen en el workspace del
  delegado; no se depositan porque `a2_por_instancia.csv.gz` es su
  concatenación exacta y cualquier bloque es regenerable desde la semilla
  (demostrado en §3). Los dos ficheros por-instancia grandes están en disco
  en este directorio pero EXCLUIDOS de git (SHA fijados aquí y en §11 del
  informe).
- El delegado ofreció empaquetar el flujo como skill reutilizable; queda
  sin responder a la espera de Raúl.

Recibido y verificado por el agente Khora, 2026-08-13.
