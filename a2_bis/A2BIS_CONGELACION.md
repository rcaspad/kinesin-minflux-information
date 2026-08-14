# A2″ — Congelación de implementación (2026-08-14, antes de ejecutar)

Criterios del gate: `prereg_wirth_v1.md` §5 (SHA `1d6ab5d5…2586a2de`),
inmutables. Este documento fija las elecciones de implementación que el
prereg delega, ANTES de la primera ejecución.

## Celdas e ids de semilla (raíz 20260812, `semilla_de`)

| celda | id | celda | id |
|---|---|---|---|
| Wolff K28C 10 µM | 3 | Wirth K28C 50 µM | 20 |
| Wolff K28C 100 µM | 4 | Wirth K28C 500 µM | 21 |
| Wolff K28C 1 mM | 5 | Wirth K28C 5 mM | 22 |
| Wolff T324C 10 µM | 6 | Wirth T324C 50 µM | 23 |
| Wolff T324C 1 mM | 8 | Wirth T324C 5 mM | 24 |

(Ids de Wolff heredados de C por trazabilidad.) Datos sintéticos:
tipo 17 (id, réplica). Permutaciones: tipo 18 (id, réplica).

## Marginales usadas (declaración del prereg §5, concretada)

1. **Rejillas temporales y σ por traza:** de las trazas crudas reales de
   cada celda, con los filtros congelados de la cadena (fotones 7/150,
   Lx>0, ≥20 localizaciones válidas); σ = std(diff(z))/√2 sobre el z real
   proyectado (construcción de C-2, precedente citado).
2. **Pools de dwell (tau) y paso (stepx, CON signo):** de las tablas
   `allsteps_reeval.xls` de los AUTORES de cada celda (marginales
   publicadas). Motivo: evita pasar nuestra segmentación por contenido
   real antes de la puerta; las tablas de los autores son insumo externo
   ya publicado y solo se usan sus marginales (sorteos i.i.d.).

## Generación (renovación, verdad conocida)

Por traza real (t, σ): construcción de `sustituta_de_traza` de C-2
(sorteo i.i.d. de dwells y pasos hasta cubrir la duración, niveles
acumulados, ruido gaussiano σ en la rejilla real), extendida para
registrar los tiempos de cambio y la **clase verdadera** de cada paso
sorteado: clase = searchsorted([6,12], |paso|) ∈ {4,8,16} nm. La
renovación i.i.d. garantiza ausencia de estructura secuencial: todo
positivo del test es falso positivo por construcción.

## A2″-1 (FPR, por celda, las 10)

Por réplica: cadena completa (KV+fusión+eventos, scripts congelados por
SHA en prereg §3) → símbolos 24-arios (octiles de dwell agrupados de la
réplica × 3 clases de paso) → IM dwell→dwell PT k=4 agrupada → test de
permutación de orden dentro de traza, **1.000 permutaciones** (tipo 18).
FP = `rechaza` a α = 0,05. Con n = 500 réplicas: banda binomial EXACTA
del 95 % → **[15, 36] falsos positivos** (binom.ppf(0,025/0,975; 500;
0,05) = 15/36). Cada celda debe caer dentro.

## A2″-2 (resolubilidad, gate solo en las 5 celdas de Wirth; Wolff se reporta como referencia)

Por réplica y traza: cada transición detectada por la cadena (frontera
fusionada interior, tiempo = t[borde]) se empareja con la transición
verdadera MÁS CERCANA en tiempo; acierto si clase(|paso estimado|) ==
clase verdadera del paso emparejado. Matriz de confusión acumulada por
celda sobre las 500 réplicas. Gate: exactitud global ≥ 80 % Y exactitud
por clase verdadera ≥ 2/3 para cada una de las tres clases, en las 5
celdas de Wirth. (La detección anclada en lo detectado penaliza las
detecciones espurias por sí sola; los pasos verdaderos no detectados
—p. ej. 4 nm fusionados— aparecen en el diagnóstico de cobertura, que se
reporta pero no es criterio del gate.)

## Réplicas y salida

500 réplicas por celda (prereg §5). Salidas: `a2bis_fpr.csv` (por celda:
FP observados, banda, veredicto), `a2bis_confusion.csv` (matriz por
celda), `a2bis_cobertura.csv` (diagnóstico), log `a2bis.log`,
`A2BIS_INFORME.md`. Script: `a2bis_analisis.py` (SHA en el informe).
Paralelización por (celda, bloque de réplicas); el orden de semillas es
por réplica, así que la paralelización no afecta a la reproducibilidad.

Congelado por el Agente Khora antes de ejecutar. 2026-08-14.
