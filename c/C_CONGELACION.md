# Etapa C — Sensibilidad sobre crudo: congelación previa (2026-08-13)

**Mandato:** Raúl, en sesión ("adelante con la deteccion de pasos
independiente sobre el crudo"), ejecutando la vía de ascenso declarada en
`../b1v2/B1V2_INFORME.md` §4.1 y contemplada en prereg v1 §6 (sensibilidad
sobre crudo). Prereg v2 FIRMADO, SHA post-firma `aa461a6e…36f22e6f`.

**Pregunta única:** el componente dwell→dwell de B1v2 (7/9 celdas,
anticorrelación lag-1), ¿es del motor o del segmentador?

## 1. Segmentador independiente (declarado antes de ejecutar)

- **Algoritmo:** Kalafut-Visscher 2008 — colocación voraz de puntos de
  cambio aceptados solo si reducen el BIC gaussiano por tramos; **sin
  parámetros de usuario** (sin ventana de suavizado, sin umbral de tamaño
  de paso, sin fusión posterior). Familia y estructura de error DISTINTAS
  del pipeline de los autores (`ischange` con umbral 12σ² + mediana móvil 9
  + fusión < 5 nm).
- **Entrada:** posiciones crudas FPGA (columnas x, y de los .txt), NO las
  posiciones SCE refinadas de los autores — independencia también en la
  etapa de localización. Filtro de validez mínimo y declarado: L_x > 0 y
  fotones totales por localización en [7, 150] (constantes tomadas de los
  autores; sin mediana móvil). Proyección al eje principal por regresión
  (mismo criterio de alineamiento que los autores: es geometría, no
  segmentación). Tiempo: suma acumulada de dt (col. 18).
- **Eventos:** dwell_j = duración del tramo j (solo tramos INTERIORES; el
  primero y el último son censurados y se excluyen); paso_j = |salto de
  nivel| al final del tramo j. Secuencia (dwell, paso) por traza, mismas
  clases {4,8,16} nm (bordes 6/12) y cuantiles agrupados por celda que B1v2.

## 2. Gate de comparabilidad (previo a cualquier veredicto)

Si el segmentador KV recupera, por celda, menos del 50 % o más del 200 %
de los pasos de las tablas de los autores, la celda se declara NO
COMPARABLE y no participa; si quedan < 5 celdas comparables, la etapa
entera es INCONCLUSA. (El KV sin umbral de tamaño detectará más pasos
pequeños; el gate acota que siga siendo el mismo fenómeno.)

## 3. Los dos contrastes (congelados)

**C-1 — Replicación con segmentador independiente:** el test de B1v2
(canal dwell→dwell, PT k=4, 1.000 permutaciones dentro de traza, nulos de
orden y paridad) sobre las secuencias KV. Semillas:
`semilla_de(9, id_celda_real, 0)` orden; `semilla_de(10, id_celda_real, 0)`
paridad. Holm entre condiciones por constructo, como B1v2.

**C-2 — Nulo de artefacto calibrado (el discriminador):** por celda, 200
réplicas sustitutas SIN memoria: trayectorias de renovación (dwells i.i.d.
remuestreados del marginal empírico KV de la celda; pasos i.i.d. con signo
del marginal empírico; camino constante a tramos) muestreadas en la MISMA
rejilla temporal real de cada traza, con ruido gaussiano de desviación
σ estimada por traza (std(diff)/√2, el estimador de los autores), pasadas
por el MISMO segmentador KV. Semillas: `semilla_de(11, id_celda_real,
replica)`. En cada sustituta se mide la IM dwell→dwell agrupada (PT k=4)
del resultado segmentado: es la "memoria" que fabrica la cadena
localización+segmentación sola.

## 4. Regla de decisión (congelada)

Por celda comparable que confirmó dwell→dwell en B1v2:
- **ARTEFACTO** si IM_obs(KV real) ≤ percentil 95 del nulo de artefacto —
  la segmentación sola fabrica tanta o más señal.
- **CANDIDATO SOBREVIVE** si IM_obs(KV real) > percentil 95 del nulo de
  artefacto Y C-1 confirma la celda (Holm < 0,05) con el mismo signo de
  correlación lag-1.
- Veredicto de etapa: el candidato dwell→dwell **PERSISTE** si sobrevive en
  ≥ la mitad de las celdas comparables que confirmaban; **ARTEFACTO
  (total o parcial)** en caso contrario, celda a celda.
- Se registra pase lo que pase: n_pares KV vs autores (comparabilidad del
  sesgo del estimador), signo y magnitud de la correlación lag-1 en real
  KV, real autores y sustitutas.

**Lo que esta etapa NO puede hacer:** probar positivamente el origen motor
(dos segmentadores sobre el mismo ruido comparten parte del error de
fronteras; el nulo de artefacto lo modela, pero con supuestos — renovación
exacta, ruido gaussiano, σ constante por traza). El veredicto externo
(gate 4) sigue siendo obligatorio antes de canonizar.

Congelado por el agente Khora, 2026-08-13, antes de segmentar el crudo.
