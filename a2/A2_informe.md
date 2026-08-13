# A2 — Validación del estimador de información mutua en datos sintéticos

**Prereg congelado v1.0** · SHA-256 `c73041ab8bdd81299f27cbdf039af6351b206171d4526bc73ace4a02cc6c71b2`  
**Semilla raíz** `20260812` · **Permutaciones** 1000 · **percentil nulo** 95  
**Datos**: exclusivamente sintéticos. Ninguna trayectoria MINFLUX real fue leída, cargada ni consultada en esta etapa.

## 1. Qué se preguntó y qué se responde

A2 pregunta si nuestros estimadores de información mutua **fabrican señal**. La respuesta se decide con dos criterios congelados, no negociables:

- **(i) Control del error de tipo I**: bajo H0 (eventos i.i.d., sin memoria), tasa de falsos positivos ≤ 0.05 **en cada celda** de la rejilla.
- **(ii) Potencia**: detectar la memoria inyectada con potencia ≥ 0.80.

> **Dictamen de esta etapa: NO MEDIBLE bajo el criterio congelado en su lectura literal (alguna celda supera FPR 0.05); MEDIBLE con reservas bajo la lectura binomial. La discrepancia se documenta en §3.1 y la decide el responsable del prereg, no esta etapa**

## 2. Diseño del experimento sintético

**Generador.** Proceso de renovación semi-Markov de eventos de paso: cada evento es el par (dwell time, tamaño de paso). Bajo H0 los pares son i.i.d. y no hay memoria alguna entre eventos. Dos variantes con memoria conocida y parametrizable:

- **limping**: alternancia par/impar de los dwells, factor (1+δ)/(1−δ); la fase se sortea por trayectoria, de modo que la memoria es de lag 1 y no un artefacto de índice absoluto. La media marginal se conserva.
- **AR(1)**: dwells acoplados por una cópula gaussiana de coeficiente ρ; el marginal gamma se conserva exactamente.

Los **tamaños de paso son siempre i.i.d.**: la memoria se inyecta sólo en los dwells. Es deliberado — diluye la señal en el canal «evento» y hace que las potencias reportadas sean una **cota conservadora**.

**Rejilla.** dwell medio [10.0, 25.0, 50.0, 100.0] ms × formas gamma [1.0, 2.0] × 3 mezclas de paso {4, 8, 16} nm × longitudes [50, 100, 200] eventos × [50, 100, 300] trayectorias = **216 celdas**, 500 réplicas independientes por celda y condición.

**Condiciones de memoria.** limping δ ∈ [0.05, 0.1, 0.2, 0.35, 0.5]; AR(1) ρ ∈ [0.05, 0.1, 0.2, 0.35, 0.5].

**Degeneración exacta del eje `dwell_medio_ms`.** El dwell medio es un parámetro de escala puro de la gamma y la discretización por cuantiles es invariante bajo transformaciones estrictamente crecientes. Dos celdas que sólo difieren en el dwell medio producen resultados **idénticos bit a bit** (verificado: máxima diferencia 0.0e+00 sobre 216 comparaciones). Se computaron por tanto las 54 celdas distintas y se expandieron a las 216 del prereg (columna `escala_dwell_degenerada` en los CSV). **Consecuencia sustantiva**: esta rejilla *no* interroga la escala temporal absoluta — un estimador basado en cuantiles no puede, por construcción, distinguir 10 ms de 100 ms.

**Estimadores** (prereg §4). Información mutua I[evento_k ; evento_{k+1}] en bits:

| estimador | descripción |
|---|---|
| plug-in | máxima verosimilitud, sin corrección (referencia de sesgo) |
| Miller-Madow | corrección analítica de primer orden (referencia) |
| **Panzeri-Treves** | plug-in con conteo bayesiano de bins relevantes (prereg §4a) |
| **NSB** | Nemenman-Shafee-Bialek, prior de entropía casi uniforme (prereg §4b) |

**Discretización.** Cuantiles del dwell con k = 4 y k = 8, **ambas siempre reportadas**. Los cuartiles son un engrosamiento exacto de los octiles. El tamaño de paso conserva sus 3 clases nativas {4, 8, 16} nm: discretizar por cuantiles una variable con 3 átomos es degenerado. Alfabeto de evento: 3k símbolos (12 ó 24), alfabeto conjunto del par: (3k)² = 144 ó 576.

**Nulo.** Permutación del orden de los eventos **dentro de cada trayectoria** (1000 permutaciones), que preserva los marginales de cada trayectoria y destruye sólo el orden temporal. Se rechaza si la IM observada supera el percentil 95 del nulo.

## 3. Criterio (i) — tasa de falsos positivos bajo H0

Canal «evento», 216 celdas × 500 réplicas. `fpr_max` es la peor celda de la rejilla; el criterio congelado exige que **ninguna** celda supere 0.05.

| estimador | k | celdas | fpr_media | fpr_p50 | fpr_max | celdas_sobre_0.05 | frac_celdas_>0.05 |
|---|---|---|---|---|---|---|---|
| Miller-Madow | 4 | 216 | 0.048 | 0.050 | 0.074 | 92 | 0.426 |
| NSB | 4 | 216 | 0.049 | 0.049 | 0.074 | 84 | 0.389 |
| Panzeri-Treves | 4 | 216 | 0.049 | 0.050 | 0.074 | 92 | 0.426 |
| plug-in | 4 | 216 | 0.049 | 0.050 | 0.074 | 92 | 0.426 |
| Miller-Madow | 8 | 216 | 0.048 | 0.050 | 0.072 | 96 | 0.444 |
| NSB | 8 | 216 | 0.050 | 0.052 | 0.074 | 112 | 0.519 |
| Panzeri-Treves | 8 | 216 | 0.048 | 0.048 | 0.072 | 92 | 0.426 |
| plug-in | 8 | 216 | 0.048 | 0.048 | 0.072 | 92 | 0.426 |

*(error estándar de una FPR de 0.05 con 500 réplicas: 0.010)*

### 3.1 Lectura binomial — diagnóstico, no sustituto del criterio

El criterio congelado es literal: FPR ≤ 0.05 **en cada celda**. Conviene saber qué exige eso en realidad. Un estimador **perfectamente calibrado** (FPR verdadera exactamente 0.05) supera 0.05 por azar de muestreo en cerca de la mitad de las celdas; con 216 celdas por par (estimador, k) — de las cuales 54 son distintas y el resto réplicas exactas del eje degenerado — la probabilidad de que *ninguna* lo supere es despreciable. La tabla siguiente separa ambas cosas: cuántas celdas superan el umbral (criterio literal) y en cuántas la FPR es **estadísticamente** mayor que 0.05 (test binomial exacto de una cola, corrección de Benjamini-Hochberg dentro de cada par). **El criterio congelado se aplica en su forma literal en la sección 6; esto es diagnóstico añadido, no una relajación.**

| estimador | k | celdas | celdas_sobre_0_05 | celdas_signif_sobre_0_05_BH | fpr_max | fpr_media |
|---|---|---|---|---|---|---|
| Miller-Madow | 4 | 54 | 23 | 0 | 0.074 | 0.048 |
| NSB | 4 | 54 | 21 | 0 | 0.074 | 0.049 |
| Panzeri-Treves | 4 | 54 | 23 | 0 | 0.074 | 0.049 |
| plug-in | 4 | 54 | 23 | 0 | 0.074 | 0.049 |
| Miller-Madow | 8 | 54 | 24 | 0 | 0.072 | 0.048 |
| NSB | 8 | 54 | 28 | 0 | 0.074 | 0.050 |
| Panzeri-Treves | 8 | 54 | 23 | 0 | 0.072 | 0.048 |
| plug-in | 8 | 54 | 23 | 0 | 0.072 | 0.048 |

## 4. Criterio (ii) — potencia frente a memoria inyectada

Potencia media sobre las 216 celdas; `potencia_min` es la peor celda y `celdas_con_potencia_80` la fracción de celdas que alcanzan el umbral congelado.

### 4.1 Limping (alternancia par/impar, magnitud δ)

| magnitud | estimador | k | potencia_media | potencia_min | celdas_con_potencia_80 |
|---|---|---|---|---|---|
| 0.050 | Miller-Madow | 4 | 0.053 | 0.030 | 0.000 |
| 0.100 | Miller-Madow | 4 | 0.078 | 0.042 | 0.000 |
| 0.200 | Miller-Madow | 4 | 0.477 | 0.052 | 0.278 |
| 0.350 | Miller-Madow | 4 | 0.936 | 0.338 | 0.833 |
| 0.500 | Miller-Madow | 4 | 1.000 | 0.996 | 1.000 |
| 0.050 | NSB | 4 | 0.053 | 0.032 | 0.000 |
| 0.100 | NSB | 4 | 0.077 | 0.040 | 0.000 |
| 0.200 | NSB | 4 | 0.474 | 0.050 | 0.278 |
| 0.350 | NSB | 4 | 0.934 | 0.336 | 0.833 |
| 0.500 | NSB | 4 | 1.000 | 0.996 | 1.000 |
| 0.050 | Panzeri-Treves | 4 | 0.053 | 0.034 | 0.000 |
| 0.100 | Panzeri-Treves | 4 | 0.078 | 0.042 | 0.000 |
| 0.200 | Panzeri-Treves | 4 | 0.478 | 0.056 | 0.278 |
| 0.350 | Panzeri-Treves | 4 | 0.937 | 0.376 | 0.833 |
| 0.500 | Panzeri-Treves | 4 | 1.000 | 0.996 | 1.000 |
| 0.050 | plug-in | 4 | 0.053 | 0.034 | 0.000 |
| 0.100 | plug-in | 4 | 0.078 | 0.042 | 0.000 |
| 0.200 | plug-in | 4 | 0.478 | 0.056 | 0.278 |
| 0.350 | plug-in | 4 | 0.938 | 0.382 | 0.833 |
| 0.500 | plug-in | 4 | 1.000 | 0.996 | 1.000 |
| 0.050 | Miller-Madow | 8 | 0.051 | 0.028 | 0.000 |
| 0.100 | Miller-Madow | 8 | 0.066 | 0.028 | 0.000 |
| 0.200 | Miller-Madow | 8 | 0.363 | 0.050 | 0.167 |
| 0.350 | Miller-Madow | 8 | 0.877 | 0.192 | 0.833 |
| 0.500 | Miller-Madow | 8 | 0.996 | 0.922 | 1.000 |
| 0.050 | NSB | 8 | 0.051 | 0.030 | 0.000 |
| 0.100 | NSB | 8 | 0.065 | 0.028 | 0.000 |
| 0.200 | NSB | 8 | 0.358 | 0.054 | 0.167 |
| 0.350 | NSB | 8 | 0.877 | 0.196 | 0.833 |
| 0.500 | NSB | 8 | 0.996 | 0.906 | 1.000 |
| 0.050 | Panzeri-Treves | 8 | 0.052 | 0.026 | 0.000 |
| 0.100 | Panzeri-Treves | 8 | 0.066 | 0.030 | 0.000 |
| 0.200 | Panzeri-Treves | 8 | 0.366 | 0.048 | 0.167 |
| 0.350 | Panzeri-Treves | 8 | 0.884 | 0.188 | 0.833 |
| 0.500 | Panzeri-Treves | 8 | 0.997 | 0.918 | 1.000 |
| 0.050 | plug-in | 8 | 0.052 | 0.026 | 0.000 |
| 0.100 | plug-in | 8 | 0.066 | 0.030 | 0.000 |
| 0.200 | plug-in | 8 | 0.366 | 0.048 | 0.167 |
| 0.350 | plug-in | 8 | 0.884 | 0.188 | 0.833 |
| 0.500 | plug-in | 8 | 0.997 | 0.918 | 1.000 |

### 4.2 AR(1) entre dwells consecutivos (coeficiente ρ)

| magnitud | estimador | k | potencia_media | potencia_min | celdas_con_potencia_80 |
|---|---|---|---|---|---|
| 0.050 | Miller-Madow | 4 | 0.457 | 0.078 | 0.222 |
| 0.100 | Miller-Madow | 4 | 0.844 | 0.262 | 0.667 |
| 0.200 | Miller-Madow | 4 | 0.998 | 0.968 | 1.000 |
| 0.350 | Miller-Madow | 4 | 1.000 | 1.000 | 1.000 |
| 0.500 | Miller-Madow | 4 | 1.000 | 1.000 | 1.000 |
| 0.050 | NSB | 4 | 0.454 | 0.078 | 0.222 |
| 0.100 | NSB | 4 | 0.839 | 0.256 | 0.667 |
| 0.200 | NSB | 4 | 0.998 | 0.954 | 1.000 |
| 0.350 | NSB | 4 | 1.000 | 1.000 | 1.000 |
| 0.500 | NSB | 4 | 1.000 | 1.000 | 1.000 |
| 0.050 | Panzeri-Treves | 4 | 0.458 | 0.082 | 0.222 |
| 0.100 | Panzeri-Treves | 4 | 0.846 | 0.262 | 0.667 |
| 0.200 | Panzeri-Treves | 4 | 0.998 | 0.976 | 1.000 |
| 0.350 | Panzeri-Treves | 4 | 1.000 | 1.000 | 1.000 |
| 0.500 | Panzeri-Treves | 4 | 1.000 | 1.000 | 1.000 |
| 0.050 | plug-in | 4 | 0.458 | 0.082 | 0.222 |
| 0.100 | plug-in | 4 | 0.846 | 0.262 | 0.667 |
| 0.200 | plug-in | 4 | 0.998 | 0.976 | 1.000 |
| 0.350 | plug-in | 4 | 1.000 | 1.000 | 1.000 |
| 0.500 | plug-in | 4 | 1.000 | 1.000 | 1.000 |
| 0.050 | Miller-Madow | 8 | 0.316 | 0.068 | 0.111 |
| 0.100 | Miller-Madow | 8 | 0.712 | 0.134 | 0.463 |
| 0.200 | Miller-Madow | 8 | 0.975 | 0.750 | 0.889 |
| 0.350 | Miller-Madow | 8 | 1.000 | 1.000 | 1.000 |
| 0.500 | Miller-Madow | 8 | 1.000 | 1.000 | 1.000 |
| 0.050 | NSB | 8 | 0.313 | 0.064 | 0.111 |
| 0.100 | NSB | 8 | 0.709 | 0.128 | 0.500 |
| 0.200 | NSB | 8 | 0.976 | 0.764 | 0.889 |
| 0.350 | NSB | 8 | 1.000 | 1.000 | 1.000 |
| 0.500 | NSB | 8 | 1.000 | 1.000 | 1.000 |
| 0.050 | Panzeri-Treves | 8 | 0.318 | 0.066 | 0.111 |
| 0.100 | Panzeri-Treves | 8 | 0.722 | 0.132 | 0.481 |
| 0.200 | Panzeri-Treves | 8 | 0.978 | 0.754 | 0.963 |
| 0.350 | Panzeri-Treves | 8 | 1.000 | 1.000 | 1.000 |
| 0.500 | Panzeri-Treves | 8 | 1.000 | 1.000 | 1.000 |
| 0.050 | plug-in | 8 | 0.318 | 0.066 | 0.111 |
| 0.100 | plug-in | 8 | 0.722 | 0.132 | 0.481 |
| 0.200 | plug-in | 8 | 0.978 | 0.754 | 0.963 |
| 0.350 | plug-in | 8 | 1.000 | 1.000 | 1.000 |
| 0.500 | plug-in | 8 | 1.000 | 1.000 | 1.000 |

## 5. Magnitud mínima detectable (potencia ≥ 0.80)

Por celda, la menor magnitud del barrido que alcanza potencia 0.80 **entre las celdas que además controlan la FPR**. `celdas_sin_deteccion` cuenta las celdas donde ninguna magnitud barrida llegó a 0.80: son el resultado honesto de la resolución del método, no un dato faltante.

| modo | n_tray | largo | estimador | k | mediana_min_detectable | peor_min_detectable | celdas_sin_deteccion | celdas |
|---|---|---|---|---|---|---|---|---|
| ar1 | 50 | 50 | Miller-Madow | 4 | 0.200 | 0.200 | 0 | 16 |
| ar1 | 50 | 100 | Miller-Madow | 4 | 0.200 | 0.200 | 0 | 12 |
| ar1 | 50 | 200 | Miller-Madow | 4 | 0.100 | 0.100 | 0 | 4 |
| ar1 | 100 | 50 | Miller-Madow | 4 | 0.200 | 0.200 | 0 | 20 |
| ar1 | 100 | 100 | Miller-Madow | 4 | 0.100 | 0.100 | 0 | 12 |
| ar1 | 100 | 200 | Miller-Madow | 4 | 0.100 | 0.100 | 0 | 12 |
| ar1 | 300 | 50 | Miller-Madow | 4 | 0.100 | 0.100 | 0 | 12 |
| ar1 | 300 | 100 | Miller-Madow | 4 | 0.050 | 0.050 | 0 | 16 |
| ar1 | 300 | 200 | Miller-Madow | 4 | 0.050 | 0.050 | 0 | 20 |
| ar1 | 50 | 50 | NSB | 4 | 0.200 | 0.200 | 0 | 16 |
| ar1 | 50 | 100 | NSB | 4 | 0.200 | 0.200 | 0 | 12 |
| ar1 | 50 | 200 | NSB | 4 | 0.100 | 0.100 | 0 | 4 |
| ar1 | 100 | 50 | NSB | 4 | 0.200 | 0.200 | 0 | 20 |
| ar1 | 100 | 100 | NSB | 4 | 0.100 | 0.100 | 0 | 12 |
| ar1 | 100 | 200 | NSB | 4 | 0.100 | 0.100 | 0 | 20 |
| ar1 | 300 | 50 | NSB | 4 | 0.100 | 0.100 | 0 | 12 |
| ar1 | 300 | 100 | NSB | 4 | 0.050 | 0.050 | 0 | 16 |
| ar1 | 300 | 200 | NSB | 4 | 0.050 | 0.050 | 0 | 20 |
| ar1 | 50 | 50 | Panzeri-Treves | 4 | 0.200 | 0.200 | 0 | 16 |
| ar1 | 50 | 100 | Panzeri-Treves | 4 | 0.200 | 0.200 | 0 | 12 |
| ar1 | 50 | 200 | Panzeri-Treves | 4 | 0.100 | 0.100 | 0 | 4 |
| ar1 | 100 | 50 | Panzeri-Treves | 4 | 0.200 | 0.200 | 0 | 20 |
| ar1 | 100 | 100 | Panzeri-Treves | 4 | 0.100 | 0.100 | 0 | 12 |
| ar1 | 100 | 200 | Panzeri-Treves | 4 | 0.100 | 0.100 | 0 | 12 |
| ar1 | 300 | 50 | Panzeri-Treves | 4 | 0.100 | 0.100 | 0 | 12 |
| ar1 | 300 | 100 | Panzeri-Treves | 4 | 0.050 | 0.050 | 0 | 16 |
| ar1 | 300 | 200 | Panzeri-Treves | 4 | 0.050 | 0.050 | 0 | 20 |
| ar1 | 50 | 50 | plug-in | 4 | 0.200 | 0.200 | 0 | 16 |
| ar1 | 50 | 100 | plug-in | 4 | 0.200 | 0.200 | 0 | 12 |
| ar1 | 50 | 200 | plug-in | 4 | 0.100 | 0.100 | 0 | 4 |
| ar1 | 100 | 50 | plug-in | 4 | 0.200 | 0.200 | 0 | 20 |
| ar1 | 100 | 100 | plug-in | 4 | 0.100 | 0.100 | 0 | 12 |
| ar1 | 100 | 200 | plug-in | 4 | 0.100 | 0.100 | 0 | 12 |
| ar1 | 300 | 50 | plug-in | 4 | 0.100 | 0.100 | 0 | 12 |
| ar1 | 300 | 100 | plug-in | 4 | 0.050 | 0.050 | 0 | 16 |
| ar1 | 300 | 200 | plug-in | 4 | 0.050 | 0.050 | 0 | 20 |
| ar1 | 50 | 50 | Miller-Madow | 8 | 0.350 | 0.350 | 0 | 8 |
| ar1 | 50 | 100 | Miller-Madow | 8 | 0.200 | 0.200 | 0 | 12 |
| ar1 | 50 | 200 | Miller-Madow | 8 | 0.200 | 0.200 | 0 | 12 |
| ar1 | 100 | 50 | Miller-Madow | 8 | 0.200 | 0.200 | 0 | 20 |
| ar1 | 100 | 100 | Miller-Madow | 8 | 0.200 | 0.200 | 0 | 20 |
| ar1 | 100 | 200 | Miller-Madow | 8 | 0.100 | 0.100 | 0 | 12 |
| ar1 | 300 | 50 | Miller-Madow | 8 | 0.100 | 0.100 | 0 | 8 |
| ar1 | 300 | 100 | Miller-Madow | 8 | 0.100 | 0.100 | 0 | 16 |
| ar1 | 300 | 200 | Miller-Madow | 8 | 0.050 | 0.050 | 0 | 12 |
| ar1 | 50 | 50 | NSB | 8 | 0.350 | 0.350 | 0 | 8 |
| ar1 | 50 | 100 | NSB | 8 | 0.200 | 0.200 | 0 | 12 |
| ar1 | 50 | 200 | NSB | 8 | 0.150 | 0.200 | 0 | 8 |
| ar1 | 100 | 50 | NSB | 8 | 0.200 | 0.200 | 0 | 12 |
| ar1 | 100 | 100 | NSB | 8 | 0.200 | 0.200 | 0 | 20 |
| ar1 | 100 | 200 | NSB | 8 | 0.100 | 0.100 | 0 | 12 |
| ar1 | 300 | 50 | NSB | 8 | 0.100 | 0.100 | 0 | 8 |
| ar1 | 300 | 100 | NSB | 8 | 0.100 | 0.100 | 0 | 16 |
| ar1 | 300 | 200 | NSB | 8 | 0.050 | 0.050 | 0 | 8 |
| ar1 | 50 | 50 | Panzeri-Treves | 8 | 0.350 | 0.350 | 0 | 12 |
| ar1 | 50 | 100 | Panzeri-Treves | 8 | 0.200 | 0.200 | 0 | 12 |
| ar1 | 50 | 200 | Panzeri-Treves | 8 | 0.200 | 0.200 | 0 | 12 |
| ar1 | 100 | 50 | Panzeri-Treves | 8 | 0.200 | 0.200 | 0 | 20 |
| ar1 | 100 | 100 | Panzeri-Treves | 8 | 0.200 | 0.200 | 0 | 20 |
| ar1 | 100 | 200 | Panzeri-Treves | 8 | 0.100 | 0.100 | 0 | 12 |
| ar1 | 300 | 50 | Panzeri-Treves | 8 | 0.100 | 0.100 | 0 | 8 |
| ar1 | 300 | 100 | Panzeri-Treves | 8 | 0.100 | 0.100 | 0 | 16 |
| ar1 | 300 | 200 | Panzeri-Treves | 8 | 0.050 | 0.050 | 0 | 12 |
| ar1 | 50 | 50 | plug-in | 8 | 0.350 | 0.350 | 0 | 12 |
| ar1 | 50 | 100 | plug-in | 8 | 0.200 | 0.200 | 0 | 12 |
| ar1 | 50 | 200 | plug-in | 8 | 0.200 | 0.200 | 0 | 12 |
| ar1 | 100 | 50 | plug-in | 8 | 0.200 | 0.200 | 0 | 20 |
| ar1 | 100 | 100 | plug-in | 8 | 0.200 | 0.200 | 0 | 20 |
| ar1 | 100 | 200 | plug-in | 8 | 0.100 | 0.100 | 0 | 12 |
| ar1 | 300 | 50 | plug-in | 8 | 0.100 | 0.100 | 0 | 8 |
| ar1 | 300 | 100 | plug-in | 8 | 0.100 | 0.100 | 0 | 16 |
| ar1 | 300 | 200 | plug-in | 8 | 0.050 | 0.050 | 0 | 12 |
| limping | 50 | 50 | Miller-Madow | 4 | 0.425 | 0.500 | 0 | 16 |
| limping | 50 | 100 | Miller-Madow | 4 | 0.500 | 0.500 | 0 | 12 |
| limping | 50 | 200 | Miller-Madow | 4 | 0.350 | 0.350 | 0 | 4 |
| limping | 100 | 50 | Miller-Madow | 4 | 0.350 | 0.500 | 0 | 20 |
| limping | 100 | 100 | Miller-Madow | 4 | 0.350 | 0.350 | 0 | 12 |
| limping | 100 | 200 | Miller-Madow | 4 | 0.350 | 0.350 | 0 | 12 |
| limping | 300 | 50 | Miller-Madow | 4 | 0.200 | 0.350 | 0 | 12 |
| limping | 300 | 100 | Miller-Madow | 4 | 0.200 | 0.350 | 0 | 16 |
| limping | 300 | 200 | Miller-Madow | 4 | 0.200 | 0.200 | 0 | 20 |
| limping | 50 | 50 | NSB | 4 | 0.425 | 0.500 | 0 | 16 |
| limping | 50 | 100 | NSB | 4 | 0.500 | 0.500 | 0 | 12 |
| limping | 50 | 200 | NSB | 4 | 0.350 | 0.350 | 0 | 4 |
| limping | 100 | 50 | NSB | 4 | 0.350 | 0.500 | 0 | 20 |
| limping | 100 | 100 | NSB | 4 | 0.350 | 0.350 | 0 | 12 |
| limping | 100 | 200 | NSB | 4 | 0.200 | 0.350 | 0 | 20 |
| limping | 300 | 50 | NSB | 4 | 0.200 | 0.350 | 0 | 12 |
| limping | 300 | 100 | NSB | 4 | 0.200 | 0.350 | 0 | 16 |
| limping | 300 | 200 | NSB | 4 | 0.200 | 0.200 | 0 | 20 |
| limping | 50 | 50 | Panzeri-Treves | 4 | 0.425 | 0.500 | 0 | 16 |
| limping | 50 | 100 | Panzeri-Treves | 4 | 0.500 | 0.500 | 0 | 12 |
| limping | 50 | 200 | Panzeri-Treves | 4 | 0.350 | 0.350 | 0 | 4 |
| limping | 100 | 50 | Panzeri-Treves | 4 | 0.350 | 0.500 | 0 | 20 |
| limping | 100 | 100 | Panzeri-Treves | 4 | 0.350 | 0.350 | 0 | 12 |
| limping | 100 | 200 | Panzeri-Treves | 4 | 0.350 | 0.350 | 0 | 12 |
| limping | 300 | 50 | Panzeri-Treves | 4 | 0.200 | 0.350 | 0 | 12 |
| limping | 300 | 100 | Panzeri-Treves | 4 | 0.200 | 0.350 | 0 | 16 |
| limping | 300 | 200 | Panzeri-Treves | 4 | 0.200 | 0.200 | 0 | 20 |
| limping | 50 | 50 | plug-in | 4 | 0.425 | 0.500 | 0 | 16 |
| limping | 50 | 100 | plug-in | 4 | 0.500 | 0.500 | 0 | 12 |
| limping | 50 | 200 | plug-in | 4 | 0.350 | 0.350 | 0 | 4 |
| limping | 100 | 50 | plug-in | 4 | 0.350 | 0.500 | 0 | 20 |
| limping | 100 | 100 | plug-in | 4 | 0.350 | 0.350 | 0 | 12 |
| limping | 100 | 200 | plug-in | 4 | 0.350 | 0.350 | 0 | 12 |
| limping | 300 | 50 | plug-in | 4 | 0.200 | 0.350 | 0 | 12 |
| limping | 300 | 100 | plug-in | 4 | 0.200 | 0.350 | 0 | 16 |
| limping | 300 | 200 | plug-in | 4 | 0.200 | 0.200 | 0 | 20 |
| limping | 50 | 50 | Miller-Madow | 8 | 0.425 | 0.500 | 0 | 8 |
| limping | 50 | 100 | Miller-Madow | 8 | 0.350 | 0.500 | 0 | 12 |
| limping | 50 | 200 | Miller-Madow | 8 | 0.350 | 0.350 | 0 | 12 |
| limping | 100 | 50 | Miller-Madow | 8 | 0.350 | 0.500 | 0 | 20 |
| limping | 100 | 100 | Miller-Madow | 8 | 0.350 | 0.350 | 0 | 20 |
| limping | 100 | 200 | Miller-Madow | 8 | 0.350 | 0.350 | 0 | 12 |
| limping | 300 | 50 | Miller-Madow | 8 | 0.350 | 0.350 | 0 | 8 |
| limping | 300 | 100 | Miller-Madow | 8 | 0.275 | 0.350 | 0 | 16 |
| limping | 300 | 200 | Miller-Madow | 8 | 0.200 | 0.350 | 0 | 12 |
| limping | 50 | 50 | NSB | 8 | 0.425 | 0.500 | 0 | 8 |
| limping | 50 | 100 | NSB | 8 | 0.350 | 0.500 | 0 | 12 |
| limping | 50 | 200 | NSB | 8 | 0.350 | 0.350 | 0 | 8 |
| limping | 100 | 50 | NSB | 8 | 0.350 | 0.500 | 0 | 12 |
| limping | 100 | 100 | NSB | 8 | 0.350 | 0.350 | 0 | 20 |
| limping | 100 | 200 | NSB | 8 | 0.350 | 0.350 | 0 | 12 |
| limping | 300 | 50 | NSB | 8 | 0.350 | 0.350 | 0 | 8 |
| limping | 300 | 100 | NSB | 8 | 0.275 | 0.350 | 0 | 16 |
| limping | 300 | 200 | NSB | 8 | 0.275 | 0.350 | 0 | 8 |
| limping | 50 | 50 | Panzeri-Treves | 8 | 0.350 | 0.500 | 0 | 12 |
| limping | 50 | 100 | Panzeri-Treves | 8 | 0.350 | 0.500 | 0 | 12 |
| limping | 50 | 200 | Panzeri-Treves | 8 | 0.350 | 0.350 | 0 | 12 |
| limping | 100 | 50 | Panzeri-Treves | 8 | 0.350 | 0.500 | 0 | 20 |
| limping | 100 | 100 | Panzeri-Treves | 8 | 0.350 | 0.350 | 0 | 20 |
| limping | 100 | 200 | Panzeri-Treves | 8 | 0.350 | 0.350 | 0 | 12 |
| limping | 300 | 50 | Panzeri-Treves | 8 | 0.350 | 0.350 | 0 | 8 |
| limping | 300 | 100 | Panzeri-Treves | 8 | 0.275 | 0.350 | 0 | 16 |
| limping | 300 | 200 | Panzeri-Treves | 8 | 0.200 | 0.350 | 0 | 12 |
| limping | 50 | 50 | plug-in | 8 | 0.350 | 0.500 | 0 | 12 |
| limping | 50 | 100 | plug-in | 8 | 0.350 | 0.500 | 0 | 12 |
| limping | 50 | 200 | plug-in | 8 | 0.350 | 0.350 | 0 | 12 |
| limping | 100 | 50 | plug-in | 8 | 0.350 | 0.500 | 0 | 20 |
| limping | 100 | 100 | plug-in | 8 | 0.350 | 0.350 | 0 | 20 |
| limping | 100 | 200 | plug-in | 8 | 0.350 | 0.350 | 0 | 12 |
| limping | 300 | 50 | plug-in | 8 | 0.350 | 0.350 | 0 | 8 |
| limping | 300 | 100 | plug-in | 8 | 0.275 | 0.350 | 0 | 16 |
| limping | 300 | 200 | plug-in | 8 | 0.200 | 0.350 | 0 | 12 |

## 6. Veredicto por estimador y k

| estimador | k | fpr_max | criterio_i_fpr | frac_celdas_con_potencia80_delta020 | frac_celdas_con_potencia80_rho020 | celdas_sin_deteccion_en_rejilla | celdas_evaluadas | criterio_i_fpr_binomial |
|---|---|---|---|---|---|---|---|---|
| Miller-Madow | 4 | 0.074 | False | 0.093 | 0.630 | 0 | 248 | True |
| NSB | 4 | 0.074 | False | 0.093 | 0.630 | 0 | 264 | True |
| Panzeri-Treves | 4 | 0.074 | False | 0.093 | 0.630 | 0 | 248 | True |
| plug-in | 4 | 0.074 | False | 0.093 | 0.630 | 0 | 248 | True |
| Miller-Madow | 8 | 0.072 | False | 0.056 | 0.488 | 0 | 240 | True |
| NSB | 8 | 0.074 | False | 0.056 | 0.500 | 0 | 208 | True |
| Panzeri-Treves | 8 | 0.072 | False | 0.056 | 0.519 | 0 | 248 | True |
| plug-in | 8 | 0.072 | False | 0.056 | 0.519 | 0 | 248 | True |

## 7. Diagnóstico: canales secundarios

El prereg fija el canal «evento» como primario. Se registran además, por marginalización exacta de los mismos conteos, los canales «dwell» (I[dwell_k ; dwell_{k+1}]), «paso» y «dwell_paso». El canal «paso» es un control interno: **por construcción los tamaños de paso son i.i.d. en toda condición**, así que cualquier rechazo ahí es un falso positivo, incluso bajo memoria inyectada.

**FPR bajo H0, canal «dwell»:**

| estimador | k | celdas | fpr_media | fpr_p50 | fpr_max | celdas_sobre_0.05 | frac_celdas_>0.05 |
|---|---|---|---|---|---|---|---|
| Miller-Madow | 4 | 216 | 0.050 | 0.049 | 0.082 | 92 | 0.426 |
| NSB | 4 | 216 | 0.050 | 0.050 | 0.082 | 88 | 0.407 |
| Panzeri-Treves | 4 | 216 | 0.050 | 0.049 | 0.082 | 92 | 0.426 |
| plug-in | 4 | 216 | 0.050 | 0.049 | 0.082 | 92 | 0.426 |
| Miller-Madow | 8 | 216 | 0.053 | 0.053 | 0.084 | 116 | 0.537 |
| NSB | 8 | 216 | 0.052 | 0.052 | 0.084 | 120 | 0.556 |
| Panzeri-Treves | 8 | 216 | 0.053 | 0.053 | 0.084 | 116 | 0.537 |
| plug-in | 8 | 216 | 0.053 | 0.053 | 0.084 | 116 | 0.537 |

**FPR bajo H0, canal «paso» (control interno):**

| estimador | k | celdas | fpr_media | fpr_p50 | fpr_max | celdas_sobre_0.05 | frac_celdas_>0.05 |
|---|---|---|---|---|---|---|---|
| Miller-Madow | 4 | 216 | 0.051 | 0.050 | 0.074 | 104 | 0.481 |
| NSB | 4 | 216 | 0.051 | 0.051 | 0.074 | 108 | 0.500 |
| Panzeri-Treves | 4 | 216 | 0.051 | 0.050 | 0.074 | 104 | 0.481 |
| plug-in | 4 | 216 | 0.051 | 0.050 | 0.074 | 104 | 0.481 |
| Miller-Madow | 8 | 216 | 0.051 | 0.050 | 0.074 | 104 | 0.481 |
| NSB | 8 | 216 | 0.051 | 0.051 | 0.074 | 108 | 0.500 |
| Panzeri-Treves | 8 | 216 | 0.051 | 0.050 | 0.074 | 104 | 0.481 |
| plug-in | 8 | 216 | 0.051 | 0.050 | 0.074 | 104 | 0.481 |

**Control interno bajo memoria inyectada — la prueba directa de que no se fabrica señal.** En las condiciones con memoria, los dwells SÍ están correlacionados pero los tamaños de paso siguen siendo i.i.d. Si los estimadores inventaran estructura, o si la memoria de los dwells se filtrara al canal equivocado, la tasa de rechazo en «paso» subiría por encima de 0.05. Tasa de rechazo observada en el canal «paso», agregando **todas** las condiciones con memoria inyectada:

| estimador | k | celdas | tasa_media | tasa_max |
|---|---|---|---|---|
| Miller-Madow | 4 | 2160 | 0.0501 | 0.0840 |
| NSB | 4 | 2160 | 0.0501 | 0.0820 |
| Panzeri-Treves | 4 | 2160 | 0.0501 | 0.0840 |
| plug-in | 4 | 2160 | 0.0501 | 0.0840 |
| Miller-Madow | 8 | 2160 | 0.0501 | 0.0840 |
| NSB | 8 | 2160 | 0.0501 | 0.0820 |
| Panzeri-Treves | 8 | 2160 | 0.0501 | 0.0840 |
| plug-in | 8 | 2160 | 0.0501 | 0.0840 |

## 8. Recomendación para la etapa B1

Regla de decisión declarada: entre los pares (estimador, k) **del prereg §4** que controlan la FPR, mayor potencia en el régimen realista (magnitud ≤ 0.20); empates por menor magnitud mínima detectable mediana y, después, por menor FPR máxima. plug-in y Miller-Madow figuran en el barrido como referencias de sesgo y no son candidatos.

- **Empate exacto en los tres criterios entre: NSB (k=4), Panzeri-Treves (k=4).** Ninguna de las cifras del barrido los separa, así que A2 no elige entre ellos: la decisión es del responsable del prereg. Si hace falta un criterio adicional, Panzeri-Treves es computacionalmente más barato y NSB tiene menos sesgo residual con alfabetos grandes — ninguno de los dos efectos es visible en esta rejilla.
- FPR máxima sobre la rejilla: 0.074 (criterio (i): NO CUMPLE)
- Fracción de celdas con potencia ≥ 0.80: δ ≤ 0.20 → 0.09; ρ ≤ 0.20 → 0.63
- Celdas de la rejilla sin ninguna detección: 0 de 264

**Ambas k se reportan siempre en B1** (requisito del prereg); la recomendación fija cuál se usa para la decisión primaria.

## 9. Limitaciones honestas

1. **El eje de escala temporal no está interrogado.** La discretización por cuantiles hace que dwell medio 10 ms y 100 ms sean el mismo problema, exactamente. A2 no dice nada sobre sensibilidad a la escala absoluta, y B1 no debe interpretar estos resultados como cobertura de ese eje.
2. **La memoria inyectada vive sólo en los dwells.** Los tamaños de paso son i.i.d. en todas las condiciones. Una memoria real que acoplara dwell y tamaño de paso sería un régimen no cubierto aquí.
3. **Las formas de memoria probadas son dos, y ambas de lag 1.** limping alterna y AR(1) decae; una memoria de mayor alcance, no estacionaria o con estructura de estados no se ha ensayado. La potencia frente a esas alternativas es desconocida, no «buena por extensión».
4. **La rejilla no está calibrada con datos reales** — por diseño y por mandato del prereg. Los rangos son plausibles según la literatura, pero si el régimen experimental real cae fuera, estas garantías no se transfieren.
5. **El test es de una cola** (IM observada > percentil 95 del nulo), coherente con que la información mutua sea no negativa; no detecta «menos estructura de la esperada».
6. **Las réplicas por celda son finitas**: con 500 réplicas, una FPR verdadera de 0.05 se estima con error estándar 0.010. Una celda con FPR observada de 0.06 no es distinguible de una de 0.05 con esta resolución; el criterio se aplicó tal como está congelado, sin margen de tolerancia.
7. **Los cuatro estimadores son casi indistinguibles en esta rejilla.** La corrección de sesgo (Panzeri-Treves, NSB) desplaza el valor de la IM pero apenas cambia la decisión, porque el nulo por permutación se recalcula con el MISMO estimador: un sesgo común al observado y al nulo se cancela en la comparación. Esto no dice que las correcciones sean inútiles — dice que en un contraste basado en permutación su aporte es marginal, y que la elección entre ellas no es el factor limitante. Lo limitante es el tamaño de muestra.
8. **Los estimadores comparten conteos y nulo.** plug-in, Miller-Madow, Panzeri-Treves y NSB se evalúan sobre exactamente los mismos datos y las mismas permutaciones. Sus resultados están correlacionados y no constituyen confirmaciones independientes entre sí (regla de eco correlacionado).

## 10. Reproducción

```bash
python a2_pruebas.py       # 13 comprobaciones del código (deben pasar)
python a2_barrido.py 9      # barrido completo
python a2_figuras.py       # figuras
python a2_informe.py       # este informe
```

`a2_pruebas.py` verifica las propiedades de las que dependen estos resultados y es la primera cosa que debe ejecutarse al reproducir: determinismo bit a bit; reanudabilidad del barrido (un bloque presente se omite sin reescribirse y su recálculo reproduce el fichero byte a byte); marginales y correlación de lag 1 del generador en las tres condiciones; que los tamaños de paso sean i.i.d. en toda condición (premisa del control interno de §7); la degeneración exacta del eje de dwell; que k = 4 sea engrosamiento exacto de k = 8; los estimadores de entropía contra la uniforme y el orden del sesgo bajo submuestreo; y la convergencia de la cuadratura NSB.

Toda instancia queda determinada por `(20260812, tipo, id_celda, id_cond, replica)` vía `numpy.random.SeedSequence` (`a2_estimadores.semilla_de`): tipo 1 genera los datos, tipo 2 las permutaciones. La ejecución es reproducible bit a bit y no depende del número de procesos.

## 11. Archivos y SHA-256

`a2_por_instancia.csv.gz` es el reporte por-instancia canónico: las 216 celdas de la rejilla del prereg (38.0 M filas). `a2_por_instancia_54celdas.parquet` contiene la misma información sin las réplicas exactas del eje degenerado (54 celdas, 9.5 M filas, ~5x menor): es el fichero práctico para reanalizar, y del que las 216 celdas se recuperan con `a2_barrido.expandir_a_rejilla_completa`.

| archivo | bytes | SHA-256 |
|---|---|---|
| `a2_estimadores.py` | 18,283 | `dd784b2c7ce7156dfa318a7317b816c0cda8c51d571fee111a1e357d6eab9936` |
| `a2_barrido.py` | 12,136 | `94b0b75236160ccdccbadfddd00a69f7a126cf894c38f125ffe5e9756d015e36` |
| `a2_figuras.py` | 6,878 | `5205d8bb84dceaf748708c0cf89e29ba7cdc498f45ed0bfd0f0631d2c9f861f9` |
| `a2_informe.py` | 29,148 | `05e2c8952eb1093bba80a569022b13888315e32d62644e05f721b733c455d7db` |
| `a2_pruebas.py` | 6,851 | `3590ea591fe03a793e9f6fa219762128cb740648e1ee8516d68729d86d319594` |
| `a2_por_instancia.csv.gz` | 548,744,455 | `38d36f49364a4da2adb2250d0cdc7041883d10f0746640aeb6ea1f88713f0839` |
| `a2_por_instancia_54celdas.parquet` | 114,491,231 | `0012b4da5062a40a3f0ee29c271822e596c9830eae111bfef36906bcd12f37c0` |
| `a2_por_celda.csv` | 10,576,284 | `e665802bd73af526dda4af7e0eb6e1f4b654f038b2ee07bed7edac767ef23fe9` |
| `a2_minima_detectable.csv` | 454,880 | `1beb7bb387715cf7d2fa659f460a59963c685113545589a0f05da14688b9d3ee` |
| `a2_tabla_fpr.csv` | 451 | `786fa596bf9460457aec338fd9a4feb5dd8d4678c1862aa6ddaa8137a64178c0` |
| `a2_tabla_potencia.csv` | 3,098 | `6e2569ec2a65ce82761644e91bd4ef82a60711605225884101da0e2ac27a24e2` |
| `a2_tabla_minima_detectable.csv` | 5,737 | `b45b6270ea9ef9506e4173e15d5a397e5972db97ad28270f6f146ecd20126bfc` |
| `a2_veredicto.csv` | 575 | `ca86884094457206a2270ef1fab1e14b09b0a94ebda0f15e81d906b68a47e50d` |
| `a2_fig_fpr.png` | 168,472 | `bcf4605321fde1943cf300536499b82254e12e2618871a50881ab89d2b296072` |
| `a2_fig_potencia.png` | 269,189 | `9a210e076f01e67db152d3bbe5462a1399c1744cdd55a78dd7519976e1b1802d` |
