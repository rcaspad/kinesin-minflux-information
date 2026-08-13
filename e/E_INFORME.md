# Etapa E — Informe del desempate (2026-08-13, noche)

**Congelación:** `E_CONGELACION.md` (SHA `187004f1…10762f5f`, con E-1c
añadido ANTES de ejecutar). Script `e_analisis.py` (SHA
`54a758ed…87df039b`), semillas tipo 13-16. Resultados:
`e1_resultados.csv`, `e2_resultados.csv`, log `e_analisis.log`.
Terminología: `../../TERMINOLOGIA.md`.

## 1. Resultados

**E-1 — trazas de centro de masas (N356C), 3/3 celdas comparables
(ratios 0,78 / 0,92 / 1,01):**

| test | DOL1 10µM | DOL1 1mM | DOL2 1mM |
|---|---|---|---|
| E-1a residuo dwell→dwell (IM, bits) | 0,0031 (p 0,55) | 0,0056 (p 0,46) | 0,0033 (p 0,37) |
| E-1b identidad de sitio (T, p) | −0,22 (0,30) | −0,26 (0,54) | +0,20 (1,0) |
| E-1c periodicidad de dos pasos (D, p) | −0,007 (0,95) | −0,063 (0,48) | −0,004 (0,96) |

**Los tres efectos son NULOS en trazas de centro de masas.** El residuo
dwell→dwell, la identidad de sitio y la organización por tipo de paso
desaparecen cuando el marcaje no está en la cabeza.

**E-2 — test cruzado de identidad de sitio entre motores ("prueba
reina"; 19 microtúbulos, muestras activamente estabilizadas, 2.566 pares
cruzados mismo-sitio):** ρ = **−0,019** (p = 0,77); estrato mismo lado
(n = 1.813): ρ = −0,040. **Motores distintos sobre el mismo microtúbulo
NO comparten sitios lentos.** Con 2.566 pares, la potencia para detectar
una identidad de sitio compartida modesta era sobrada.

## 2. VEREDICTO: ATRIBUCIÓN RESUELTA

Aplicando las predicciones congeladas:

- E-1a nula Y E-1b nula → predicción de **dinámica de la cabeza marcada**
  (congelación E-1: "el COM no parpadea").
- E-2 nula → descarta **heterogeneidad del sustrato** con identidad de
  sitio reproducible (la vía que D-1 dejaba abierta).
- E-1a nula descarta también **memoria del motor** de ciclo químico (su
  predicción era que el residuo persistiera en centro de masas).
- El **artefacto de instrumento** ya estaba acotado por C-2 y D-4.

**El residuo dwell→dwell de B1v2/C queda atribuido a la fuente 1 de la
taxonomía: el ciclo de paso, en su manifestación de dinámica de la cabeza
marcada** (subpasos, retrocesos y fragmentación de pausas del fluoróforo
en cabeza). No hay evidencia de memoria del motor más allá del ciclo, ni
de heterogeneidad del sustrato con efecto medible, en este dataset.

Cierre del arco B1v2→C→D→E: la información confirmada en B1v2 es
enteramente atribuible al ciclo de paso (incluida su dinámica de cabeza);
la cuantificación (0,06-0,51 bits/evento) y el nulo de artefacto
calibrado quedan como las dos aportaciones del frente.

## 3. Lo que este veredicto NO dice

1. No dice que el microtúbulo sea homogéneo — dice que su heterogeneidad
   no deja identidad de sitio medible en ESTE ensayo (in vitro,
   reconstituido, sin MAPs, N de este dataset).
2. No dice que la memoria del motor no exista — dice que no es medible
   por encima del ciclo de paso a esta resolución y N. Las predicciones
   direccionales para el dataset in vivo (Wirth) del debate quedan
   vigentes.
3. Los cuatro tests de D quedan reinterpretados coherentemente: la
   identidad de sitio de D-1 era intra-traza (la misma cabeza volviendo
   al mismo sitio en la misma pasada) y desaparece entre motores (E-2) y
   sin cabeza marcada (E-1b): era la cabeza, no el sitio.

## 4. Consecuencias

- La carta del gate 4 se actualiza (v1.2): el "residuo sin explicar" pasa
  a "residuo atribuido a dinámica de la cabeza marcada mediante control
  de centro de masas y test cruzado entre motores"; las preguntas al
  verificador pasan a auditar la cadena de atribución completa.
- El frente kinesina in vitro queda **cerrado en resultados**: pendientes
  solo gate 4 (envío de Raúl) y la decisión sobre Wirth (in vivo).
- Sin etapa F sobre este dataset (regla de E: se acabaron las rendijas).

Ejecutado y firmado: Agente Khora, 2026-08-13.
