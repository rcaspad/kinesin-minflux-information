# Etapa D — Informe de atribución (2026-08-13, noche)

**Congelación:** `D_CONGELACION.md` (SHA `ee0c0b72…27c45615`), previa a la
ejecución. Script `d_analisis.py` (SHA `bbe80ad2…14192a2b63`), semillas
tipo 12. Resultados: `d_resultados.csv`, log completo `d_analisis.log`.
Nota de ejecución: D-4 apuntaba a `N356C/DOL*` por error de ruta; el
control estacionario real está en `stationary_N356/` y se ejecutó ahí
(corrección de ruta, no de regla). El hallazgo colateral de ese error se
registra en §4.

## 1. Resultados por test

**D-1 (revisitas — ¿el sitio tiene identidad?):** T < 0 en las 8 celdas;
significativo en 4: E215C 10µM (T=−0,50, p=0,001), E215C 1mM (−0,43,
p=0,023), T324C 10µM (−0,23, p=0,007), T324C 1mM (−0,36, p=0,001); al
borde en E215C 100µM (p=0,052) y K28C 10µM (p=0,051). **Las pausas en el
mismo sitio se parecen más** — dirección H_raíl, de forma consistente.

**D-2 (decaimiento espacial):** D < 0 fuerte en 7 de 8 celdas (−0,24 a
−0,61, p≈1 en la dirección raíl): **lo contrario de lo predicho por
H_raíl** — la similitud de pausas es MAYOR a través de saltos de 16 nm
que de saltos cortos. Lectura: la similitud no está organizada por
distancia sino por TIPO de paso — pares que flanquean una zancada
completa se parecen (ambos son dwells "de zancada"); pares que flanquean
subpasos alternan. Es la fase del ciclo, otra vez.

**D-3 (rachas):** sin exceso de rachas lentas en 7 de 8 celdas; en varias
las rachas son MÁS CORTAS que el nulo (anti-agrupamiento → alternancia).
**No hay barro**: ni escala espacial de parche ni régimen de caminante.
Única excepción: E215C 100µM (exceso, p=0,002) — justo la celda cuyo
lag-1 era positivo en C. Celda anómala coherente consigo misma.

**D-4 (control estacionario, `stationary_N356/`):** 209 trazas de motor
clavado por la cadena completa → **75 eventos fabricados en total (0,36
por traza)**: la cadena NO fabrica secuencias de eventos sobre ruido real
estacionario (con <2 eventos/traza no hay ni pares). Cota de instrumento
sobre ruido real: despreciable para la estructura secuencial. Complementa
y refuerza C-2.

## 2. VEREDICTO DE ETAPA (regla congelada): ATRIBUCIÓN NO RESUELTA

D-1 apunta a raíl; D-2 y D-3 contradicen la firma de parche y muestran
organización por fase de ciclo; D-4 exonera al instrumento como fabricante
de secuencias. La regla congelada es explícita: tests contradictorios →
NO RESUELTA, y el paquete pasa tal cual al gate 4.

**Síntesis interpretativa (no veredicto):** el cuadro conjunto es
coherente con estructura ligada al CICLO de paso y a retornos al mismo
sitio (dinámica de vaivén de la cabeza marcada), no con parches extensos
del raíl ni con regímenes del caminante. La identidad de sitio de D-1 es
lo más intrigante que queda: o hay sitios con carácter (raíl puntual —
"árboles", no "barro"), o el vaivén fragmenta pausas físicas únicas en
varios eventos al mismo nivel (confounder de fragmentación que ningún
test actual separa). El desempate natural existe en el propio dataset:
las trazas N356C de tallo/centro-de-masas (§4), sin vaivén de cabeza.

## 3. Auditoría de la revisión del Nodo MNS-IC (pinneada contra disco)

1. **"Hueco no señalado" del barajado:** INCORRECTO como cargo — el
   registro del debate ya decía que in vitro sobreviven "defectos
   agrupados, costura o islas de reparación" (DEBATE §2). La expansión
   del Nodo (taxol/GMPCPP como fuente extra) sí añade valor y se acepta.
2. **§3.2 del Nodo ("N356 está en el switch II del dominio motor…"):**
   **FALSO.** Contra el texto del paper (errata v1.1): N356C es una
   cisteína del COILED-COIL (tallo) introducida para MARCAJE
   (center-of-mass), no una mutación cinética del sitio catalítico; todos
   los constructos son variantes de marcaje. El test "WT vs mutante"
   propuesto parte de una premisa inventada. **Se rescata la versión
   correcta**: comparar cabeza-marcada vs tallo-marcada (COM) — el COM no
   tiene vaivén de cabeza, luego separa fragmentación/ciclo de
   identidad-de-sitio. Tercera errata factual del Nodo registrada; la
   regla de pinnear sigue pagándose sola.
3. **§3.3 (polaridad plus/minus):** INAPLICABLE — la kinesina-1 es
   unidireccional (+); no existen caminatas minus-end en el dataset.
4. **§3.1 (química de estabilización con `stabilized_microtubules/`):**
   ACEPTADO como candidato a test, condicionado a verificar contra el
   paper qué contiene ese directorio (sin abrirlo para dwells hasta
   congelar).
5. **§3.4 (protocolo de normalización para Wirth):** ACEPTADO — se
   congelará dentro del prereg de Wirth (longitud de serie, ruido, [ATP],
   temperatura).
6. **§4 (terminología sobria):** ACEPTADO — adoptamos "heterogeneidad del
   sustrato", "desorden congelado (quenched disorder)", "código de
   tubulina"; se retira "el raíl como canal de información" de cualquier
   texto de cara externa.
7. **"La etapa D debe congelarse con estos tests añadidos antes de
   ejecutar":** EXTEMPORÁNEO — D estaba congelada y ejecutada por orden
   de Raúl antes de recibir la revisión. Los tests nuevos aceptados (2 y
   4) usan DATOS NO ABIERTOS (N356C-COM caminando; estabilizados), luego
   caben como **etapa E** con congelación propia, sin violar el "sin
   D-v2" (que protege los cuatro tests de D, no prohíbe evidencia nueva).

## 4. Hallazgo colateral: la errata v1.1 contiene un error factual propio

La errata decía que N356C "en el dataset solo aparece como control
estacionario". FALSO: existen `N356C/DOL1/{10uM,1mM}` y `N356C/DOL2/1mM`
con trazas CAMINANDO y tablas de pasos — datos de centro-de-masas en
marcha. Se emite `prereg_kinesina_ERRATA_v1.2.md`. Consecuencia feliz: el
desempate cabeza-vs-COM (§2-§3.2) es ejecutable con el dataset actual.

Ejecutado y firmado: Agente Khora, 2026-08-13.
