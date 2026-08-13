# Preregistro v2 — Estructura informacional del caminar de la kinesina-1 (F1-B, estadístico agrupado)

**Estado:** `FIRMADO — 2026-08-13. Borrador congelado con SHA
e4aebb2a6511e291e8a988ad8f237aece1643352cdb64229d47363506b756842 ANTES de
ejecutar A2′; A2′ = PASA (resultados/a2_prima/A2P_INFORME.md); firma de
Raúl en sesión ("firmo, ejecuta B1v2") posterior al veredicto de A2′,
conforme a §5. El SHA post-firma se registra en memoria/decisiones.md.`
**Fecha:** 2026-08-13. **Redactor:** Agente Khora. **Gobernador:** Raúl.

## 0. Relación con v1 (pública e íntegra)

- v1 (SHA `c73041ab8bdd81299f27cbdf039af6351b206171d4526bc73ace4a02cc6c71b2`,
  + errata v1.1) permanece inmutable. Su etapa B1 terminó **INCONCLUSO** por
  el criterio de muerte de elegibilidad: 0 trazas ≥ 50 eventos en las 9
  celdas (máximo real: 35; mediana 8; fotoblanqueo MINFLUX). Informe:
  `resultados/b1/B1_INFORME.md`.
- **Qué cambia en v2 y por qué (declaración pública de cambios):**
  1. El estadístico por-trayectoria de v1 pasa a **estadístico agrupado por
     celda** (constructo × [ATP]): conteos de pares consecutivos sumados
     sobre trazas, permutación DENTRO de cada traza. Es exactamente la
     variante que A2 validó (instancias agrupadas).
  2. La elegibilidad "≥ 50 eventos" (inalcanzable) pasa a "**≥ 2 eventos**"
     (mínimo que forma un par) — el criterio menos arbitrario posible.
  3. Se añade la etapa **A2′**: revalidación del estimador en sintético con
     la rejilla de LARGOS REALES antes de tocar contenido informacional.
- **Conocimiento del dato real usado para redactar v2 (declarado):**
  longitudes de traza por celda (`resultados/a2_prima/longitudes_reales.csv`,
  SHA `ee4afd8ee758288b0ef23ed25a6bc7629493c6ad525dc49ce340ed30131957a8`)
  y los agregados de A1 (medianas de paso/dwell). **Ningún estadístico de
  orden temporal ha sido computado sobre datos reales** (B1-v1 murió en la
  puerta de elegibilidad). El contenido que H1 interroga sigue virgen.
- Criterio (i) de validación: se adopta desde el inicio la **lectura
  binomial** (decisión de Raúl sobre A2, 2026-08-12/13): calibración limpia
  = 0 celdas con FPR significativamente > 0.05 (test binomial exacto de una
  cola, corrección BH dentro de cada par estimador × k). Se reporta también
  el conteo literal, como diagnóstico.

## 1. Estimador (congelado, heredado de A2 + B1_CONGELACION)

**Panzeri-Treves, k = 4** para la decisión primaria; k = 8 SIEMPRE
reportada; plug-in, Miller-Madow y NSB registrados como diagnóstico sin
papel decisorio. Discretización y clases de paso idénticas a
`resultados/b1/B1_CONGELACION.md` §2 (cuantiles agrupados por celda, octiles
nativos con k=4 por engrosamiento exacto; |stepx| → {4,8,16} nm por centro
más cercano, bordes 6 y 12 nm; sensibilidad por terciles).

## 2. Etapa A2′ — revalidación en sintético con largos reales (SIN datos reales)

**Rejilla:** las 9 celdas reales, cada una con su **multiconjunto empírico
de longitudes** (trazas ≥ 2 eventos, del CSV congelado) × formas gamma
{1, 2} × mezclas de paso {equilibrada, dominante_8nm, sesgada_4nm} = 54
celdas sintéticas; dwell medio fijo en 25 ms (parámetro de escala puro,
degeneración exacta demostrada en A2). Condiciones: H0 + limping
δ ∈ {0.05, 0.1, 0.2, 0.35, 0.5} + AR(1) ρ ∈ {0.05, 0.1, 0.2, 0.35, 0.5}.
500 réplicas por celda × condición; 1.000 permutaciones; percentil 95.

**Generador:** `a2_estimadores.generar_instancia` REUTILIZADO SIN CAMBIOS,
generando (n_trazas, L_max) y truncando cada traza a su longitud real
(prefijo de proceso estacionario: válido para H0, limping y AR(1)).

**Semillas:** `semilla_de(5, id_sint, id_cond, replica)` para los datos,
`semilla_de(6, id_sint, id_cond, replica)` para las permutaciones, con
id_sint = id_celda_real·6 + id_forma·3 + id_mezcla ∈ [0, 53].

**Criterios de A2′ (congelados):**
- **PASA** si (i) calibración binomial limpia bajo H0 (0 celdas
  significativamente sobre 0.05 tras BH, para PT con k=4 y k=8) Y
  (ii) para PT k=4 existe ≥ 1 celda sintética derivada de cada constructo
  con potencia ≥ 0.80 para alguna magnitud δ ≤ 0.35 **o** ρ ≤ 0.20
  (racional: δ = 0.35 ≈ razón de dwells 2:1, extremo plausible del limping
  publicado sin carga; ρ = 0.20 correlación moderada).
- **NO MEDIBLE** si (i) falla (los estimadores fabrican señal a estos
  largos) o si (ii) falla en TODOS los constructos (la potencia alcanzable
  a los N y largos reales no cubre magnitudes plausibles). En ese caso el
  frente F1-B sobre este dataset SE CIERRA sin tocar los datos y sin
  tercera reformulación (tripwire).
- La **magnitud mínima detectable** por celda se reporta siempre y acompaña
  a cualquier resultado de B1v2 como contexto de sensibilidad.

## 3. Etapa B1v2 — CONFIRMATORIA: memoria más allá de la renovación (agrupada)

**H1 (idéntica a v1):** la secuencia de eventos (dwell, paso) contiene
información predictiva por encima del nulo de renovación en ≥ 1 condición
de ATP, en algún constructo.

- **Estadístico por celda:** I[evento_k ; evento_{k+1}] (PT, k=4) sobre los
  conteos de pares agrupados de TODAS las trazas elegibles (≥ 2 eventos) de
  la celda; nulo primario: 1.000 permutaciones del orden dentro de cada
  traza (marginales por traza preservados). p de celda = p de permutación.
- **Semillas (datos reales):** `semilla_de(7, id_celda_real, 0)` nulo de
  orden; `semilla_de(8, id_celda_real, 0)` nulo de paridad.
- **Decisión:** Holm entre las 3 condiciones DENTRO de cada constructo
  (los constructos se reportan por separado, nunca se mezclan).
  **CONFIRMA** si alguna condición rechaza (I_obs > p95 del nulo) con
  p_Holm < 0.05. **REFUTA** si ninguna condición de ningún constructo
  rechaza. **INCONCLUSO** si alguna celda queda con < 10 trazas elegibles
  (no aplicable según el CSV congelado: mínimo real 52) o si A2′ = NO
  MEDIBLE (en cuyo caso B1v2 no se ejecuta).
- **Guard de limping (innegociable, heredado):** recomputación íntegra con
  el nulo de paridad (permutación solo dentro de posiciones pares/impares
  de cada traza). Si el rechazo desaparece bajo ese nulo → se reporta como
  **replicación de limping** (Asbury 2003, corpus #23), no como hallazgo.
- Sensibilidades obligatorias: k = 8; clases de paso por terciles. Se
  reportan siempre; si cambian el veredicto, el resultado se declara FRÁGIL
  y no asciende.

## 4. B2 y el resto del diseño de v1

B2 (nostalgia_op vs [ATP]) queda como en v1 §3-B2, condicionada a que A2′
valide un estimador de I_mem con estos largos — validación que NO está
incluida en esta A2′ y exigiría su propia adenda pública. Confounders (v1
§6), límites de conclusión (v1 §7) y atribución (v1 §9) se heredan sin
cambios. Gate 4 (verificador externo humano) sigue pendiente y es previo a
B2 y a cualquier ascenso a canon.

## 5. Secuencia de ejecución y firmas

1. Congelar este documento por SHA-256 y registrarlo en
   `memoria/decisiones.md` → **hecho antes de ejecutar A2′**.
2. Ejecutar A2′ (solo sintético; autorizada por Raúl en sesión 2026-08-13).
3. Presentar el veredicto de A2′ a Raúl.
4. **B1v2 sobre datos reales SOLO con firma expresa de Raúl** posterior al
   veredicto de A2′. Sin firma, los datos no se tocan.

## 6. Firma del gobernador

- [x] **FIRMADO** — Raúl, 2026-08-13, en sesión ("firmo, ejecuta B1v2"),
  con el veredicto de A2′ (PASA) y su mapa de sensibilidad a la vista.
  Habilita la ejecución de B1v2 sobre los datos reales.
