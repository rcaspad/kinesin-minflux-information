# Etapa C — Informe: ¿motor o segmentador? (2026-08-13)

**Congelación:** `C_CONGELACION.md` + `C_CONGELACION_v2.md` (comprometidas
antes de los resultados, commit 9ef8fd4a). Scripts: `c_kv_pipeline.py`
(SHA `be3756b8…8c4d84d6`), `c_analisis.py` (SHA `0efdfb50…53081b44`).
Semillas tipo 9/10/11, raíz 20260812. Salidas: `c1_replicacion.csv`,
`c2_nulo_artefacto.csv`.

## 1. El hallazgo central: la firma de anticorrelación ES del segmentador

El nulo de artefacto (200 sustitutas de renovación SIN memoria, con el
ruido σ y la rejilla temporal reales, por la misma cadena KV) demuestra
**cuantitativamente el mecanismo sospechado**: las sustitutas sin memoria
salen del segmentador con correlación lag-1 mediana de **−0,075 a −0,107**
— el mismo orden que los datos reales (−0,08 a −0,19). **La anticorrelación
de dwells adyacentes, como firma, queda explicada por el corte de
fronteras** y no puede volver a usarse como evidencia de memoria.

## 2. Pero la información mutua de 3 celdas EXCEDE lo que el artefacto fabrica

La IM dwell→dwell que las sustitutas fabrican es minúscula (p95 entre
0,0013 y 0,0092 bits). En tres celdas — **una por constructo** — la IM real
KV la supera con claridad y además replica con el segmentador independiente
(C-1, Holm) con el mismo signo:

| celda (B1v2 confirmaba) | C-1 KV replica | C-2 sobre nulo artefacto | signo | veredicto celda |
|---|---|---|---|---|
| E215C 10µM | ✅ p_Holm 0,003 | ✅ IM 0,019 vs p95 0,002 (p 0,005) | −0,13 ✓ | **SOBREVIVE** |
| E215C 100µM | ✅ p_Holm 0,003 | ✅ IM 0,037 vs p95 0,007 (p 0,005) | **+0,03 ✗** | NO (regla de signo) |
| K28C 10µM | ✅ p_Holm 0,015 | ✅ IM 0,006 vs p95 0,003 (p 0,005) | −0,10 ✓ | **SOBREVIVE** |
| K28C 100µM | ✗ p 0,34 | ✗ p_artefacto 0,40 | — | **ARTEFACTO** |
| K28C 1mM | ✗ p 0,25 | ✗ p_artefacto 0,065 | −0,08 ≈ surr −0,075 | **ARTEFACTO** |
| T324C 10µM | ✗ p 0,99 | ✅ (justo: 0,0015 vs 0,0013) | — | NO REPLICA con KV |
| T324C 1mM | ✅ p_Holm 0,002 | ✅ IM 0,0137 vs p95 0,0013 (p 0,005) | −0,16 ✓ | **SOBREVIVE** |

(T324C 100µM excluida por el gate de comparabilidad; E215C 1mM no
confirmaba en B1v2 — con KV sí da señal, p_Holm 0,031 y artefacto
descartado: se registra como observación, sin estatuto confirmatorio.)

## 3. VEREDICTO DE ETAPA (regla congelada §4, aplicada sin retoque)

Sobreviven **3 de 7** celdas evaluables; la regla exigía ≥ la mitad (≥ 4)
para declarar PERSISTE. Por tanto: **ARTEFACTO PARCIAL, celda a celda**:

- **Explicado por artefacto o no robusto (4 celdas):** K28C 100µM y 1mM
  (dentro del nulo de artefacto), T324C 10µM (no replica con KV),
  E215C 100µM (signo inconsistente entre segmentadores).
- **Residuo que sobrevive ambas pruebas (3 celdas, una por constructo):**
  E215C 10µM, K28C 10µM, T324C 1mM — IM dwell→dwell 4-10× el techo del
  artefacto, replicada por dos segmentadores de familia distinta, mismo
  signo. Es un **candidato real pero no generalizado**; no alcanza el
  listón congelado de la etapa y NO asciende.

## 4. Lectura honesta

1. La historia simple "las pausas se recuerdan (anticorrelación)" está
   **muerta como evidencia**: el nulo calibrado fabrica esa firma sin
   memoria alguna. Este resultado es en sí mismo una contribución
   metodológica (cuánta memoria espuria fabrica una cadena
   MINFLUX→segmentación, medida con números).
2. Queda un residuo informacional real en 3 celdas grandes que NINGUNO de
   nuestros dos nulos de artefacto explica. Con el diseño actual no se
   puede subir más: los supuestos del nulo (renovación exacta, ruido
   gaussiano σ constante por traza) son ahora el límite.
3. **El siguiente juez no es otro análisis nuestro: es el gate 4.** Este
   informe, con sus dos congelaciones y sus CSV, es exactamente el paquete
   que debe ver un verificador externo (Ariga/Osaka o escuela Sivak).
   Ninguna carta se envía sin decisión de Raúl.

Ejecutado y firmado: Agente Khora, 2026-08-13.
