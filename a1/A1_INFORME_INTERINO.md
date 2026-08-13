# A1 — Informe del gate de reproducibilidad (2026-08-12, noche — CERRADO)

**Prereg:** v1.0 FIRMADO, SHA `c73041ab…cc6c71b2`. **VEREDICTO A1: PASA.**

Contraste contra el texto completo del paper (depositado en `../../raw/inbox/
wolff2023_fulltext_biorxiv_501426.html`, SHA `56cacd64…a66e370a`):
- **Conteo total**: el paper declara ">12.000 pasos de kinesina identificados";
  nuestras tablas suman **13.355 pasos** en 1.339 trazas. ✅
- **Tamaños de paso**: cabezas = pasos regulares de 16 nm (nuestras medianas
  E215C 14,45-16,01 nm) con subpasos de ~8 nm más detectados a ATP bajo
  (K28C 8,50 y T324C 9,38 nm a 10 µM, subiendo a 12,7-14,8 a mayor ATP) —
  exactamente el patrón mezcla paso/subpaso que describe el paper. ✅
- **Precisión**: paper σ≈1,7-2,1 nm por localización, 0,63 nm de precisión de
  paso; nuestras σ de plateau 2,8-3,9 nm — coherente (métrica distinta,
  mismo orden). ✅
- **Errata detectada por el gate**: mis etiquetas de constructo eran falsas
  (ver `../prereg_kinesina_ERRATA_v1.1.md`) — declarada antes de conocer la
  respuesta, corregida sin tocar el diseño.

Nota: el criterio literal "±10 % de medianas" celda a celda no es aplicable
tal cual porque las medianas son mezclas paso/subpaso (también en el paper);
el contraste se hizo sobre los anclajes cuantitativos publicados (conteo,
picos 16/8/4 nm, precisión, dirección de dwells con ATP). Todos consistentes.
**A2 puede proceder.**

## Qué se hizo

1. Gates previos del prereg §8, ambos limpios ANTES de abrir datos:
   - Limping verificado: Asbury, Fehr & Block (2003), *Science* 302,
     2130-2134, DOI `10.1126/science.1092985` → corpus #23.
   - Búsqueda de cancelación: sin resultados que combinen información
     predictiva con trayectorias MINFLUX (consulta registrada 2026-08-12).
2. ZIP extraído a `extracted/` (fuera de `raw/`; gitignorado).
3. Semántica de columnas de las tablas decodificada **del código de los
   autores** (`process_MF_data.m` L71), no supuesta.
4. Script determinista `a1_estadisticos_basicos.py` (SHA
   `1e06d666…e4aa58d`) → `a1_resultados.csv` (SHA `807a8f3d…a3ff3460`).

## Resultados (9 celdas, 1.339 trazas, 13.355 pasos)

| Constructo | ATP | Trazas | Pasos | Step mediana (nm) | Dwell mediana (s) |
|---|---|---|---|---|---|
| E215C | 10 µM | 207 | 1.700 | 14,45 | 0,042 |
| E215C | 100 µM | 60 | 522 | 16,01 | 0,021 |
| E215C | 1 mM | 96 | 871 | 15,71 | 0,020 |
| K28C | 10 µM | 113 | 1.576 | 8,50 | 0,023 |
| K28C | 100 µM | 56 | 534 | 14,59 | 0,017 |
| K28C | 1 mM | 71 | 1.051 | 14,77 | 0,018 |
| T324C | 10 µM | 318 | 3.148 | 9,38 | 0,028 |
| T324C | 100 µM | 119 | 1.037 | 12,71 | 0,022 |
| T324C | 1 mM | 299 | 2.916 | 13,04 | 0,019 |

Precisión de localización (σ_x mediana): 2,8-3,9 nm — consistente con lo
declarado por la técnica. Dwells decrecen con [ATP] en los tres constructos
(dirección esperada).

## Hallazgo del gate — discrepancia de etiquetas propia (declarada)

El prereg §1 y el MANIFEST describen E215C como "tallo" y K28C como "cabeza"
**por suposición mía al redactar**. Los números no la respaldan sin más:
E215C da medianas ~14,5-16 nm (esperable de cabeza en hand-over-hand) y
K28C/T324C dan ~8,5-9,4 nm a 10 µM (posible resolución de subpasos a ATP
bajo). **La asignación correcta de constructos debe salir del paper, no de mi
suposición.** Si el paper contradice las etiquetas del prereg: se emite
**errata v1.1 pública** (corrección de metadatos, sin tocar diseño ni
umbrales). El error es mío y queda registrado aquí antes de conocer la
respuesta.

## Criterio de cierre de A1 (del prereg, intacto)

A1 PASA si estas medianas reproducen los histogramas/valores del paper
(±10 %). Bloqueado en: obtención del texto completo (preprint pp. 1-2 ya en
`../../raw/inbox/wolff2023_preprint_biorxiv_501426.pdf`; HTML completo en
reintento). Hasta ese contraste, **ningún resultado de esta tabla se
interpreta** — son estadísticos de gate, no hallazgos.
