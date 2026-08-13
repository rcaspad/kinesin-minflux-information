# Errata v1.1 al prereg de kinesina (pública, 2026-08-12 noche)

**El prereg v1.0 (SHA congelado `c73041ab…cc6c71b2`) permanece inmutable; esta
errata lo acompaña sin modificarlo.** Corrige METADATOS erróneos; no toca
hipótesis, estimadores, nulos, semillas ni umbrales.

## Corrección

El §1 del prereg describía los constructos como "tallo (E215C) y cabeza
(K28C)". **Era una suposición mía (agente Khora), y es falsa.** El texto
completo del paper (depositado:
`../raw/inbox/wolff2023_fulltext_biorxiv_501426.html`, SHA-256
`56cacd64…a66e370a`) establece:

| Constructo | Posición real del marcaje |
|---|---|
| E215C | **Cabeza** — extremo C-terminal de la lámina β6 del dominio motor |
| K28C | **Cabeza** — cisteína en posición 28 del dominio motor |
| T324C | **Cabeza** — extremo C-terminal de la hélice α6, junto al neck linker |
| N356C | **Coiled-coil (tallo)** — usado para center-of-mass; en el dataset solo aparece como control estacionario |

Consecuencia interpretativa (no de diseño): los tres constructos con
trayectorias de caminar son de cabeza → pasos regulares de 16 nm con subpasos
de ~8 nm, más detectables a ATP bajo. Las medianas de A1 (8,5-16 nm según
celda) son mezclas paso/subpaso, como los propios histogramas del paper.

## Registro de honestidad

La discrepancia fue detectada por el propio gate A1 y declarada en
`analisis/A1_INFORME_INTERINO.md` **antes** de conocer la respuesta del paper.
El diseño del prereg (B1/B2 sobre secuencias de eventos por constructo ×
condición) no depende de qué parte está marcada, y queda intacto.
