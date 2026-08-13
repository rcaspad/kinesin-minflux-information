# B1 — Informe de ejecución (2026-08-13)

**Prereg:** v1.0 FIRMADO, SHA `c73041ab…cc6c71b2` (+ errata v1.1).
**Congelación previa:** `B1_CONGELACION.md` (lectura binomial adoptada por
Raúl → A2 PASA; estimador congelado: Panzeri-Treves k=4).
**Script:** `b1_analisis.py` (determinista, semillas tipo 3/4 derivadas de
20260812; reutiliza `../a2/a2_estimadores.py` sin modificar).

## VEREDICTO: INCONCLUSO (criterio de muerte del prereg §3-B1, forma extrema)

**Cero trazas elegibles (≥ 50 eventos) en las 9 celdas** (constructo ×
[ATP]). El prereg declara INCONCLUSO si N elegible < 10 en todas las
condiciones; aquí N elegible = 0 en todas.

### La causa, con los números

Longitudes de traza en las tablas de pasos de los autores (eventos válidos
`end_flag==0, tau>0` por traza):

| celda | trazas | mediana | p90 | máx | ≥20 | ≥30 | ≥50 |
|---|---|---|---|---|---|---|---|
| E215C 10µM | 207 | 6 | 14 | 28 | 3 | 0 | 0 |
| E215C 100µM | 60 | 7 | 14 | 19 | 0 | 0 | 0 |
| E215C 1mM | 96 | 7 | 14 | 18 | 0 | 0 | 0 |
| K28C 10µM | 113 | 11 | 24 | 35 | 23 | 5 | 0 |
| K28C 100µM | 56 | 9 | 14 | 16 | 0 | 0 | 0 |
| K28C 1mM | 71 | 14 | 19 | 24 | 6 | 0 | 0 |
| T324C 10µM | 318 | 8 | 17 | 26 | 16 | 0 | 0 |
| T324C 100µM | 119 | 6 | 14 | 26 | 6 | 0 | 0 |
| T324C 1mM | 299 | 9 | 15 | 22 | 2 | 0 | 0 |
| **TOTAL** | **1.339** | **8** | **16** | **35** | 56 | 5 | 0 |

La traza más larga del dataset entero tiene 35 eventos. El umbral de 50 no
falla por margen: es **estructuralmente inalcanzable** — la fotoblanqueo
limita las caminatas MINFLUX de un fluoróforo a ~10 pasos observados.

### Por qué el umbral era 50 y por qué eso no se corrige ahora

El prereg §4 congeló "≥ 50 eventos" ANTES de abrir ningún archivo de
trayectorias (declaración de integridad §0), calibrándolo sobre la rejilla
sintética plausible (largos 50-200). El precio de firmar a ciegas es este:
el umbral resultó incompatible con la física del instrumento. Bajarlo ahora,
con las longitudes ya vistas, sería el jardín de senderos que el método
prohíbe. **El umbral no se toca; el veredicto es INCONCLUSO.**

### Lo que queda intacto (y es valioso)

1. **Ningún estadístico de información fue computado sobre datos reales.**
   El pipeline murió en la puerta de elegibilidad; de las trayectorias
   reales solo se han leído longitudes y (en A1) medianas agregadas de
   paso/dwell. El contenido informacional del dataset sigue **virgen** para
   cualquier prereg v2.
2. La cadena A2 completa (estimadores calibrados, control interno,
   mínima detectable) es reutilizable tal cual para una revalidación con
   rejilla realista.
3. INCONCLUSO ≠ REFUTA: B1 no dice nada sobre si el caminar de la kinesina
   tiene memoria. Dice que ESTE diseño sobre ESTE dataset no puede medirla.

### Hecho adicional registrado (sin abrir contenido)

La rejilla de A2 (largos 50-200, y por tanto sus potencias) queda FUERA del
régimen real (largos 5-35): incluso un rediseño v2 con estadístico agrupado
por celda (la variante que A2 sí validó con n_tray 50-300) exigiría una
**A2' recalibrada a las longitudes reales** antes de tocar los datos. Las
longitudes necesarias para esa rejilla ya son conocimiento legítimo (esta
tabla).

## Opciones para el replanteo (decisión de Raúl; prereg §3: "el frente se
replantea sin tocar los datos")

- **(1) Prereg v2 pública** con estadístico agrupado por celda (validado en
  A2 para n_tray≥50) + A2' con rejilla de largos reales {5, 8, 16, 35} para
  medir FPR y potencia alcanzable ANTES de decidir si ejecutar. Si A2' da
  potencia inservible a estos largos, se declara NO MEDIBLE y se cierra.
- **(2) Otro dataset con caminatas largas:** candidato B (Wirth 2024,
  kinesina in vivo en neuritas, 7,1 GB, Zenodo 10718784) — procesividades
  in vivo ~1 µm ≈ 100+ pasos; o candidato D (Ariga, forzado estocástico —
  sinergia con gate 4). Ambos exigirían su propio gate 2 + prereg.
- **(3) Cerrar el frente F1-B** sobre este dataset y evaluar el tripwire
  coste/valor del frente completo.

Nada se ejecuta hasta esa decisión. Gate 4 (verificador externo humano)
sigue pendiente en cualquier rama.

Ejecutado y firmado: Agente Khora, 2026-08-13.
