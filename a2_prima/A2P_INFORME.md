# A2′ — Informe: revalidación del estimador con largos reales (2026-08-13)

**Prereg v2** (congelado ANTES de ejecutar): SHA
`e4aebb2a6511e291e8a988ad8f237aece1643352cdb64229d47363506b756842`.
**Longitudes reales** (único conocimiento del dato usado): CSV congelado SHA
`ee4afd8ee758288b0ef23ed25a6bc7629493c6ad525dc49ce340ed30131957a8`.
**Datos:** exclusivamente sintéticos; ninguna trayectoria real leída.
**Ejecución:** 594 tareas (54 celdas × 11 condiciones × 500 réplicas,
1.000 permutaciones), semillas tipo 5/6 de raíz 20260812, 8 procesos
locales, 209 min. Estimadores de A2 reutilizados sin cambios
(`../a2/a2_estimadores.py`, SHA `dd784b2c…6eab9936`).

## VEREDICTO: PASA (criterios congelados de v2 §2)

**(i) Calibración bajo H0** (canal evento, PT): FPR media 0,0508 (k=4) y
0,0519 (k=8); máximos 0,076/0,078; **0 celdas de 54 significativamente
sobre 0,05 tras BH en ambas k** (lectura binomial adoptada). El conteo
literal (27/54 y 29/54 sobre 0,05) es el esperado por azar de muestreo con
500 réplicas (EE 0,010).

**(ii) Potencia en magnitudes plausibles** (PT k=4; umbral δ≤0,35 ó ρ≤0,20):
los 3 constructos ALCANZAN potencia ≥0,80 —
E215C y K28C vía limping δ=0,35 a 10 µM; T324C vía AR(1) ρ=0,20 y limping
δ=0,35 a 10 µM y 1 mM.

## Mapa de sensibilidad (mínima magnitud detectable, PT k=4, rango entre formas × mezclas)

| constructo | modo | 10 µM | 100 µM | 1 mM |
|---|---|---|---|---|
| E215C | limping δ | 0,35–0,50 | 0,50 | 0,35–0,50 |
| K28C  | limping δ | 0,35–0,50 | 0,50 | 0,35–0,50 |
| T324C | limping δ | 0,35 | 0,35–0,50 | 0,35 |
| E215C | AR(1) ρ | 0,20 | 0,35 | 0,20 |
| K28C  | AR(1) ρ | 0,20 | 0,35 | 0,20 |
| T324C | AR(1) ρ | 0,10–0,20 | 0,20 | 0,10–0,20 |

Lectura honesta del mapa, registrada ANTES de tocar datos reales:

1. **B1v2 puede ver correlaciones AR(1) moderadas (ρ ≥ 0,1–0,2) en las
   celdas grandes** (10 µM y 1 mM; 900–2.500 pares) y **solo limping fuerte
   (razón de dwells ≥ 2:1)**. Una memoria más débil que eso quedaría
   invisible: un REFUTA de B1v2 acota, no exonera.
2. **Las celdas de 100 µM son las menos sensibles** (400–800 pares). Un
   rechazo solo en 100 µM sería sospechoso (menos potencia ⇒ más
   probablemente artefacto que señal débil detectada justo ahí).
3. T324C es el constructo más informativo (mayor N: 305/111/287 trazas).

## Archivos y SHA-256

| archivo | SHA-256 |
|---|---|
| `a2p_barrido.py` | ver `git log` (comprometido pre-ejecución, commit 6267bde4) |
| `a2p_por_celda.csv` | (se registra en el commit de este informe) |
| `a2p_minima_detectable.csv` | (ídem) |
| `a2p_barrido.log` | (ídem) |

Los bloques crudos (`bloques_a2p/`, ~760 MB) quedan en disco fuera de git:
cada bloque es regenerable bit a bit desde la semilla (propiedad demostrada
para esta maquinaria en `../a2/RECEPCION_A2.md` §3).

## Qué habilita y qué no

- **Habilita B1v2** (estadístico agrupado, PT k=4, nulos de orden y paridad)
  **CONDICIONADA A LA FIRMA EXPRESA de Raúl** (prereg v2 §5.4). Sin firma,
  los datos reales no se tocan.
- No valida I_mem (B2 exigiría su propia adenda pública, v2 §4).
- No cubre memorias de alcance > lag 1 ni acopladas dwell↔paso (herencia de
  las limitaciones 2-3 del informe A2).

Ejecutado y firmado: Agente Khora, 2026-08-13.
