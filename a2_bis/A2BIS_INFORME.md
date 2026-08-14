# A2″ — Informe del gate (2026-08-14): NO-GO. Muerte honrosa (a) del frente Wirth

**Congelación:** criterios en `prereg_wirth_v1.md` §5 (SHA
`1d6ab5d5…2586a2de`); implementación en `A2BIS_CONGELACION.md` (SHA
`e1f39023…63de337b`), comprometida ANTES de resultados (commit
`083e0d06`). Script final `a2bis_analisis.py` SHA
`41609aa8c5a8ebf00eab425e1fbfcda7445bdd51a7c0a662ffe7b943686ecf32`.
**Transparencia de cambio:** tras el commit de congelación, el script
recibió UNA corrección (el cargador asumía un `.txt` por celda; las
celdas de Wolff tienen varios por fechas/muestras — misma resolución que
`segmentar_celda` de C). El cambio afecta solo a la carga de rejillas;
ningún criterio, semilla ni umbral fue tocado. Semillas tipos 17/18,
raíz 20260812. 500 réplicas × 10 celdas; **ningún estadístico de
contenido real computado**.

## 1. Resultados

Banda binomial exacta 95 % para FPR a α = 0,05 con n = 500: **[16, 35]**.

| celda | FP/500 | FPR | A2″-1 | exact. global | c4 | c8 | c16 | A2″-2 |
|---|---|---|---|---|---|---|---|---|
| Wolff K28C 10µM | 34 | 6,8 % | OK | 0,872 | 0,502 | 0,890 | 0,868 | (ref) NO |
| Wolff K28C 100µM | 38 | 7,6 % | FUERA | 0,831 | 0,405 | 0,847 | 0,839 | (ref) NO |
| Wolff K28C 1mM | 31 | 6,2 % | OK | 0,858 | 0,487 | 0,860 | 0,871 | (ref) NO |
| Wolff T324C 10µM | 49 | 9,8 % | FUERA | 0,806 | 0,494 | 0,839 | 0,846 | (ref) NO |
| Wolff T324C 1mM | 75 | 15,0 % | FUERA | 0,810 | 0,441 | 0,825 | 0,836 | (ref) NO |
| **Wirth K28C 50µM** | **102** | **20,4 %** | FUERA | 0,785 | 0,361 | 0,807 | 0,804 | **NO** |
| **Wirth K28C 500µM** | **217** | **43,4 %** | FUERA | 0,772 | 0,307 | 0,752 | 0,785 | **NO** |
| **Wirth K28C 5mM** | **97** | **19,4 %** | FUERA | 0,777 | 0,339 | 0,773 | 0,805 | **NO** |
| **Wirth T324C 50µM** | **87** | **17,4 %** | FUERA | 0,743 | 0,334 | 0,755 | 0,789 | **NO** |
| **Wirth T324C 5mM** | **99** | **19,8 %** | FUERA | 0,768 | 0,359 | 0,766 | 0,798 | **NO** |

- **A2″-1 (calibración de FPR): FALLA** — 8 de 10 celdas fuera de banda;
  en Wirth, el test de permutación dispara al 17-43 % con α nominal 5 %.
- **A2″-2 (resolubilidad a σ de Wirth): FALLA en 5/5 celdas de Wirth** —
  exactitud global 0,74-0,79 (< 0,80) y clase de 4 nm en 0,31-0,36
  (≪ 2/3). En Wolff (referencia, no gate) la clase de 4 nm tampoco llega
  (0,41-0,50).
- Cobertura (diagnóstico): 0,96-1,16 detecciones/verdad — el fallo no es
  de conteo, es de identidad.

**VEREDICTO: NO-GO. Se aplica la muerte honrosa (a) del prereg §9: el
frente Wirth muere sin haber abierto un solo dato real, y sin etapa de
rescate ("lo que no salga de aquí, no sale").**

## 2. Diagnóstico (por qué muere, y qué significa)

1. **El fallo de FPR es el artefacto de C-2 escalando con el ruido.** Los
   sintéticos son renovación pura (sin memoria por construcción); toda la
   estructura secuencial que el test detecta la FABRICA la cadena de
   segmentación. A σ de Wolff (~1-2 nm) esa fabricación infla el FPR a
   7-15 %; a σ de Wirth (~4,5 nm) lo infla a 17-43 %. **El estadístico
   primario congelado (IM agrupada vs permutación de orden) no está
   calibrado a nivel de cadena, y a ruido de Wirth está roto sin
   paliativos.** Es la generalización cuantitativa exacta del hallazgo de
   C-2 — y su confirmación más dura hasta la fecha.
2. **La clase de 4 nm no existe para esta cadena.** La fusión < 5 nm
   (congelada desde C-v2) se come los pasos pequeños en ambos datasets;
   a σ de Wirth, además, las clases 8/16 pierden filo (0,75-0,81). El
   alfabeto de 3 clases {4,8,16} nm no es transportable a 4,5 nm de
   ruido.
3. **Nota sobre el barrido:** el peor FPR (43 %) es Wirth K28C 500 µM,
   la celda con dwells más cortos respecto a su rejilla — coherente con
   que el artefacto crece cuando la segmentación fragmenta más.

## 3. Lo que esta muerte NO toca (y lo que retro-ilumina)

- **Los veredictos del frente Wolff (B1v2, C, D, E) quedan como están.**
  Sus reclamaciones finales nunca descansaron en el p de permutación a
  secas: el residuo de C superaba 4-10× el techo del artefacto calibrado
  (C-2), y E atribuyó todo al ciclo de paso por controles inmunes al
  artefacto (COM, prueba reina). A2″ retro-confirma que ese diseño era
  necesario: las lecturas "solo permutación" (C-1) eran anticonservadoras
  a nivel de cadena, exactamente como C ya declaró al fallar su nulo.
- **La aportación №2 de la carta del gate 4 sale REFORZADA:** la memoria
  espuria inducida por segmentación tiene ahora una cuantificación de
  FPR en función del ruido (7-15 % a 1-2 nm; 17-43 % a 4,5 nm; α nominal
  5 %). Es el primer resultado de este tipo que conocemos y es
  exportable a cualquier análisis secuencial de datos MINFLUX/steppers
  segmentados.
- El dataset Wirth queda descargado, verificado y VIRGEN (ningún
  estadístico de contenido). Si alguna vez existiera un replanteo, sería
  un preregistro NUEVO que debería confrontar explícitamente la cláusula
  de muerte de `prereg_wirth_v1.md` §9 — decisión que no pertenece al
  agente.

## 4. Registro

- La muerte se registra en `memoria/decisiones.md` y en `PLAN.md`
  (frente Wirth: CERRADO por A2″).
- Las predicciones direccionales del debate (H_raíl vs H_motor para el
  raíl nativo) quedan SIN RESOLVER — no falsadas, no confirmadas: el
  instrumento estadístico disponible no está calibrado al ruido del
  único dataset in situ existente. Que quede escrito con esa precisión.

Ejecutado y firmado: Agente Khora, 2026-08-14.
