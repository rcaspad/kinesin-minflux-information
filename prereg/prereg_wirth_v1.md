# Preregistro — Frente Wirth: el motor como sonda del entorno (v1)

**Estado:** CONGELADO por orden de Raúl ("congela", 2026-08-14). SHA-256
de este fichero registrado en `memoria/decisiones.md` en el momento de la
congelación. Ningún estadístico de contenido de Wirth ha sido computado
antes de esta congelación (scouting de inventario aparte, con alcance
congelado previo: `resultados/scouting_wirth/SCOUTING_WIRTH_INFORME.md`).

**Diseño de origen:** adenda 2 de `DEBATE_RAIL_20260813.md` (convergencia
agente ↔ Nodo en tres rondas auditadas, commit `3fcd80e5`), con
predicciones direccionales y apuestas selladas en ese registro ANTES del
scouting.

---

## 1. Estimando y alcance de la reclamación (sonda, no mapa)

Pregunta única: **¿imprime el raíl nativo (microtúbulos axonales
preservados in situ, con código de tubulina y MAPs) una firma medible en
la estadística secuencial del paso de un motor calibrado, por encima de
la que el mismo motor muestra sobre raíl barajado (in vitro
reconstituido)?**

- El motor (constructos K28C y T324C, idénticos en ambos datasets) actúa
  como **sensor calibrado**: su respuesta sobre raíl barajado quedó
  caracterizada y atribuida en el frente Wolff (A1→E).
- El estimando es la **respuesta informacional del sensor** al raíl
  nativo frente al raíl barajado: reclamación AMPLIA sobre la estadística
  secuencial del paso.
- **Renuncia explícita**: ninguna reclamación espacial anclada. La
  palabra "sitio" queda prohibida en asociación con Wirth (sin constructo
  COM y sin marco de coordenadas común, E-1/E-2 no son ejecutables; la
  distinción raíl-puro vs raíl-vía-cabeza-marcada NO es separable aquí,
  y la reclamación lo declara).

## 2. Datos (congelados, ya verificados)

- **Wirth**: Zenodo 10718784 (MD5 ZIP `0fc8f961d74e3d1a93b3a944a686f515`
  idéntico al publicado; manifiesto SHA-256 en
  `raw/inbox/kinesin_minflux_wirth2024_zenodo10718784/`). 5 celdas:
  K28C × {50 µM, 500 µM, 5 mM}, T324C × {50 µM, 5 mM}.
- **Wolff** (regeneración): las 5 celdas EMPAREJADAS del dataset ya
  depositado (Zenodo 7565676): K28C × {10 µM, 100 µM, 1 mM},
  T324C × {10 µM, 1 mM}. E215C no tiene contraparte y queda fuera de la
  sustracción. T324C 100 µM sigue excluida (gate de C, sin cambio).
- **Emparejamiento** (rejillas exactamente ×5):
  bajo = 10↔50 µM; medio = 100↔500 µM (solo K28C); alto = 1↔5 mM.

## 3. Cadena congelada (sin cambios, por SHA)

- Segmentación: `resultados/c_crudo/c_kv_pipeline.py`, SHA
  `be3756b8ddadc91f274c0abaf75143ce5327a9f6fa76e72e65bbb1608c4d84d6`
  (parser 19-col, filtros de fotones 7/150 de los autores, proyección al
  eje principal, Kalafut-Visscher BIC, fusión <5 nm, eventos interiores
  con paso CON SIGNO). Aplica a ambos datasets sin modificación
  (verificado en scouting).
- Estimadores: `resultados/a2/a2_estimadores.py`, SHA
  `dd784b2c7ce7156dfa318a7317b816c0cda8c51d571fee111a1e357d6eab9936`.
  Primario **Panzeri-Treves k=4** (octiles con coarsening exacto);
  k=8 siempre reportado; NSB/plug-in/Miller-Madow como diagnósticos.
- Estadístico: **IM dwell→dwell agrupada por celda** sobre pares
  consecutivos dentro de traza (diseño B1v2), nulo de renovación
  (permutación de orden dentro de traza, 1.000 permutaciones, umbral
  p95 y p de permutación).
- **Simetría total**: ambos datasets pasan por ESTA congelación con
  semillas nuevas. **La cifra de Wolff usada en la sustracción es la
  regenerada aquí; C-2 y B1v2 quedan como referencia histórica, no como
  insumo.**

## 4. Semillas (raíz 20260812, `semilla_de`, tipos nuevos)

| tipo | uso |
|---|---|
| 17 | A2″ — datos sintéticos |
| 18 | A2″ — permutaciones |
| 19 | sustitutas del techo sintético — Wirth |
| 20 | sustitutas del techo sintético — Wolff (regeneración) |
| 21 | permutaciones sobre datos reales — Wirth |
| 22 | permutaciones sobre datos reales — Wolff (regeneración) |
| 23 | submuestreo de longitudes igualadas |

## 5. Gate A2″ (go/no-go DURO; primera puerta, antes de todo dato real)

Sintéticos de renovación (sin memoria secuencial por construcción) por
celda, **500 réplicas por celda** (10 celdas), construidos con las
**marginales reales** de cada dataset: multiset de longitudes de traza,
σ por traza (std(diff)/√2), rejillas temporales reales y clases de paso
nativas {4, 8, 16} nm con etiqueta de verdad conocida, pasados por la
cadena completa del §3. **Declaración de marginales:** las marginales no
contienen la pregunta (que es secuencial); precedente A2′ (multisets de
longitudes reales) y C-2 (σ y rejillas reales). Criterios congelados:

- **A2″-1 (calibración de FPR, lectura binomial como A2):** a α = 0,05,
  la fracción de réplicas con falso positivo debe caer dentro de la
  banda binomial del 95 % en cada celda.
- **A2″-2 (resolubilidad de clases a σ de Wirth):** exactitud de
  asignación de clase (fronteras 6/12 nm sobre |paso| estimado por la
  cadena vs clase verdadera) **global ≥ 80 % y por clase ≥ 2/3** en las
  celdas de Wirth.
- **Muerte honrosa (a):** si A2″-1 o A2″-2 fallan, el frente muere sin
  tocar contenido real y se registra públicamente.

## 6. Diseño primario — doble sustracción por bandas

Por cada celda D de cada dataset:

1. **Techo sintético**: **200 sustitutas de renovación** por celda
   (marginales reales, cadena completa) → distribución IM_sur;
   techo = p95(IM_sur).
2. **exceso_D** = IM_real_D − mediana(IM_sur_D).
3. **Criterio 1 (por celda de Wirth):** IM_real > techo p95 de su propia
   celda.
4. **Criterio 2 (efecto raíl, por par emparejado):**
   efecto_raíl = exceso_Wirth − exceso_Wolff debe superar el p95 de la
   **banda combinada**: distribución de diferencias de sustitutas
   centradas, (IM_sur_Wirth − mediana_Wirth) − (IM_sur_Wolff −
   mediana_Wolff), con 10.000 emparejamientos aleatorios (semillas §4).
5. **Confirmación requiere criterio 1 Y criterio 2** en el test primario.

## 7. Jerarquía de evidencia por [ATP] (asimetría de K_M, congelada)

v/V_max = [ATP]/(K_M+[ATP]); K_M nominal 50 µM, robustez verificada en
28-100 µM (razones cinéticas del par: alto ×1,02-1,08; medio ×1,21-1,67;
bajo ×2,4-3,7).

- **Primario: K28C alto (1 mM ↔ 5 mM)** — par cinéticamente igualado.
  Sin corrección múltiple (test único).
- **Apoyo: K28C medio y T324C alto** — Holm dentro de la familia (m=2).
- **Exploratorio: K28C bajo y T324C bajo** — reportados y etiquetados
  como hipótesis; un efecto presente SOLO en pares bajos no asciende a
  resultado (ambigüedad de régimen cinético declarada).
- Si el paper de Wirth reporta K_M efectiva propia, se cita como
  contraste; no altera esta jerarquía.

## 8. Robustez congelada

- **Longitudes igualadas (confundidor longitud↔[ATP]):** 200
  submuestreos (semilla tipo 23) igualando la distribución de longitudes
  de traza entre condiciones de cada constructo; la decisión del
  criterio 1 del test primario debe mantenerse en **≥ 80 %** de los
  submuestreos. Si no, la comparación se reagrupa por longitud y se
  reporta como limitación.
- **Comparabilidad (gate formal en ejecución):** eventos KV / pasos de
  los autores dentro de 50-200 % por celda (pinneo previo: Wolff
  1,23-1,61, Wirth 1,35-1,63 — propiedad del par de segmentadores).
- **Guard de paridad** (limping): nulo de paridad reportado junto al de
  renovación en el test primario, como en B1v2.

## 9. Muertes honrosas (dos caminos legítimos de cierre)

- **(a)** A2″ falla (§5) → cierre sin tocar datos reales.
- **(b)** El test primario no supera el criterio 1 → **resultado nulo
  limpio**: "el raíl nativo no deja firma medible en la estadística del
  paso a esta resolución y N" — publicable, cierra el frente.
- Sin etapa de rescate: no habrá C-v3 ni umbrales retocados. Lo que no
  salga de aquí, no sale.

## 10. Terminología (entra en `TERMINOLOGIA.md` con esta congelación)

**sonda del entorno** (el motor calibrado como sensor del raíl);
**raíl nativo** (microtúbulos axonales preservados in situ, Wirth);
**raíl barajado** (microtúbulos repolimerizados in vitro, Wolff);
**doble sustracción** (diseño del §6); **techo sintético** (p95 de las
sustitutas de renovación por la cadena completa); **A2″** (gate de
resolubilidad a ruido de Wirth). Prohibido: "sitio" asociado a Wirth;
"in vivo" para Wirth (es in situ).

## 11. Reproducibilidad

Todo determinista desde la raíz 20260812 vía `semilla_de` (§4); scripts
y tablas se depositarán en `resultados/a2_bis/` (A2″) y
`resultados/wirth_sustraccion/` (diseño primario); los informes citarán
los SHA de este prereg y de los scripts ejecutados. Idioma: español
verbatim (regla del repo público).

Redactado por el Agente Khora por orden de congelación de Raúl.
2026-08-14.
