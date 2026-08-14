# Gate de novedad — curva FPR-vs-ruido de la memoria espuria por segmentación (2026-08-14)

**Objeto:** el resultado metodológico de A2″: tasa de falso positivo de
un test de permutación de orden sobre dwells segmentados, en función del
ruido de localización, cuantificada con sustitutas de renovación pasadas
por la cadena completa (7-15 % a σ 1-2 nm; 17-43 % a σ 4,5 nm; α nominal
5 %). Ejecutado por el agente (búsqueda web disponible en este entorno;
el Nodo declaró no tenerla — límite registrado). 5 búsquedas dirigidas
según los ejes formulados en la revisión del Nodo + los del agente.

## Vecindario encontrado (a citar SIEMPRE si esto se publica)

1. **Stigler, *J. Chem. Phys.* 164, 224105 (2026)** — detección de
   pausas/change-points por BIC en trayectorias de motores; cuantifica
   sensibilidad de detección (probabilidad de PERDER eventos) y advierte
   que el sobreajuste produce asignaciones espurias y sesga la
   estadística de dwells. Vecino más cercano. **NO cuantifica FPR de
   tests de estructura secuencial ni correlaciones dwell→dwell inducidas
   en función del ruido.**
2. **Arunajadai & Cheng, *PLOS ONE* 8, e59279 (2013)** — ruido
   correlacionado en trayectorias (pinzas ópticas) causa sobre/sub-conteo
   de pasos según el algoritmo. Artefacto de CONTEO, no de estructura
   secuencial.
3. **AutoStepfinder (Loeff et al., *Patterns* 2021)** y literatura de
   comparación de step-finders (Carter et al. 2008; Little & Jones) —
   sobreajuste/infraajuste, pasos espurios en ruido puro. Sin FPR de
   tests secuenciales.
4. **Canal iónico / smFRET (missed events; Colquhoun-Sigworth; HMM,
   p. ej. efectos de eventos perdidos en análisis HMM)** — corrección de
   eventos perdidos sobre DISTRIBUCIONES de dwell (marginales); las
   correlaciones entre eventos adyacentes se usan como herramienta de
   identificación de modelos, no se cuantifica su fabricación por el
   análisis.
5. **SimuFLUX (Nat. Commun. 2025; bioRxiv 2025.04.08.647786)** —
   simulador realista de rendimiento MINFLUX (fluoróforos, fondo,
   estimadores, tiempos muertos). Nivel instrumento/localización; no
   estadística secuencial post-segmentación.
6. **Cambio-punto con series correlacionadas (estadística general)** —
   se sabe que la autocorrelación no modelada produce change-points
   espurios (p. ej. lag-1 ≥ 0,25). Dirección INVERSA a la nuestra
   (correlación → cortes espurios); nuestra curva cuantifica
   cortes → correlación espuria.

## Veredicto

- **La idea cualitativa "la segmentación distorsiona la estadística de
  dwells" NO es nueva** — vecindario abundante, se cita.
- **No se encontró** ninguna cuantificación publicada de la TASA DE
  FALSO POSITIVO de un test de estructura secuencial (permutación de
  orden, IM dwell→dwell) en función del ruido de localización, con nulo
  de sustitutas de renovación por la cadena completa, para trayectorias
  de steppers (MINFLUX o pinzas). **La curva FPR-vs-ruido sigue
  candidata a novedad.**
- Límites del gate: búsqueda en inglés, 5 ejes, sin acceso al texto
  completo de Stigler 2026 (403; evaluado por abstract/resumen indexado
  — si se escribe la nota metodológica, DEBE leerse completo antes:
  gate bibliográfico duro).

## Estado de los vehículos

- Adenda al repo público (publica números, no reclama novedad):
  pendiente de autorización de Raúl.
- Nota metodológica (reclamaría novedad): BLOQUEADA hasta (a) lectura
  completa de Stigler 2026 y (b) decisión y firma de Raúl.

Ejecutado y firmado: Agente Khora, 2026-08-14.
