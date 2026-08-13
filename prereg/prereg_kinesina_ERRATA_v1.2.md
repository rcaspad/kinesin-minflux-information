# Errata v1.2 al prereg de kinesina (pública, 2026-08-13 noche)

**El prereg v1.0 (SHA `c73041ab…cc6c71b2`) y la errata v1.1 permanecen
inmutables; esta errata corrige un error factual DE LA PROPIA ERRATA
v1.1.** No toca hipótesis, estimadores, nulos, semillas ni umbrales de
ningún preregistro.

## Corrección

La errata v1.1 afirmaba que N356C (marcaje en coiled-coil/tallo,
center-of-mass) "en el dataset solo aparece como control estacionario".
**Es falso.** El dataset contiene, además del control estacionario
(`stationary_N356/`, 3 ficheros), trazas de N356C **caminando** con sus
tablas de pasos: `N356C/DOL1/10uM`, `N356C/DOL1/1mM` y `N356C/DOL2/1mM`
(ficheros crudos + `allsteps_reeval.xls`).

Origen del error: la v1.1 se redactó desde el texto del paper (que usa
N356C para COM y como control estacionario en las figuras revisadas) sin
inventariar el directorio N356C del depósito. Detectado el 2026-08-13 al
corregir la ruta del test D-4.

## Consecuencia (positiva, sin efecto retroactivo)

Ningún análisis ejecutado usó trazas N356C caminando, así que ningún
resultado cambia. Hacia delante: las trazas COM de N356C habilitan el
desempate cabeza-vs-tallo (fragmentación/vaivén de cabeza vs identidad de
sitio) propuesto en `resultados/d_atribucion/D_INFORME.md` §2-§3, como
etapa E si el gobernador la ordena.

Registro de honestidad: tercera errata pública del proyecto; las tres
detectadas por nuestros propios gates.
