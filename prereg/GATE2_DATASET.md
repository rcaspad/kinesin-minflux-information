# Gate 2 — Caza del dataset (2026-08-12)

**Estado:** ✅ **COMPLETADO** (2026-08-12). Raúl autorizó la descarga
("autorizo la descarga"); el candidato A está depositado en
`../raw/inbox/kinesin_minflux_wolff2023_zenodo7565676/` con MD5 verificado
IDÉNTICO al publicado (`473e28f9…350e`), SHA-256 propio
`69a335b54e732b8496850874ab6fa02877dc81a5f4a955f12399238dfb9a1930`, e
inventario de 194 archivos / 1.03 GB sin extraer (manifiesto completo:
`MANIFEST.md` del depósito). El candidato D (Ariga) sigue abierto como vía
para F1-B con señal de forzado + gate 4.

## Candidato A — PRINCIPAL (verificado completo)

| Campo | Valor |
|---|---|
| Depósito | "MINFLUX dissects the unimpeded walking of kinesin-1" |
| Repositorio | Zenodo, DOI `10.5281/zenodo.7565676` (2023-01-24) |
| Autores del depósito | Jan Otto Wolff, Lukas Scheiderer (MPI Medical Research; ORCIDs verificados) |
| **Licencia** | **CC BY 4.0** — uso permitido con atribución |
| Archivo | `Repository for MINFLUX dissects the unimpeded walking of kinesin-1.zip` |
| Tamaño | **38.213.434 bytes (36,4 MB)** — un único ZIP |
| **MD5 publicado** | `473e28f9444427c55540362f84df350e` |
| Contenido declarado | Datos MINFLUX crudos + scripts Matlab; stepping de kinesina-1 sobre microtúbulos a **1,7 nm / 1 ms**; rotaciones de tallo y cabezas; dinámica de unión de ATP |
| Paper compañero | Wolff, Scheiderer, Engelhardt, Engelhardt, Matthias & **Hell** (2023). *Science* 379, 1004-1010. DOI `10.1126/science.ade2650` — **verificado CrossRef, revisado por pares, NO retirado** |
| Depósitos hermanos | Scripts: `10.5281/zenodo.7442902`; código+secuencias: `10.5281/zenodo.7837326` (contiene `kinesin sequences ade2650.zip`, 143.972 bytes) |

**Advertencia de idoneidad (honesta):** las trayectorias son de una cabeza
marcada (un fluoróforo) — resolución extraordinaria pero SIN señal de forzado
externo registrada. Para F1-B (descomposición de Still, que necesita señal
estocástica + respuesta) el mapeo no es directo; el prereg deberá definir qué
juega el papel de x_t (p. ej. dinámica de unión de ATP detectada, o reformular
sobre la estructura interna del paso). Esta limitación se declara ANTES de
descargar, no después.

## Candidato B — alternativa publicada del mismo linaje

"Uncovering kinesin dynamics in neurites with MINFLUX" (Wirth et al., *Comm.
Biology* 7, 661, 2024; DOI datos `10.5281/zenodo.10718784`, CC-BY-4.0,
~7,1 GB): kinesina **in vivo** en neuritas. Más masa, menos control; segunda
opción.

## Candidato C — CUARENTENA (presa del gate)

"Dual-color MINFLUX: Kinesin-1 takes Chassé-Inchworm steps" (Edmond,
`10.17617/3.2QKPKH`, CC0, 101,4 MB, dos cabezas simultáneas — sobre el papel
el ideal para flujos entre cabezas). **DESCARTADO**: su preprint compañero
(bioRxiv `10.1101/2024.03.05.583551`) fue **RETIRADO por los autores** —
"experimentos mejorados no respaldan el mecanismo propuesto; no deseamos que
este trabajo sea citado". Un dataset cuyo companion cayó puede llevar el
artefacto dentro. El gate existió para esto. (Nota: la revisión #20 citaba el
chassé-inchworm como hallazgo 2024 — las revisiones también arrastran
retiradas.)

## Candidato D — el mejor mapeo conceptual a Still, pendiente de averiguar

Ariga et al. 2021, "Noise-induced acceleration of single molecule kinesin-1"
(arXiv `2012.04214`): aplican **fuerza estocástica externa** al motor y miden
la respuesta — literalmente el montaje del teorema de Still (señal x_t +
sistema s_t, ambos registrados por construcción). Disponibilidad de los datos:
desconocida (era PRL, sin depósito localizado). **Contactar a Ariga serviría
doble: acceso a datos + candidato natural a verificador externo humano
(gate 4).** [PENDIENTE]

## Siguiente paso (bloqueado en autorización)

Con el OK de Raúl: descargar el ZIP del candidato A (36,4 MB, Zenodo) a
`../raw/inbox/`, verificar MD5 contra el publicado (`473e28f9…350e`), calcular
SHA-256 propio, y registrar el depósito. Sin el OK, nada se descarga.
