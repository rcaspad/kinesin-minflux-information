# -*- coding: utf-8 -*-
"""
B1v2 — Memoria más allá de la renovación, estadístico agrupado por celda
(prereg v2 FIRMADO, SHA post-firma aa461a6e28b4aaa8198b3c44fd480ce0d7d14f6169e95a887460ef6836f22e6f;
borrador pre-firma e4aebb2a…b756842; A2' = PASA).

PRIMERA lectura de la estructura temporal de las trayectorias reales.
Determinista: nulo de orden semilla_de(7, id_celda_real, 0); nulo de
paridad semilla_de(8, id_celda_real, 0). Estimadores de A2 sin cambios.

Decisión primaria: canal evento, Panzeri-Treves, k=4, nulo orden,
clases de paso por centros {4,8,16}. Sensibilidades obligatorias: k=8 y
terciles. Guard de limping: nulo de paridad.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(AQUI, "..", "a2"))
sys.path.insert(0, os.path.join(AQUI, "..", "b1"))
import a2_estimadores as A          # noqa: E402
from b1_analisis import (cargar_trazas, CONSTRUCTOS, CONDICIONES,  # noqa: E402
                         BORDES_PASO_NM)

N_PERM = A.N_PERMUTACIONES
MIN_EVENTOS = 2                     # prereg v2 (>= 1 par)
CANALES = ["evento", "dwell", "paso", "dwell_paso"]
ALFA = 0.05


def conteos_pares_variable(trazas_sim, alfabeto, rng, paridad=False):
    """(1+N_PERM, A²) agrupado; permutación dentro de cada traza.

    paridad=True: permuta SOLO dentro de las posiciones pares e impares de
    cada traza (nulo de paridad, guard de limping — prereg v2 §3).
    """
    a2 = alfabeto * alfabeto
    salida = np.zeros((1 + N_PERM, a2), dtype=np.int64)
    for sim in trazas_sim:
        cod = sim[:-1].astype(np.int32) * alfabeto + sim[1:]
        salida[0] += np.bincount(cod, minlength=a2)
        M = np.tile(sim, (N_PERM, 1))
        if paridad:
            M[:, 0::2] = rng.permuted(M[:, 0::2], axis=1)
            M[:, 1::2] = rng.permuted(M[:, 1::2], axis=1)
        else:
            M = rng.permuted(M, axis=1)
        codp = M[:, :-1].astype(np.int32) * alfabeto + M[:, 1:]
        desplaz = (np.arange(N_PERM, dtype=np.int32) * a2)[:, None]
        salida[1:] += np.bincount((codp + desplaz).ravel(),
                                  minlength=N_PERM * a2).reshape(N_PERM, a2)
    return salida


def holm(pvals: np.ndarray) -> np.ndarray:
    orden = np.argsort(pvals)
    m = len(pvals)
    adj = np.empty(m)
    corr_max = 0.0
    for rango, i in enumerate(orden):
        corr = min(1.0, (m - rango) * pvals[i])
        corr_max = max(corr_max, corr)
        adj[i] = corr_max
    return adj


def main():
    filas = []
    for ci, constructo in enumerate(CONSTRUCTOS):
        for di, cond in enumerate(CONDICIONES):
            id_real = ci * 3 + di
            trazas = [t for t in cargar_trazas(constructo, cond)
                      if len(t) >= MIN_EVENTOS]
            dw_pool = np.concatenate([[e[0] for e in t] for t in trazas])
            sx_pool = np.concatenate([[e[1] for e in t] for t in trazas])
            bordes_dwell = np.quantile(dw_pool, np.arange(1, 8) / 8)
            variantes = {
                "centros_4_8_16": BORDES_PASO_NM,
                "terciles": np.quantile(sx_pool, [1 / 3, 2 / 3]),
            }
            print(f"[B1v2] {constructo} {cond}: {len(trazas)} trazas, "
                  f"{len(dw_pool)} eventos, "
                  f"{sum(len(t) - 1 for t in trazas)} pares", flush=True)
            for variante, bordes_paso in variantes.items():
                trazas_sim = []
                for t in trazas:
                    dw = np.array([e[0] for e in t])
                    sx = np.array([e[1] for e in t])
                    b = np.searchsorted(bordes_dwell, dw,
                                        side="right").astype(np.int8)
                    c = np.searchsorted(bordes_paso, sx,
                                        side="right").astype(np.int8)
                    trazas_sim.append((b * A.N_CLASES_PASO + c).astype(np.int8))
                for nulo, tipo in (("orden", 7), ("paridad", 8)):
                    rng = A.semilla_de(tipo, id_real, 0)
                    conteos8 = conteos_pares_variable(
                        trazas_sim, 8 * A.N_CLASES_PASO, rng,
                        paridad=(nulo == "paridad"))
                    for k in (4, 8):
                        canales = A.derivar_canales(conteos8, k)
                        for canal, (cj, Ax, Ay) in canales.items():
                            ims = A.estimar_im(cj, Ax, Ay)
                            for est, vals in ims.items():
                                r = A.test_permutacion(vals)
                                filas.append(dict(
                                    constructo=constructo, cond=cond,
                                    id_celda_real=id_real, variante=variante,
                                    nulo=nulo, canal=canal, k=k,
                                    estimador=est, n_trazas=len(trazas),
                                    n_pares=int(sum(len(t) - 1
                                                    for t in trazas)), **r))

    df = pd.DataFrame(filas)
    claves = ["variante", "nulo", "canal", "k", "estimador", "constructo"]
    df["p_holm"] = np.nan
    for _, idx in df.groupby(claves, sort=False).groups.items():
        df.loc[idx, "p_holm"] = holm(df.loc[idx, "p_perm"].to_numpy())
    df["confirma"] = df["rechaza"].astype(bool) & (df["p_holm"] < ALFA)
    df.to_csv(os.path.join(AQUI, "b1v2_resultados.csv"), index=False,
              float_format="%.6g")

    def tabla(variante, nulo, k, est="pt", canal="evento"):
        s = df[(df.variante == variante) & (df.nulo == nulo) & (df.k == k)
               & (df.estimador == est) & (df.canal == canal)]
        return s[["constructo", "cond", "n_trazas", "n_pares", "im_obs",
                  "umbral_p95", "p_perm", "p_holm", "confirma"]]

    prim = tabla("centros_4_8_16", "orden", 4)
    par = tabla("centros_4_8_16", "paridad", 4)
    print("\n=== PRIMARIO (evento, PT, k=4, nulo orden, centros 4/8/16) ===")
    print(prim.to_string(index=False))
    print("\n=== GUARD DE PARIDAD (nulo paridad) ===")
    print(par.to_string(index=False))
    print("\n=== SENSIBILIDAD k=8 (nulo orden) ===")
    print(tabla("centros_4_8_16", "orden", 8).to_string(index=False))
    print("\n=== SENSIBILIDAD terciles (k=4, nulo orden) ===")
    print(tabla("terciles", "orden", 4).to_string(index=False))

    if prim["confirma"].any():
        celdas = prim[prim.confirma][["constructo", "cond"]].values.tolist()
        sobrevive = par["confirma"].any()
        celdas_par = par[par.confirma][["constructo", "cond"]].values.tolist()
        print(f"\nVEREDICTO B1v2: CONFIRMA en {celdas}.")
        print("Guard de limping: " + (
            f"la señal SOBREVIVE al nulo de paridad en {celdas_par} — "
            f"memoria más allá de la alternancia par/impar."
            if sobrevive else
            "la señal DESAPARECE bajo el nulo de paridad — se reporta como "
            "REPLICACIÓN DE LIMPING (Asbury 2003), no como hallazgo."))
    else:
        print("\nVEREDICTO B1v2: REFUTA — ninguna condición rechaza con "
              "Holm en ningún constructo. (Acotado por el mapa de "
              "sensibilidad de A2': no exonera memorias débiles.)")


if __name__ == "__main__":
    main()
