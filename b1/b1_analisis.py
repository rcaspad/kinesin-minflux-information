# -*- coding: utf-8 -*-
"""
B1 — Memoria más allá de la renovación en el caminar real de la kinesina-1
(prereg v1.0 §3-B1, congelación previa: B1_CONGELACION.md).

Determinista: semillas semilla_de(3|4, constructo_id, cond_id, tray_id).
Reutiliza los estimadores VALIDADOS en A2 (resultados/a2/a2_estimadores.py,
SHA dd784b2c…6eab9936) sin modificarlos.

Salidas: b1_por_traza.csv (por-instancia), b1_resumen.csv (por condición),
veredicto por stdout.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import binom

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(AQUI, "..", "a2"))
import a2_estimadores as A  # noqa: E402

DATA = os.path.join(AQUI, "..", "..", "analisis",
                    "extracted",
                    "Repository for MINFLUX dissects the unimpeded walking "
                    "of kinesin-1", "Data repository", "KinesinDataFiles")
COLS = ["stepx", "stepy", "tau", "sig_x", "sig_y", "phot_x", "phot_y",
        "end_flag", "transitions"]
CONSTRUCTOS = ["E215C", "K28C", "T324C"]
CONDICIONES = ["10uM", "100uM", "1mM"]
MIN_EVENTOS = 50           # prereg §4, congelado
MIN_ELEGIBLES = 10         # prereg §3-B1 (INCONCLUSO si < 10 en todas)
UMBRAL_FRACCION = 0.60     # prereg §3-B1
BORDES_PASO_NM = np.array([6.0, 12.0])   # centro más cercano de {4,8,16} nm
N_PERM = A.N_PERMUTACIONES               # 1000
ALFA = 0.05


def cargar_trazas(constructo: str, cond: str):
    """Trazas como listas de eventos (tau, |stepx|) en orden de fichero."""
    f = os.path.join(DATA, constructo, cond, "allsteps_reeval.xls")
    df = pd.read_excel(f, header=0)
    df.columns = COLS
    trazas, actual = [], []
    for fila in df.itertuples(index=False):
        if fila.end_flag == 1:
            trazas.append(actual)
            actual = []
        elif fila.tau > 0:
            actual.append((float(fila.tau), abs(float(fila.stepx))))
    if actual:
        trazas.append(actual)
    return trazas


def conteos_pares_matriz(M: np.ndarray, alfabeto: int) -> np.ndarray:
    """Conteos de pares consecutivos por fila de M (P, L). Devuelve (P, A²)."""
    P = M.shape[0]
    a2 = alfabeto * alfabeto
    cod = M[:, :-1].astype(np.int32) * alfabeto + M[:, 1:]
    desplaz = (np.arange(P, dtype=np.int32) * a2)[:, None]
    return np.bincount((cod + desplaz).ravel(),
                       minlength=P * a2).reshape(P, a2)


def permutaciones(sim: np.ndarray, rng, paridad: bool) -> np.ndarray:
    """(1+N_PERM, L): fila 0 = observado; resto, orden permutado.

    paridad=False: permutación libre del orden (nulo primario, prereg §5).
    paridad=True : permutación SOLO dentro de posiciones pares e impares
                   (nulo de paridad — guard de limping, prereg §5).
    """
    M = np.tile(sim, (1 + N_PERM, 1))
    if paridad:
        M[1:, 0::2] = rng.permuted(M[1:, 0::2], axis=1)
        M[1:, 1::2] = rng.permuted(M[1:, 1::2], axis=1)
    else:
        M[1:] = rng.permuted(M[1:], axis=1)
    return M


def analizar_celda(ci: int, constructo: str, di: int, cond: str,
                   variante_paso: str, filas: list):
    trazas = cargar_trazas(constructo, cond)
    elegibles = [(i, t) for i, t in enumerate(trazas) if len(t) >= MIN_EVENTOS]
    n_eleg = len(elegibles)
    print(f"[B1] {constructo} {cond} ({variante_paso}): {len(trazas)} trazas, "
          f"{n_eleg} elegibles (>= {MIN_EVENTOS} eventos)", flush=True)
    if n_eleg == 0:
        return 0

    dwells_pool = np.concatenate([[e[0] for e in t] for _, t in elegibles])
    bordes_dwell = np.quantile(dwells_pool, np.arange(1, 8) / 8)
    if variante_paso == "terciles":
        pasos_pool = np.concatenate([[e[1] for e in t] for _, t in elegibles])
        bordes_paso = np.quantile(pasos_pool, [1 / 3, 2 / 3])
    else:
        bordes_paso = BORDES_PASO_NM

    for tray_id, t in elegibles:
        dw = np.array([e[0] for e in t])
        sx = np.array([e[1] for e in t])
        bins = np.searchsorted(bordes_dwell, dw, side="right").astype(np.int8)
        clase = np.searchsorted(bordes_paso, sx, side="right").astype(np.int8)
        sim = (bins * A.N_CLASES_PASO + clase).astype(np.int8)

        for nulo, tipo_semilla in (("orden", 3), ("paridad", 4)):
            rng = A.semilla_de(tipo_semilla, ci, di, tray_id)
            M = permutaciones(sim, rng, paridad=(nulo == "paridad"))
            conteos8 = conteos_pares_matriz(M, 8 * A.N_CLASES_PASO)
            for k in (4, 8):
                canales = A.derivar_canales(conteos8, k)
                for canal, (cj, Ax, Ay) in canales.items():
                    ims = A.estimar_im(cj, Ax, Ay)
                    for est, vals in ims.items():
                        r = A.test_permutacion(vals)
                        filas.append(dict(
                            constructo=constructo, cond=cond,
                            variante_paso=variante_paso, nulo=nulo,
                            tray_id=tray_id, n_eventos=len(t),
                            canal=canal, k=k, estimador=est, **r))
    return n_eleg


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
    for variante in ("centros_4_8_16", "terciles"):
        for ci, constructo in enumerate(CONSTRUCTOS):
            for di, cond in enumerate(CONDICIONES):
                analizar_celda(ci, constructo, di, cond, variante, filas)

    if not filas:
        print(f"\nVEREDICTO B1: INCONCLUSO — 0 trazas elegibles "
              f"(>= {MIN_EVENTOS} eventos) en TODAS las celdas; el criterio "
              f"de muerte del prereg §3-B1 (N elegible < {MIN_ELEGIBLES} en "
              f"todas las condiciones) se cumple en su forma extrema. No se "
              f"computó ningún estadístico de información sobre datos reales.")
        return

    df = pd.DataFrame(filas)
    df.to_csv(os.path.join(AQUI, "b1_por_traza.csv"), index=False,
              float_format="%.6g")

    # Resumen por condición (todas las combinaciones registradas; la decisión
    # primaria es canal=evento, estimador=pt, k=4, variante centros, nulo orden)
    g = (df.groupby(["variante_paso", "nulo", "canal", "k", "estimador",
                     "constructo", "cond"], sort=False)
           .agg(n_eleg=("rechaza", "size"), n_rechaza=("rechaza", "sum"))
           .reset_index())
    g["fraccion"] = g["n_rechaza"] / g["n_eleg"]
    g["p_binomial"] = binom.sf(g["n_rechaza"] - 1, g["n_eleg"], ALFA)
    # Holm entre las 3 condiciones dentro de cada (todo lo demás, constructo)
    g["p_holm"] = np.nan
    claves = ["variante_paso", "nulo", "canal", "k", "estimador", "constructo"]
    for _, idx in g.groupby(claves, sort=False).groups.items():
        g.loc[idx, "p_holm"] = holm(g.loc[idx, "p_binomial"].to_numpy())
    g["confirma"] = (g["fraccion"] >= UMBRAL_FRACCION) & (g["p_holm"] < ALFA)
    g.to_csv(os.path.join(AQUI, "b1_resumen.csv"), index=False,
             float_format="%.6g")

    # Veredicto primario congelado
    prim = g[(g.variante_paso == "centros_4_8_16") & (g.nulo == "orden")
             & (g.canal == "evento") & (g.k == 4) & (g.estimador == "pt")]
    par = g[(g.variante_paso == "centros_4_8_16") & (g.nulo == "paridad")
            & (g.canal == "evento") & (g.k == 4) & (g.estimador == "pt")]
    print("\n=== PRIMARIO (evento, PT, k=4, nulo orden) ===")
    print(prim[["constructo", "cond", "n_eleg", "n_rechaza", "fraccion",
                "p_binomial", "p_holm", "confirma"]].to_string(index=False))
    print("\n=== GUARD DE PARIDAD (mismo estadístico, nulo paridad) ===")
    print(par[["constructo", "cond", "n_eleg", "n_rechaza", "fraccion",
               "p_binomial", "p_holm", "confirma"]].to_string(index=False))

    max_eleg = prim["n_eleg"].max() if len(prim) else 0
    if len(prim) == 0 or max_eleg < MIN_ELEGIBLES:
        print(f"\nVEREDICTO B1: INCONCLUSO — N elegible < {MIN_ELEGIBLES} "
              f"en todas las condiciones (máx = {max_eleg}).")
    elif prim["confirma"].any():
        conds_ok = prim[prim.confirma][["constructo", "cond"]].to_records(index=False)
        sobrevive = par["confirma"].any()
        print(f"\nVEREDICTO B1: CONFIRMA en {list(conds_ok)}.")
        print("Guard de limping: " + (
            "la señal SOBREVIVE al nulo de paridad — memoria más allá de la "
            "alternancia." if sobrevive else
            "la señal DESAPARECE bajo el nulo de paridad — se reporta como "
            "REPLICACIÓN DE LIMPING (Asbury 2003), no como hallazgo."))
    else:
        print("\nVEREDICTO B1: REFUTA — ninguna condición alcanza el criterio.")


if __name__ == "__main__":
    main()
