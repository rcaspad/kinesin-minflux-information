# -*- coding: utf-8 -*-
"""
A2″ — gate go/no-go del frente Wirth (prereg_wirth_v1.md §5, SHA
1d6ab5d5…2586a2de; implementación congelada en A2BIS_CONGELACION.md).

Determinista: semillas tipo 17 (datos), 18 (permutaciones), raíz 20260812.
"""
from __future__ import annotations

import os
import sys
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy.stats import binom

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(AQUI, "..", "a2"))
sys.path.insert(0, os.path.join(AQUI, "..", "c_crudo"))
import a2_estimadores as A                      # noqa: E402
from c_kv_pipeline import (cargar_fichero_crudo, proyectar_eje_principal,  # noqa: E402
                           kv_bic, fusionar_tramos)

BASE_WOLFF = os.path.join(AQUI, "..", "..", "analisis", "extracted",
                          "Repository for MINFLUX dissects the unimpeded "
                          "walking of kinesin-1", "Data repository",
                          "KinesinDataFiles")
BASE_WIRTH = os.path.join(AQUI, "..", "..", "..", "raw", "inbox",
                          "kinesin_minflux_wirth2024_zenodo10718784",
                          "Uncovering kinesin dynamics in neurites with "
                          "MINFLUX - Repository", "Data Repository")

# (dataset, constructo, cond, id_semilla) — A2BIS_CONGELACION.md
CELDAS = [
    ("wolff", "K28C", "10uM", 3), ("wolff", "K28C", "100uM", 4),
    ("wolff", "K28C", "1mM", 5), ("wolff", "T324C", "10uM", 6),
    ("wolff", "T324C", "1mM", 8),
    ("wirth", "K28C", "50uM", 20), ("wirth", "K28C", "500uM", 21),
    ("wirth", "K28C", "5mM", 22), ("wirth", "T324C", "50uM", 23),
    ("wirth", "T324C", "5mM", 24),
]
BORDES_PASO = np.array([6.0, 12.0])
N_REPLICAS = 500
N_PERM = A.N_PERMUTACIONES         # 1000
ALFA = 0.05


def dir_celda(dataset, constructo, cond):
    base = BASE_WOLFF if dataset == "wolff" else BASE_WIRTH
    return os.path.join(base, constructo, cond)


def cargar_marginales(dataset, constructo, cond):
    """(t, sigma) por traza real + pools (tau, stepx) de los autores."""
    d = dir_celda(dataset, constructo, cond)
    txts = sorted(f for f in os.listdir(d) if f.endswith(".txt"))
    assert txts, d
    rejillas = []
    for txt in txts:                 # varias fechas/muestras por celda (Wolff)
        for tr in cargar_fichero_crudo(os.path.join(d, txt)):
            if len(tr["t"]) < 20:
                continue
            z = proyectar_eje_principal(tr["x"], tr["y"])
            sigma = float(np.std(np.diff(z)) / np.sqrt(2.0))
            rejillas.append((tr["t"], sigma))
    xls = [f for f in os.listdir(d) if f.endswith("allsteps_reeval.xls")]
    assert len(xls) == 1, (d, xls)
    df = pd.read_excel(os.path.join(d, xls[0]), header=0)
    df.columns = ["stepx", "stepy", "tau", "sig_x", "sig_y", "phot_x",
                  "phot_y", "end_flag", "transitions"]
    tau = pd.to_numeric(df["tau"], errors="coerce").to_numpy(float)
    stepx = pd.to_numeric(df["stepx"], errors="coerce").to_numpy(float)
    ok = np.isfinite(tau) & np.isfinite(stepx) & (tau > 0)
    return rejillas, tau[ok], stepx[ok]


def generar_sustituta(t, sigma, dwell_pool, paso_pool, rng):
    """sustituta_de_traza de C-2 + registro de la verdad (tiempos, clases)."""
    T = t[-1] - t[0]
    dws, pss = [], []
    acum = 0.0
    while acum < T:
        d = float(dwell_pool[rng.integers(len(dwell_pool))])
        dws.append(d)
        pss.append(float(paso_pool[rng.integers(len(paso_pool))]))
        acum += d
    tiempos_cambio = t[0] + np.cumsum(dws)[:-1]
    niveles = np.concatenate([[0.0], np.cumsum(pss[:-1])])
    idx = np.searchsorted(tiempos_cambio, t, side="right")
    z = niveles[idx] + rng.normal(0.0, sigma, size=len(t))
    clases_verdad = np.searchsorted(BORDES_PASO, np.abs(pss[:-1]),
                                    side="right")
    return z, tiempos_cambio, clases_verdad


def segmentar_con_fronteras(t, z):
    """Cadena congelada; devuelve eventos interiores y TODAS las
    transiciones fusionadas (tiempo, paso estimado)."""
    cortes = kv_bic(z)
    if len(cortes) < 3:
        return [], []
    bordes = np.concatenate([[0], cortes, [len(z)]])
    niveles = np.array([z[i:j].mean() for i, j in zip(bordes[:-1], bordes[1:])])
    bordes, niveles = fusionar_tramos(bordes, niveles)
    if len(niveles) < 4:
        return [], []
    eventos = []
    for s in range(1, len(niveles) - 1):
        dwell = t[bordes[s + 1] - 1] - t[bordes[s]]
        paso = niveles[s + 1] - niveles[s]
        if dwell > 0:
            eventos.append((float(dwell), float(paso)))
    transiciones = [(float(t[bordes[s]]), float(niveles[s] - niveles[s - 1]))
                    for s in range(1, len(niveles))]
    return eventos, transiciones


def simbolos_de_eventos(trazas_eventos):
    dw_pool = np.concatenate([[e[0] for e in ev] for ev in trazas_eventos])
    bordes = np.quantile(dw_pool, np.arange(1, 8) / 8)
    salida = []
    for ev in trazas_eventos:
        dw = np.array([e[0] for e in ev])
        sx = np.abs(np.array([e[1] for e in ev]))
        b = np.searchsorted(bordes, dw, side="right").astype(np.int8)
        c = np.searchsorted(BORDES_PASO, sx, side="right").astype(np.int8)
        salida.append((b * A.N_CLASES_PASO + c).astype(np.int8))
    return salida


def conteos_permutados(trazas_sim, rng):
    a2 = 24 * 24
    salida = np.zeros((1 + N_PERM, a2), dtype=np.int64)
    for sim in trazas_sim:
        if len(sim) < 2:
            continue
        cod = sim[:-1].astype(np.int32) * 24 + sim[1:]
        salida[0] += np.bincount(cod, minlength=a2)
        M = np.tile(sim, (N_PERM, 1))
        M = rng.permuted(M, axis=1)
        codp = M[:, :-1].astype(np.int32) * 24 + M[:, 1:]
        desplaz = (np.arange(N_PERM, dtype=np.int32) * a2)[:, None]
        salida[1:] += np.bincount((codp + desplaz).ravel(),
                                  minlength=N_PERM * a2).reshape(N_PERM, a2)
    return salida


def im_dwell_pt_k4(conteos24):
    canales = A.derivar_canales(conteos24, 4)
    cj, Ax, Ay = canales["dwell"]
    return A.estimar_im(cj, Ax, Ay)["pt"]


def replica(args):
    """Una réplica de una celda: FP (0/1) + confusión 3x3 + cobertura."""
    id_celda, rep, rejillas, dwell_pool, paso_pool = args
    rng_d = A.semilla_de(17, id_celda, rep)
    ev_reps, confusion = [], np.zeros((3, 3), dtype=np.int64)
    n_verdad, n_detect = 0, 0
    for t, sigma in rejillas:
        z, t_cambio, clases_v = generar_sustituta(
            t, sigma, dwell_pool, paso_pool, rng_d)
        if len(t_cambio) == 0:
            continue
        n_verdad += len(t_cambio)
        eventos, transiciones = segmentar_con_fronteras(t, z)
        if len(eventos) >= 2:
            ev_reps.append(eventos)
        n_detect += len(transiciones)
        for t_det, paso_est in transiciones:
            j = int(np.argmin(np.abs(t_cambio - t_det)))
            c_est = int(np.searchsorted(BORDES_PASO, abs(paso_est),
                                        side="right"))
            confusion[clases_v[j], c_est] += 1
    fp = 0
    n_pares = 0
    if ev_reps:
        sims = simbolos_de_eventos(ev_reps)
        n_pares = sum(len(s) - 1 for s in sims if len(s) >= 2)
        rng_p = A.semilla_de(18, id_celda, rep)
        conteos = conteos_permutados(sims, rng_p)
        r = A.test_permutacion(im_dwell_pt_k4(conteos))
        fp = int(bool(r["rechaza"]))
    return id_celda, rep, fp, confusion, n_verdad, n_detect, n_pares


def main():
    lo = int(binom.ppf(0.025, N_REPLICAS, ALFA))
    hi = int(binom.ppf(0.975, N_REPLICAS, ALFA))
    print(f"[A2\"] banda binomial exacta 95%: [{lo}, {hi}] FP "
          f"de {N_REPLICAS} a alfa={ALFA}", flush=True)

    tareas, meta = [], {}
    for dataset, constructo, cond, idc in CELDAS:
        rejillas, tau, stepx = cargar_marginales(dataset, constructo, cond)
        meta[idc] = (dataset, constructo, cond)
        print(f"[A2\"] {dataset} {constructo} {cond}: {len(rejillas)} "
              f"rejillas, pools {len(tau)}", flush=True)
        for rep in range(N_REPLICAS):
            tareas.append((idc, rep, rejillas, tau, stepx))

    filas_fpr, filas_conf, filas_cob = [], [], []
    acc = {idc: dict(fp=0, conf=np.zeros((3, 3), dtype=np.int64),
                     nv=0, nd=0, pares=[]) for idc in meta}
    with Pool(processes=8) as pool:
        for idc, rep, fp, conf, nv, nd, npares in pool.imap_unordered(
                replica, tareas, chunksize=4):
            a = acc[idc]
            a["fp"] += fp
            a["conf"] += conf
            a["nv"] += nv
            a["nd"] += nd
            a["pares"].append(npares)

    print("", flush=True)
    for idc, a in sorted(acc.items()):
        dataset, constructo, cond = meta[idc]
        dentro = lo <= a["fp"] <= hi
        conf = a["conf"]
        tot = conf.sum()
        glob = float(np.trace(conf)) / tot if tot else np.nan
        por_clase = [float(conf[c, c]) / conf[c].sum() if conf[c].sum()
                     else np.nan for c in range(3)]
        pasa_clases = (glob >= 0.80) and all(x >= 2 / 3 for x in por_clase)
        filas_fpr.append(dict(
            dataset=dataset, constructo=constructo, cond=cond, id=idc,
            fp=a["fp"], banda_lo=lo, banda_hi=hi, fpr_ok=dentro,
            pares_p50=float(np.median(a["pares"]))))
        filas_conf.append(dict(
            dataset=dataset, constructo=constructo, cond=cond, id=idc,
            exact_global=glob, exact_c4=por_clase[0],
            exact_c8=por_clase[1], exact_c16=por_clase[2],
            resoluble=pasa_clases,
            **{f"conf_{i}{j}": int(conf[i, j])
               for i in range(3) for j in range(3)}))
        filas_cob.append(dict(
            dataset=dataset, constructo=constructo, cond=cond, id=idc,
            transiciones_verdad=int(a["nv"]),
            transiciones_detectadas=int(a["nd"]),
            cobertura=float(a["nd"]) / a["nv"] if a["nv"] else np.nan))
        print(f"[A2\"] {dataset:5s} {constructo:5s} {cond:5s}  "
              f"FP {a['fp']:3d} {'OK ' if dentro else 'FUERA'}  "
              f"exact {glob:.3f} clases "
              f"{por_clase[0]:.3f}/{por_clase[1]:.3f}/{por_clase[2]:.3f} "
              f"{'RESOLUBLE' if pasa_clases else 'NO'}", flush=True)

    pd.DataFrame(filas_fpr).to_csv(os.path.join(AQUI, "a2bis_fpr.csv"),
                                   index=False, float_format="%.6g")
    pd.DataFrame(filas_conf).to_csv(os.path.join(AQUI, "a2bis_confusion.csv"),
                                    index=False, float_format="%.6g")
    pd.DataFrame(filas_cob).to_csv(os.path.join(AQUI, "a2bis_cobertura.csv"),
                                   index=False, float_format="%.6g")

    fpr_all = all(f["fpr_ok"] for f in filas_fpr)
    res_wirth = all(f["resoluble"] for f in filas_conf
                    if f["dataset"] == "wirth")
    print(f"\n[A2\"] A2\"-1 (FPR, 10 celdas): "
          f"{'PASA' if fpr_all else 'FALLA'}", flush=True)
    print(f"[A2\"] A2\"-2 (resolubilidad, 5 celdas Wirth): "
          f"{'PASA' if res_wirth else 'FALLA'}", flush=True)
    print(f"[A2\"] GATE: {'GO' if (fpr_all and res_wirth) else 'NO-GO'}",
          flush=True)


if __name__ == "__main__":
    main()
