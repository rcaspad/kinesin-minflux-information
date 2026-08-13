# -*- coding: utf-8 -*-
"""
Etapa C — Contrastes C-1 (replicación con segmentador independiente) y
C-2 (nulo de artefacto calibrado). Congelación: C_CONGELACION.md +
C_CONGELACION_v2.md (fusión < 5 nm; gate re-aplicado: 8/9 comparables,
excluida T324C 100uM).

Determinista: semillas tipo 9 (perms orden), 10 (perms paridad),
11 (sustitutas), raíz 20260812.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(AQUI, "..", "a2"))
sys.path.insert(0, AQUI)
import a2_estimadores as A                      # noqa: E402
from c_kv_pipeline import (segmentar_celda, proyectar_eje_principal,  # noqa: E402
                           cargar_fichero_crudo, kv_bic, fusionar_tramos,
                           eventos_de_traza)

BASE = os.path.join(AQUI, "..", "..", "analisis", "extracted",
                    "Repository for MINFLUX dissects the unimpeded walking "
                    "of kinesin-1", "Data repository", "KinesinDataFiles")
CELDAS = [("E215C", "10uM", 0), ("E215C", "100uM", 1), ("E215C", "1mM", 2),
          ("K28C", "10uM", 3), ("K28C", "100uM", 4), ("K28C", "1mM", 5),
          ("T324C", "10uM", 6), ("T324C", "1mM", 8)]   # sin T324C 100uM (gate)
BORDES_PASO = np.array([6.0, 12.0])
N_PERM = A.N_PERMUTACIONES
N_SUSTITUTAS = 200
ALFA = 0.05


def simbolos_de_eventos(trazas_eventos, k8_bordes=None):
    """Símbolos 24-arios por traza; bordes de dwell agrupados de la celda."""
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


def conteos_permutados(trazas_sim, rng, paridad):
    a2 = 24 * 24
    salida = np.zeros((1 + N_PERM, a2), dtype=np.int64)
    for sim in trazas_sim:
        if len(sim) < 2:
            continue
        cod = sim[:-1].astype(np.int32) * 24 + sim[1:]
        salida[0] += np.bincount(cod, minlength=a2)
        M = np.tile(sim, (N_PERM, 1))
        if paridad:
            M[:, 0::2] = rng.permuted(M[:, 0::2], axis=1)
            M[:, 1::2] = rng.permuted(M[:, 1::2], axis=1)
        else:
            M = rng.permuted(M, axis=1)
        codp = M[:, :-1].astype(np.int32) * 24 + M[:, 1:]
        desplaz = (np.arange(N_PERM, dtype=np.int32) * a2)[:, None]
        salida[1:] += np.bincount((codp + desplaz).ravel(),
                                  minlength=N_PERM * a2).reshape(N_PERM, a2)
    return salida


def im_dwell_pt_k4(conteos24) -> np.ndarray:
    canales = A.derivar_canales(conteos24, 4)
    cj, Ax, Ay = canales["dwell"]
    return A.estimar_im(cj, Ax, Ay)["pt"]


def corr_lag1_dwells(trazas_eventos):
    cs = []
    for ev in trazas_eventos:
        dw = np.array([e[0] for e in ev])
        if len(dw) >= 3 and dw.std() > 0:
            cs.append(np.corrcoef(dw[:-1], dw[1:])[0, 1])
    return float(np.nanmedian(cs)) if cs else np.nan


def holm(p):
    orden = np.argsort(p)
    m = len(p)
    adj = np.empty(m)
    cmax = 0.0
    for r, i in enumerate(orden):
        cmax = max(cmax, min(1.0, (m - r) * p[i]))
        adj[i] = cmax
    return adj


def sustituta_de_traza(t, sigma, dwell_pool, paso_pool, rng):
    """Trayectoria de renovación SIN memoria en la rejilla temporal real."""
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
    return z


def main():
    filas_c1, filas_c2 = [], []
    for constructo, cond, id_real in CELDAS:
        trazas = segmentar_celda(os.path.join(BASE, constructo, cond))
        trazas_ev = [tr["eventos"] for tr in trazas if len(tr["eventos"]) >= 2]
        n_pares = sum(len(ev) - 1 for ev in trazas_ev)
        sims = simbolos_de_eventos(trazas_ev)
        print(f"[C] {constructo} {cond}: {len(trazas_ev)} trazas KV, "
              f"{n_pares} pares", flush=True)

        # ---- C-1: replicación con permutaciones
        res = {}
        for nulo, tipo, par in (("orden", 9, False), ("paridad", 10, True)):
            rng = A.semilla_de(tipo, id_real, 0)
            conteos = conteos_permutados(sims, rng, par)
            vals = im_dwell_pt_k4(conteos)
            r = A.test_permutacion(vals)
            res[nulo] = r
            filas_c1.append(dict(constructo=constructo, cond=cond,
                                 nulo=nulo, n_pares=n_pares, **r))
        im_real = res["orden"]["im_obs"]

        # ---- C-2: nulo de artefacto calibrado
        dwell_pool = np.concatenate([[e[0] for e in ev] for ev in trazas_ev])
        paso_pool = np.concatenate([[e[1] for e in ev] for ev in trazas_ev])
        ims_surr, corr_surr, pares_surr = [], [], []
        for rep in range(N_SUSTITUTAS):
            rng = A.semilla_de(11, id_real, rep)
            ev_surr = []
            for tr in trazas:
                t, z = tr["t"], tr["z"]
                if len(t) < 20:
                    continue
                sigma = np.std(np.diff(z)) / np.sqrt(2.0)
                zs = sustituta_de_traza(t, sigma, dwell_pool, paso_pool, rng)
                ev = eventos_de_traza(t, zs)
                if len(ev) >= 2:
                    ev_surr.append(ev)
            if not ev_surr:
                continue
            sims_s = simbolos_de_eventos(ev_surr)
            conteos = np.zeros((1, 24 * 24), dtype=np.int64)
            for sim in sims_s:
                if len(sim) >= 2:
                    cod = sim[:-1].astype(np.int32) * 24 + sim[1:]
                    conteos[0] += np.bincount(cod, minlength=24 * 24)
            ims_surr.append(float(im_dwell_pt_k4(conteos)[0]))
            corr_surr.append(corr_lag1_dwells(ev_surr))
            pares_surr.append(sum(len(e) - 1 for e in ev_surr))

        ims_surr = np.array(ims_surr)
        p_art = (1.0 + np.sum(ims_surr >= im_real)) / (1.0 + len(ims_surr))
        filas_c2.append(dict(
            constructo=constructo, cond=cond, im_real=im_real,
            im_surr_p50=float(np.median(ims_surr)),
            im_surr_p95=float(np.percentile(ims_surr, 95)),
            p_artefacto=p_art,
            corr_lag1_real=corr_lag1_dwells(trazas_ev),
            corr_lag1_surr_p50=float(np.nanmedian(corr_surr)),
            n_pares_real=n_pares,
            n_pares_surr_p50=float(np.median(pares_surr))))
        print(f"    C-2: IM_real {im_real:.4f} vs surr p50 "
              f"{np.median(ims_surr):.4f} p95 {np.percentile(ims_surr,95):.4f}"
              f"  p_artefacto {p_art:.4f}", flush=True)

    c1 = pd.DataFrame(filas_c1)
    c1["p_holm"] = np.nan
    for _, idx in c1.groupby(["nulo", "constructo"], sort=False).groups.items():
        c1.loc[idx, "p_holm"] = holm(c1.loc[idx, "p_perm"].to_numpy())
    c1["confirma"] = c1["rechaza"].astype(bool) & (c1["p_holm"] < ALFA)
    c2 = pd.DataFrame(filas_c2)
    c2["artefacto_descartado"] = c2["p_artefacto"] < ALFA
    c1.to_csv(os.path.join(AQUI, "c1_replicacion.csv"), index=False,
              float_format="%.6g")
    c2.to_csv(os.path.join(AQUI, "c2_nulo_artefacto.csv"), index=False,
              float_format="%.6g")

    print("\n=== C-1 (dwell→dwell, PT k=4, segmentador KV) ===")
    print(c1[c1.nulo == "orden"][["constructo", "cond", "n_pares", "im_obs",
                                  "umbral_p95", "p_perm", "p_holm",
                                  "confirma"]].to_string(index=False))
    print("\n=== C-1 bajo nulo de paridad ===")
    print(c1[c1.nulo == "paridad"][["constructo", "cond", "p_perm", "p_holm",
                                    "confirma"]].to_string(index=False))
    print("\n=== C-2 (nulo de artefacto calibrado) ===")
    print(c2.to_string(index=False))


if __name__ == "__main__":
    main()
