# -*- coding: utf-8 -*-
"""
Etapa E — Desempate de la identidad de sitio (congelación: E_CONGELACION.md).
E-1: N356C COM caminando (cabeza-vs-COM). E-2: prueba reina (trazas K28C
sobre el mismo microtúbulo, muestras activamente estabilizadas).
Determinista: semillas tipo 13-16, raíz 20260812.
"""
from __future__ import annotations

import glob
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(AQUI, "..", "a2"))
sys.path.insert(0, os.path.join(AQUI, "..", "c_crudo"))
import a2_estimadores as A                                   # noqa: E402
from c_kv_pipeline import (cargar_fichero_crudo, kv_bic,     # noqa: E402
                           fusionar_tramos, segmentar_celda,
                           proyectar_eje_principal)

BASE = os.path.join(AQUI, "..", "..", "analisis", "extracted",
                    "Repository for MINFLUX dissects the unimpeded walking "
                    "of kinesin-1", "Data repository", "KinesinDataFiles")
N_PERM = 1000
COLS = ["stepx", "stepy", "tau", "sig_x", "sig_y", "phot_x", "phot_y",
        "end_flag", "transitions"]
CELDAS_E1 = [("DOL1", "10uM", 9), ("DOL1", "1mM", 10), ("DOL2", "1mM", 11)]


def eventos_autores(dol, cond):
    f = os.path.join(BASE, "N356C", dol, cond, "allsteps_reeval.xls")
    df = pd.read_excel(f, header=0)
    df.columns = COLS
    return int(((df.end_flag == 0) & (df.tau > 0)).sum())


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


def im_dwell_pt_k4(conteos24):
    cj, Ax, Ay = A.derivar_canales(conteos24, 4)["dwell"]
    return A.estimar_im(cj, Ax, Ay)["pt"]


def d1_revisitas(trazas, rng):
    """Idéntico a d_analisis.d1_revisitas (copiado para autonomía del script)."""
    same, diff = [], {}
    for ldw, ps, pos in trazas:
        n = len(ldw)
        for i in range(n):
            for j in range(i + 3, n):
                d = abs(pos[i] - pos[j])
                if d < 4.0:
                    same.append((id(ldw), j - i, abs(ldw[j] - ldw[i])))
                elif d >= 8.0:
                    diff.setdefault((id(ldw), j - i), []).append(
                        abs(ldw[j] - ldw[i]))
    pares = [(s[2], np.array(diff[(s[0], s[1])]))
             for s in same if diff.get((s[0], s[1]))]
    if len(same) < 30 or len(pares) < 30:
        return dict(n_same=len(same), T=np.nan, p=np.nan,
                    veredicto="INCONCLUSO")

    def stat(eleccion):
        vs, vd = [], []
        for k, (dsame, ctrl) in enumerate(pares):
            todos = np.concatenate([[dsame], ctrl])
            e = eleccion[k] % len(todos)
            vs.append(todos[e])
            vd.append(np.delete(todos, e).mean())
        return float(np.median(vs) - np.median(vd))

    T = stat(np.zeros(len(pares), dtype=int))
    nulos = np.array([stat(rng.integers(0, 10 ** 9, size=len(pares)))
                      for _ in range(N_PERM)])
    p = (1.0 + np.sum(nulos <= T)) / (1.0 + N_PERM)
    return dict(n_same=len(same), T=round(T, 4), p=round(float(p), 4),
                veredicto="IDENTIDAD DE SITIO" if p < 0.05 else "nulo")


def e1():
    print("== E-1: N356C (COM, tallo) ==", flush=True)
    filas = []
    for dol, cond, id_c in CELDAS_E1:
        trs = segmentar_celda(os.path.join(BASE, "N356C", dol, cond))
        ev_kv = sum(len(t["eventos"]) for t in trs)
        ev_aut = eventos_autores(dol, cond)
        ratio = ev_kv / max(ev_aut, 1)
        comparable = 0.5 <= ratio <= 2.0
        print(f"[E-1] N356C {dol} {cond}: {len(trs)} trazas, {ev_kv} eventos "
              f"KV vs {ev_aut} autores (ratio {ratio:.2f}) → "
              f"{'COMPARABLE' if comparable else 'NO COMPARABLE'}", flush=True)
        if not comparable:
            filas.append(dict(dol=dol, cond=cond, comparable=False))
            continue
        trazas = []
        trazas_sim_src = [t["eventos"] for t in trs if len(t["eventos"]) >= 2]
        dw_pool = np.concatenate([[e[0] for e in ev] for ev in trazas_sim_src])
        bordes = np.quantile(dw_pool, np.arange(1, 8) / 8)
        sims = []
        for ev in trazas_sim_src:
            dw = np.array([e[0] for e in ev])
            ps = np.array([e[1] for e in ev])
            pos = np.concatenate([[0.0], np.cumsum(ps)[:-1]])
            trazas.append((np.log(dw), ps, pos))
            b = np.searchsorted(bordes, dw, side="right").astype(np.int8)
            c = np.searchsorted([6.0, 12.0], np.abs(ps)).astype(np.int8)
            sims.append((b * 3 + c).astype(np.int8))
        res = dict(dol=dol, cond=cond, comparable=True,
                   n_trazas=len(sims),
                   n_pares=int(sum(len(s) - 1 for s in sims)))
        for nulo, tipo, par in (("orden", 13, False), ("paridad", 14, True)):
            rng = A.semilla_de(tipo, id_c, 0)
            r = A.test_permutacion(im_dwell_pt_k4(
                conteos_permutados(sims, rng, par)))
            res[f"im_{nulo}"] = round(r["im_obs"], 5)
            res[f"p_{nulo}"] = round(r["p_perm"], 4)
        rev = d1_revisitas(trazas, A.semilla_de(15, id_c, 0))
        res.update({f"rev_{k}": v for k, v in rev.items()})
        # E-1c: periodicidad de dos pasos (estadístico D-2) en COM
        pares_e1c = []
        for ldw, ps, pos in trazas:
            m = ldw - ldw.mean()
            for k in range(len(ldw) - 1):
                pares_e1c.append((m[k], m[k + 1], abs(ps[k]) <= 12.0))
        arr = np.array([(a, b) for a, b, _ in pares_e1c])
        lab = np.array([c for _, _, c in pares_e1c])
        if len(arr) >= 100 and 30 <= lab.sum() <= len(lab) - 30:
            def D_de(labels):
                rc = spearmanr(arr[labels, 0], arr[labels, 1]).statistic
                rl = spearmanr(arr[~labels, 0], arr[~labels, 1]).statistic
                return float(rc - rl)
            D_obs = D_de(lab)
            rng = A.semilla_de(14, id_c, 1)
            nulos = np.array([D_de(rng.permutation(lab))
                              for _ in range(N_PERM)])
            p_bilateral = (1.0 + np.sum(np.abs(nulos) >= abs(D_obs))) / (
                1.0 + N_PERM)
            res["e1c_D"] = round(D_obs, 4)
            res["e1c_p"] = round(float(p_bilateral), 4)
        else:
            res["e1c_D"], res["e1c_p"] = np.nan, np.nan
        print(f"      E-1a IM dwell→dwell: {res['im_orden']} "
              f"(p_orden {res['p_orden']}, p_paridad {res['p_paridad']}); "
              f"E-1b revisitas: {rev}", flush=True)
        filas.append(res)
    return filas


def e2():
    print("\n== E-2: PRUEBA REINA (K28C, mismos microtúbulos) ==", flush=True)
    pares_x, pares_y, pares_lat = [], [], []
    eventos_por_roi = []
    datos_roi = []
    for f in sorted(glob.glob(os.path.join(BASE, "stabilized_microtubules",
                                           "*.txt"))):
        trs = cargar_fichero_crudo(f)
        if len(trs) < 2:
            continue
        X = np.concatenate([t["x"] for t in trs])
        Y = np.concatenate([t["y"] for t in trs])
        xc, yc = X - X.mean(), Y - Y.mean()
        ang = 0.5 * np.arctan2(2 * (xc * yc).sum(), (xc * xc - yc * yc).sum())
        u, v = np.array([np.cos(ang), np.sin(ang)]), np.array([-np.sin(ang),
                                                               np.cos(ang)])
        eventos = []           # (traza_id, axial, lateral, logdwell centrado)
        for tid, tr in enumerate(trs):
            z = tr["x"] * u[0] + tr["y"] * u[1]
            lat = tr["x"] * v[0] + tr["y"] * v[1]
            cortes = kv_bic(z)
            if len(cortes) < 3:
                continue
            bordes = np.concatenate([[0], cortes, [len(z)]])
            niveles = np.array([z[i:j].mean()
                                for i, j in zip(bordes[:-1], bordes[1:])])
            bordes, niveles = fusionar_tramos(bordes, niveles)
            if len(niveles) < 4:
                continue
            ldws, axs, lats = [], [], []
            for s in range(1, len(niveles) - 1):
                dwell = tr["t"][bordes[s + 1] - 1] - tr["t"][bordes[s]]
                if dwell <= 0:
                    continue
                ldws.append(np.log(dwell))
                axs.append(niveles[s])
                lats.append(float(lat[bordes[s]:bordes[s + 1]].mean()))
            if len(ldws) < 2:
                continue
            m = float(np.mean(ldws))
            eventos += [(tid, a, l, d - m)
                        for a, l, d in zip(axs, lats, ldws)]
        eventos_por_roi.append((os.path.basename(f), len(trs), len(eventos)))
        datos_roi.append(eventos)

    # pares cruzados entre trazas: evento vs media de otros en su sitio
    def emparejar(datos, rng=None):
        xs, ys, lats = [], [], []
        for eventos in datos:
            arr = eventos
            for k, (tid, a, l, d) in enumerate(arr):
                otros = [d2 for (t2, a2, l2, d2) in arr
                         if t2 != tid and abs(a2 - a) < 4.0]
                if otros:
                    lat_min = min(abs(l2 - l) for (t2, a2, l2, d2) in arr
                                  if t2 != tid and abs(a2 - a) < 4.0)
                    xs.append(d)
                    ys.append(float(np.mean(otros)))
                    lats.append(lat_min)
        return np.array(xs), np.array(ys), np.array(lats)

    xs, ys, lats = emparejar(datos_roi)
    print(f"[E-2] {len(datos_roi)} ROIs útiles; eventos/ROI: "
          f"{[e[2] for e in eventos_por_roi]}")
    print(f"[E-2] pares cruzados mismo-sitio: {len(xs)}", flush=True)
    if len(xs) < 50:
        print("[E-2] INCONCLUSO (<50 pares cruzados)")
        return dict(n_pares=len(xs), veredicto="INCONCLUSO")

    rho_obs = float(spearmanr(xs, ys).statistic)
    rng = A.semilla_de(16, 0, 0)
    nulos = np.empty(N_PERM)
    for b in range(N_PERM):
        perm = []
        for eventos in datos_roi:
            por_traza = {}
            for e in eventos:
                por_traza.setdefault(e[0], []).append(e)
            ev_p = []
            for tid, lista in por_traza.items():
                ds = rng.permutation([e[3] for e in lista])
                ev_p += [(tid, e[1], e[2], d) for e, d in zip(lista, ds)]
            perm.append(ev_p)
        xp, yp, _ = emparejar(perm)
        nulos[b] = spearmanr(xp, yp).statistic
    p = (1.0 + np.sum(nulos >= rho_obs)) / (1.0 + N_PERM)

    cerca = lats < 6.0
    rho_lado = (float(spearmanr(xs[cerca], ys[cerca]).statistic)
                if cerca.sum() >= 30 else np.nan)
    print(f"[E-2] rho = {rho_obs:.4f}  p = {p:.4f}  "
          f"(estrato mismo-lado n={int(cerca.sum())}, rho={rho_lado})")
    return dict(n_pares=len(xs), rho=round(rho_obs, 4), p=round(float(p), 4),
                n_mismo_lado=int(cerca.sum()),
                rho_mismo_lado=(round(rho_lado, 4)
                                if not np.isnan(rho_lado) else None),
                veredicto=("RAIL: sitios compartidos entre motores"
                           if p < 0.05 else
                           "sin identidad de sitio compartida entre motores"))


def main():
    filas_e1 = e1()
    res_e2 = e2()
    pd.DataFrame(filas_e1).to_csv(os.path.join(AQUI, "e1_resultados.csv"),
                                  index=False)
    pd.DataFrame([res_e2]).to_csv(os.path.join(AQUI, "e2_resultados.csv"),
                                  index=False)


if __name__ == "__main__":
    main()
