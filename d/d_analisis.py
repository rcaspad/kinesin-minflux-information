# -*- coding: utf-8 -*-
"""
Etapa D — Atribución del residuo (congelación: D_CONGELACION.md).
Determinista: semilla_de(12, test_id, id_celda_real), raíz 20260812.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(AQUI, "..", "a2"))
sys.path.insert(0, os.path.join(AQUI, "..", "c_crudo"))
import a2_estimadores as A                       # noqa: E402
from c_kv_pipeline import segmentar_celda        # noqa: E402

BASE = os.path.join(AQUI, "..", "..", "analisis", "extracted",
                    "Repository for MINFLUX dissects the unimpeded walking "
                    "of kinesin-1", "Data repository", "KinesinDataFiles")
CELDAS = [("E215C", "10uM", 0, True), ("E215C", "100uM", 1, False),
          ("E215C", "1mM", 2, False), ("K28C", "10uM", 3, True),
          ("K28C", "100uM", 4, False), ("K28C", "1mM", 5, False),
          ("T324C", "10uM", 6, False), ("T324C", "1mM", 8, True)]
N_PERM = 1000
MISMO_SITIO, OTRO_SITIO = 4.0, 8.0
CORTE_NM = 12.0


def trazas_kv(constructo, cond):
    out = []
    for tr in segmentar_celda(os.path.join(BASE, constructo, cond)):
        ev = tr["eventos"]
        if len(ev) < 2:
            continue
        dw = np.array([e[0] for e in ev])
        ps = np.array([e[1] for e in ev])          # con signo
        pos = np.concatenate([[0.0], np.cumsum(ps)[:-1]])
        out.append((np.log(dw), ps, pos))
    return out


def d1_revisitas(trazas, rng):
    same, diff_por_estrato = [], {}
    for ldw, ps, pos in trazas:
        n = len(ldw)
        for i in range(n):
            for j in range(i + 3, n):
                d = abs(pos[i] - pos[j])
                sep = j - i
                if d < MISMO_SITIO:
                    same.append((id(ldw), sep, abs(ldw[j] - ldw[i]), True))
                elif d >= OTRO_SITIO:
                    diff_por_estrato.setdefault((id(ldw), sep), []).append(
                        abs(ldw[j] - ldw[i]))
    if len(same) < 30:
        return dict(n_same=len(same), T=np.nan, p=np.nan,
                    veredicto="INCONCLUSO (<30 pares)")
    # controles: mismos estratos (traza, separación) que los pares same
    pares = []                                    # (delta, es_same) por estrato
    for tid, sep, dsame, _ in same:
        ctrl = diff_por_estrato.get((tid, sep), [])
        if ctrl:
            pares.append((dsame, np.array(ctrl)))
    if len(pares) < 30:
        return dict(n_same=len(same), T=np.nan, p=np.nan,
                    veredicto="INCONCLUSO (<30 estratos con control)")

    def estadistico(eleccion):
        # eleccion[k]: índice del valor tratado como "same" en el estrato k
        vs, vd = [], []
        for k, (dsame, ctrl) in enumerate(pares):
            todos = np.concatenate([[dsame], ctrl])
            e = eleccion[k] % len(todos)
            vs.append(todos[e])
            vd.append(np.delete(todos, e).mean())
        return float(np.median(vs) - np.median(vd))

    T_obs = estadistico(np.zeros(len(pares), dtype=int))
    nulos = np.empty(N_PERM)
    for b in range(N_PERM):
        nulos[b] = estadistico(rng.integers(0, 10 ** 9, size=len(pares)))
    p = (1.0 + np.sum(nulos <= T_obs)) / (1.0 + N_PERM)   # H_rail: T<0
    return dict(n_same=len(same), n_estratos=len(pares), T=round(T_obs, 4),
                p=round(p, 4),
                veredicto="RAIL (sitio con identidad)" if p < 0.05
                else "sin identidad de sitio detectable")


def d2_decaimiento(trazas, rng):
    x_c, y_c, x_l, y_l, etiquetas, pares = [], [], [], [], [], []
    for ldw, ps, pos in trazas:
        m = ldw - ldw.mean()
        for k in range(len(ldw) - 1):
            cerca = abs(ps[k]) <= CORTE_NM
            pares.append((m[k], m[k + 1], cerca))
    if len(pares) < 100:
        return dict(n_pares=len(pares), D=np.nan, p=np.nan,
                    veredicto="INCONCLUSO")
    arr = np.array([(a, b) for a, b, _ in pares])
    lab = np.array([c for _, _, c in pares])

    def D_de(labels):
        rc = spearmanr(arr[labels, 0], arr[labels, 1]).statistic
        rl = spearmanr(arr[~labels, 0], arr[~labels, 1]).statistic
        return float(rc - rl)

    D_obs = D_de(lab)
    nulos = np.empty(N_PERM)
    for b in range(N_PERM):
        nulos[b] = D_de(rng.permutation(lab))
    p = (1.0 + np.sum(nulos >= D_obs)) / (1.0 + N_PERM)   # H_rail: D>0
    return dict(n_pares=len(pares), n_cerca=int(lab.sum()),
                D=round(D_obs, 4), p=round(p, 4),
                veredicto="decaimiento ESPACIAL (raíl)" if p < 0.05
                else "sin decaimiento espacial detectable")


def d3_rachas(trazas, rng):
    def rachas_de(lentos, ps):
        L_ev, L_nm = [], []
        k = 0
        while k < len(lentos):
            if lentos[k]:
                j = k
                while j + 1 < len(lentos) and lentos[j + 1]:
                    j += 1
                L_ev.append(j - k + 1)
                L_nm.append(float(np.abs(ps[k:j + 1]).sum()))
                k = j + 1
            else:
                k += 1
        return L_ev, L_nm

    obs_ev, obs_nm = [], []
    listas = []
    for ldw, ps, pos in trazas:
        lentos = ldw > np.median(ldw)
        le, ln = rachas_de(lentos, ps)
        obs_ev += le
        obs_nm += ln
        listas.append((lentos, ps))
    M_obs = float(np.mean(obs_ev))
    nulos = np.empty(N_PERM)
    for b in range(N_PERM):
        acc = []
        for lentos, ps in listas:
            acc += rachas_de(rng.permutation(lentos), ps)[0]
        nulos[b] = np.mean(acc)
    p = (1.0 + np.sum(nulos >= M_obs)) / (1.0 + N_PERM)
    return dict(n_rachas=len(obs_ev), L_ev_media=round(M_obs, 3),
                L_ev_nulo=round(float(nulos.mean()), 3),
                L_nm_mediana=round(float(np.median(obs_nm)), 1),
                p=round(p, 4),
                veredicto="rachas en EXCESO" if p < 0.05 else "sin exceso")


def d4_n356c():
    filas = []
    for dol in ["DOL1", "DOL2"]:
        ruta = os.path.join(BASE, "..", "KinesinDataFiles", "N356C", dol)
        ruta = os.path.normpath(os.path.join(BASE, "N356C", dol))
        if not os.path.isdir(ruta):
            filas.append((dol, "no existe", 0, 0))
            continue
        trs = segmentar_celda(ruta)
        n_ev = [len(t["eventos"]) for t in trs]
        filas.append((dol, len(trs), int(np.sum(n_ev)),
                      round(float(np.mean(n_ev)), 2) if n_ev else 0))
    return filas


def main():
    filas = []
    for constructo, cond, id_real, foco in CELDAS:
        trazas = trazas_kv(constructo, cond)
        r1 = d1_revisitas(trazas, A.semilla_de(12, 1, id_real))
        r2 = d2_decaimiento(trazas, A.semilla_de(12, 2, id_real))
        r3 = d3_rachas(trazas, A.semilla_de(12, 3, id_real))
        print(f"[D] {constructo} {cond}{' *FOCO*' if foco else ''}")
        print(f"    D-1 revisitas: {r1}")
        print(f"    D-2 decaimiento: {r2}")
        print(f"    D-3 rachas: {r3}", flush=True)
        filas.append(dict(constructo=constructo, cond=cond, foco=foco,
                          **{f"d1_{k}": v for k, v in r1.items()},
                          **{f"d2_{k}": v for k, v in r2.items()},
                          **{f"d3_{k}": v for k, v in r3.items()}))
    pd.DataFrame(filas).to_csv(os.path.join(AQUI, "d_resultados.csv"),
                               index=False)
    print("\n[D-4] Control estacionario N356C (eventos = fabricación "
          "del instrumento):")
    for dol, n_tr, ev_tot, ev_medio in d4_n356c():
        print(f"    {dol}: {n_tr} trazas, {ev_tot} eventos totales, "
              f"{ev_medio} eventos/traza")


if __name__ == "__main__":
    main()
