# -*- coding: utf-8 -*-
"""
A2′ — Revalidación del estimador con la rejilla de LARGOS REALES
(prereg v2 §2, SHA e4aebb2a6511e291e8a988ad8f237aece1643352cdb64229d47363506b756842).

Solo datos sintéticos. Las longitudes reales entran como multiconjunto
empírico congelado (longitudes_reales.csv, SHA ee4afd8e…131957a8).

Reutiliza SIN CAMBIOS los estimadores validados de A2
(resultados/a2/a2_estimadores.py, SHA dd784b2c…6eab9936): generador,
simbolización, canales, estimadores y test de permutación.

Salidas: bloques_a2p/<id_sint>_<id_cond>.npy (checkpoint atómico),
a2p_por_celda.csv, a2p_minima_detectable.csv, veredicto por stdout.
"""
from __future__ import annotations

import itertools
import multiprocessing as mp
import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = "1"

import numpy as np
import pandas as pd

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(AQUI, "..", "a2"))
import a2_estimadores as A  # noqa: E402

DWELL_CANONICO = 25.0        # escala pura: degeneración exacta probada en A2
FORMAS_GAMMA = [1.0, 2.0]
MEZCLAS_PASO = {
    "equilibrada":   (1 / 3, 1 / 3, 1 / 3),
    "dominante_8nm": (0.15, 0.70, 0.15),
    "sesgada_4nm":   (0.60, 0.30, 0.10),
}
CONSTRUCTOS = ["E215C", "K28C", "T324C"]
CONDICIONES_ATP = ["10uM", "100uM", "1mM"]
DELTAS = [0.05, 0.10, 0.20, 0.35, 0.50]
RHOS = [0.05, 0.10, 0.20, 0.35, 0.50]
N_REPLICAS = 500
N_PERM = A.N_PERMUTACIONES
ESTIMADORES = ["plugin", "mm", "pt", "nsb"]
CANALES = ["evento", "dwell", "paso", "dwell_paso"]
DIR_BLOQUES = os.path.join(AQUI, "bloques_a2p")
MIN_EVENTOS = 2              # prereg v2 §0.2: mínimo que forma un par


def cargar_longitudes():
    """Multiconjunto de longitudes reales por celda (trazas >= 2 eventos)."""
    df = pd.read_csv(os.path.join(AQUI, "longitudes_reales.csv"))
    df = df[df.n_eventos >= MIN_EVENTOS]
    largos = {}
    for (c, d), g in df.groupby(["constructo", "cond"], sort=False):
        largos[(c, d)] = np.sort(g.n_eventos.to_numpy())[::-1].copy()
    return largos


LARGOS_REALES = cargar_longitudes()


def construir_celdas():
    """54 celdas sintéticas: 9 reales x 2 formas x 3 mezclas.

    id_sint = id_celda_real*6 + id_forma*3 + id_mezcla (prereg v2 §2).
    """
    celdas = []
    for (ci, constructo), (di, cond) in itertools.product(
            enumerate(CONSTRUCTOS), enumerate(CONDICIONES_ATP)):
        for fi, forma in enumerate(FORMAS_GAMMA):
            for mi, mezcla in enumerate(MEZCLAS_PASO):
                id_real = ci * 3 + di
                celdas.append(dict(
                    id_sint=id_real * 6 + fi * 3 + mi,
                    id_celda_real=id_real, constructo=constructo, cond=cond,
                    forma_gamma=forma, mezcla=mezcla))
    return celdas


def construir_condiciones():
    conds = [dict(id_cond=0, modo="H0", delta=0.0, rho=0.0)]
    j = 1
    for d in DELTAS:
        conds.append(dict(id_cond=j, modo="limping", delta=d, rho=0.0)); j += 1
    for r in RHOS:
        conds.append(dict(id_cond=j, modo="ar1", delta=0.0, rho=r)); j += 1
    return conds


def conteos_pares_variable(trazas_sim: list[np.ndarray], alfabeto: int,
                           rng: np.random.Generator) -> np.ndarray:
    """Conteos agrupados de pares con permutación DENTRO de cada traza.

    Devuelve (1 + N_PERM, alfabeto²): fila 0 = observado; filas 1..N_PERM =
    permutaciones del orden de eventos de cada traza (marginales por traza
    preservados). Las trazas se recorren en orden fijo (determinismo).
    """
    a2 = alfabeto * alfabeto
    salida = np.zeros((1 + N_PERM, a2), dtype=np.int64)
    for sim in trazas_sim:
        cod = sim[:-1].astype(np.int32) * alfabeto + sim[1:]
        salida[0] += np.bincount(cod, minlength=a2)
        M = np.tile(sim, (N_PERM, 1))
        M = rng.permuted(M, axis=1)
        codp = M[:, :-1].astype(np.int32) * alfabeto + M[:, 1:]
        desplaz = (np.arange(N_PERM, dtype=np.int32) * a2)[:, None]
        salida[1:] += np.bincount((codp + desplaz).ravel(),
                                  minlength=N_PERM * a2).reshape(N_PERM, a2)
    return salida


def _tarea(arg):
    celda, cond = arg
    ruta = os.path.join(DIR_BLOQUES,
                        f"{celda['id_sint']:04d}_{cond['id_cond']:02d}.npy")
    if os.path.exists(ruta):
        return ruta

    largos = LARGOS_REALES[(celda["constructo"], celda["cond"])]
    n_tray, l_max = len(largos), int(largos.max())
    pesos = np.array(MEZCLAS_PASO[celda["mezcla"]])
    filas = []
    for rep in range(N_REPLICAS):
        rng_dat = A.semilla_de(5, celda["id_sint"], cond["id_cond"], rep)
        dwells, clases = A.generar_instancia(
            n_tray=n_tray, largo=l_max, dwell_medio_ms=DWELL_CANONICO,
            forma_gamma=celda["forma_gamma"], pesos_paso=pesos,
            modo_memoria=cond["modo"], delta=cond["delta"], rho=cond["rho"],
            rng=rng_dat)
        # truncado a la longitud real de cada traza (prefijo estacionario)
        dw_pool = np.concatenate([dwells[i, :largos[i]] for i in range(n_tray)])
        bordes = np.quantile(dw_pool, np.arange(1, 8) / 8)
        trazas_sim = []
        for i in range(n_tray):
            b = np.searchsorted(bordes, dwells[i, :largos[i]],
                                side="right").astype(np.int8)
            trazas_sim.append((b * A.N_CLASES_PASO
                               + clases[i, :largos[i]]).astype(np.int8))

        rng_perm = A.semilla_de(6, celda["id_sint"], cond["id_cond"], rep)
        conteos8 = conteos_pares_variable(trazas_sim, 8 * A.N_CLASES_PASO,
                                          rng_perm)
        for k in (4, 8):
            canales = A.derivar_canales(conteos8, k)
            for canal, (cj, Ax, Ay) in canales.items():
                ims = A.estimar_im(cj, Ax, Ay)
                for est, vals in ims.items():
                    r = A.test_permutacion(vals)
                    filas.append((celda["id_sint"], cond["id_cond"], rep,
                                  CANALES.index(canal), k,
                                  ESTIMADORES.index(est),
                                  r["im_obs"], r["umbral_p95"],
                                  r["p_perm"], float(r["rechaza"])))
    arr = np.array(filas, dtype=np.float64)
    tmp = ruta + f".tmp{os.getpid()}"
    np.save(tmp, arr)
    os.replace(tmp + ".npy", ruta)
    return ruta


COLUMNAS = ["id_sint", "id_cond", "replica", "canal_id", "k", "estimador_id",
            "im_obs", "umbral_p95", "p_perm", "rechaza"]


def ejecutar(n_procesos: int):
    celdas = construir_celdas()
    conds = construir_condiciones()
    tareas = [(c, cond) for c in celdas for cond in conds]
    print(f"[A2'] {len(celdas)} celdas x {len(conds)} condiciones = "
          f"{len(tareas)} tareas; {n_procesos} procesos", flush=True)
    os.makedirs(DIR_BLOQUES, exist_ok=True)
    t0 = time.time()
    rutas, hechas = [], 0
    with mp.get_context("fork").Pool(n_procesos) as ex:
        for ruta in ex.imap_unordered(_tarea, tareas, chunksize=1):
            rutas.append(ruta)
            hechas += 1
            if hechas % 50 == 0 or hechas == len(tareas):
                el = time.time() - t0
                print(f"[A2'] {hechas}/{len(tareas)}  {el/60:.1f} min  "
                      f"(ETA {el/hechas*(len(tareas)-hechas)/60:.1f} min)",
                      flush=True)
    return sorted(rutas)


def agregar(rutas):
    arr = np.concatenate([np.load(r) for r in rutas], axis=0)
    df = pd.DataFrame(arr, columns=COLUMNAS)
    for c in ("id_sint", "id_cond", "replica", "k", "rechaza"):
        df[c] = df[c].astype(np.int32)
    df["canal"] = pd.Categorical.from_codes(df.pop("canal_id").astype(int), CANALES)
    df["estimador"] = pd.Categorical.from_codes(df.pop("estimador_id").astype(int), ESTIMADORES)
    df = (df.merge(pd.DataFrame(construir_celdas()), on="id_sint")
            .merge(pd.DataFrame(construir_condiciones()), on="id_cond"))

    g = (df.groupby(["id_sint", "id_celda_real", "constructo", "cond",
                     "forma_gamma", "mezcla", "id_cond", "modo", "delta",
                     "rho", "canal", "k", "estimador"], observed=True)
           .agg(n_replicas=("rechaza", "size"), tasa=("rechaza", "mean"))
           .reset_index())
    g.to_csv(os.path.join(AQUI, "a2p_por_celda.csv"), index=False,
             float_format="%.6g")

    # mínima detectable por celda (PT; potencia >= 0.80)
    pot = g[(g.modo != "H0") & (g.estimador == "pt")]
    filas = []
    for (ids, k, modo), sub in pot.groupby(["id_sint", "k", "modo"],
                                           observed=True):
        mag = sub["delta"] + sub["rho"]
        ok = sub[sub.tasa >= 0.80]
        filas.append(dict(
            id_sint=ids, k=k, modo=modo,
            constructo=sub.constructo.iloc[0], cond=sub.cond.iloc[0],
            forma_gamma=sub.forma_gamma.iloc[0], mezcla=sub.mezcla.iloc[0],
            min_detectable=(ok["delta"] + ok["rho"]).min() if len(ok) else np.nan,
            potencia_max=sub.tasa.max(), mag_max=mag.max()))
    md = pd.DataFrame(filas)
    md.to_csv(os.path.join(AQUI, "a2p_minima_detectable.csv"), index=False,
              float_format="%.6g")
    return g, md


def veredicto(g, md):
    from scipy.stats import binomtest
    print("\n=== A2' criterio (i): FPR bajo H0, canal evento, PT ===")
    h0 = g[(g.modo == "H0") & (g.canal == "evento") & (g.estimador == "pt")]
    for k in (4, 8):
        sub = h0[h0.k == k]
        pvals = np.array([binomtest(int(round(t * n)), int(n), 0.05,
                                    alternative="greater").pvalue
                          for t, n in zip(sub.tasa, sub.n_replicas)])
        orden = np.argsort(pvals)
        m = len(pvals)
        bh_sig = 0
        for rank, i in enumerate(np.sort(pvals)):
            if i <= 0.05 * (rank + 1) / m:
                bh_sig = rank + 1
        print(f"  k={k}: {m} celdas, FPR media {sub.tasa.mean():.4f}, "
              f"max {sub.tasa.max():.4f}, literal>0.05: {(sub.tasa > 0.05).sum()}, "
              f"significativas BH: {bh_sig}")
    crit_i = True
    for k in (4, 8):
        sub = h0[h0.k == k]
        pvals = np.array([binomtest(int(round(t * n)), int(n), 0.05,
                                    alternative="greater").pvalue
                          for t, n in zip(sub.tasa, sub.n_replicas)])
        psort = np.sort(pvals)
        m = len(psort)
        if any(psort[r] <= 0.05 * (r + 1) / m for r in range(m)):
            crit_i = False

    print("\n=== A2' criterio (ii): potencia PT k=4 en magnitudes plausibles ===")
    pot = g[(g.modo != "H0") & (g.canal == "evento")
            & (g.estimador == "pt") & (g.k == 4)]
    plaus = pot[((pot.modo == "limping") & (pot.delta <= 0.35))
                | ((pot.modo == "ar1") & (pot.rho <= 0.20))]
    crit_ii_por_constructo = {}
    for constructo in CONSTRUCTOS:
        sub = plaus[plaus.constructo == constructo]
        alcanza = sub[sub.tasa >= 0.80]
        crit_ii_por_constructo[constructo] = len(alcanza) > 0
        detalle = (alcanza.groupby(["cond", "modo"], observed=True)
                   .apply(lambda s: (s.delta + s.rho).min(), include_groups=False)
                   .to_dict() if len(alcanza) else {})
        print(f"  {constructo}: {'ALCANZA' if len(alcanza) else 'NO ALCANZA'} "
              f"potencia 0.80 en magnitud plausible; detalle {detalle}")
    crit_ii = all(crit_ii_por_constructo.values())

    if crit_i and crit_ii:
        print("\nVEREDICTO A2': PASA — calibración limpia y potencia "
              "alcanzable en magnitudes plausibles en los 3 constructos. "
              "B1v2 queda habilitada A FALTA DE LA FIRMA de Raúl (prereg v2 §5).")
    elif not crit_i:
        print("\nVEREDICTO A2': NO MEDIBLE — la calibración falla a los "
              "largos reales. El frente F1-B sobre este dataset SE CIERRA.")
    else:
        faltan = [c for c, ok in crit_ii_por_constructo.items() if not ok]
        parcial = [c for c, ok in crit_ii_por_constructo.items() if ok]
        print(f"\nVEREDICTO A2': criterio (ii) NO se cumple en TODOS los "
              f"constructos (fallan: {faltan}; alcanzan: {parcial}). Según "
              f"prereg v2 §2, NO MEDIBLE exige fallo en TODOS: "
              f"{'NO MEDIBLE — SE CIERRA' if not parcial else 'PASA PARCIAL — decisión documentada del gobernador sobre restringir B1v2 a los constructos que alcanzan'}.")


if __name__ == "__main__":
    n_proc = int(sys.argv[1]) if len(sys.argv) > 1 else max(1, (os.cpu_count() or 2) - 2)
    rutas = ejecutar(n_proc)
    g, md = agregar(rutas)
    veredicto(g, md)
