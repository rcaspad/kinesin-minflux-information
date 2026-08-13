# -*- coding: utf-8 -*-
"""
A2 — Barrido de la rejilla de regímenes (prereg v1.0 §4, criterios congelados).

Ejecuta:
  FASE H0       : tasa de falsos positivos (FPR) en la rejilla completa.
  FASE MEMORIA  : potencia frente a limping (delta) y AR(1) (rho).

Salidas (todas exportables, sin datos reales):
  a2_por_instancia.csv.gz  — una fila por (celda x réplica x canal x k x estimador)
  a2_por_celda.csv         — agregado por celda x canal x k x estimador (FPR/potencia)
  a2_minima_detectable.csv — magnitud mínima detectable con potencia >= 0.80 por celda

Semilla raíz 20260812; toda instancia queda determinada por
(SEMILLA_RAIZ, tipo, id_celda, id_cond, replica) — ver a2_estimadores.semilla_de.
"""

from __future__ import annotations

import itertools
import multiprocessing as mp
import os
import sys
import time

# Un solo hilo BLAS por proceso: el paralelismo de este barrido es por tareas,
# y dejar que cada uno de los 9 procesos abra 8 hilos sobresuscribe la máquina
# (carga observada 26 sobre 10 núcleos) y frena el conjunto.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = "1"

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import a2_estimadores as A

# --------------------------------------------------------------------------
# Rejilla de regímenes (prereg §4.2) — rangos plausibles de literatura,
# NO calibrados con datos reales.
# --------------------------------------------------------------------------

DWELL_MEDIO_MS = [10.0, 25.0, 50.0, 100.0]
DWELL_CANONICO = 25.0        # representante calculado del eje degenerado
FORMAS_GAMMA = [1.0, 2.0]
MEZCLAS_PASO = {                       # pesos sobre {4, 8, 16} nm
    "equilibrada":   (1 / 3, 1 / 3, 1 / 3),
    "dominante_8nm": (0.15, 0.70, 0.15),
    "sesgada_4nm":   (0.60, 0.30, 0.10),
}
LARGOS = [50, 100, 200]                # eventos por trayectoria
N_TRAYECTORIAS = [50, 100, 300]

DELTAS = [0.05, 0.10, 0.20, 0.35, 0.50]     # limping
RHOS = [0.05, 0.10, 0.20, 0.35, 0.50]       # AR(1)

# Al calcular sólo las 54 celdas distintas (ver construir_celdas_distintas; la
# degeneración que las reduce está documentada en construir_celdas) el presupuesto
# permite 500 réplicas: el error estándar de una FPR nominal de 0.05 baja de
# 0.015 (200 rep) a 0.010 (500 rep), lo que importa porque el criterio congelado
# se decide exactamente en 0.05.
N_REPLICAS_H0 = 500
N_REPLICAS_MEM = 500

ESTIMADORES = ["plugin", "mm", "pt", "nsb"]
CANALES = ["evento", "dwell", "paso", "dwell_paso"]
CANAL_PRIMARIO = "evento"


def construir_celdas():
    """Rejilla completa de regímenes (216 celdas). id_celda es determinista.

    IMPORTANTE — degeneración exacta del eje `dwell_medio_ms`:
    `dwell_medio_ms` es un parámetro de ESCALA puro de la distribución gamma, y
    la discretización por cuantiles (a2_estimadores.simbolizar) es invariante
    bajo cualquier transformación estrictamente creciente. Por tanto dos celdas
    que sólo difieren en `dwell_medio_ms` producen la MISMA secuencia de
    símbolos y resultados idénticos bit a bit (verificado sobre 216 pares de
    comparación: max|diferencia| = 0). Las 216 celdas del prereg contienen
    54 celdas distintas, cada una repetida en las 4 escalas de dwell.

    Se calculan las 54 distintas y se expanden a las 216 (columna
    `escala_dwell_degenerada`), de modo que el reporte por-celda cubre la
    rejilla completa del prereg sin gastar 4x cómputo en réplicas exactas.
    """
    celdas = []
    for i, (dm, fg, mz, lg, nt) in enumerate(itertools.product(
            DWELL_MEDIO_MS, FORMAS_GAMMA, sorted(MEZCLAS_PASO), LARGOS, N_TRAYECTORIAS)):
        celdas.append(dict(id_celda=i, dwell_medio_ms=dm, forma_gamma=fg,
                           mezcla=mz, largo=lg, n_tray=nt))
    return celdas


def construir_celdas_distintas():
    """Las 54 celdas efectivamente distintas (dwell_medio_ms fijado al canónico).

    id_celda es el índice de la celda de la rejilla completa con
    dwell_medio_ms = DWELL_CANONICO, para que las semillas derivadas sean las
    mismas que usaría el barrido de 216 celdas.
    """
    return [c for c in construir_celdas() if c["dwell_medio_ms"] == DWELL_CANONICO]


def expandir_a_rejilla_completa(df: pd.DataFrame) -> pd.DataFrame:
    """Replica los resultados de las 54 celdas distintas a las 216 del prereg."""
    completa = pd.DataFrame(construir_celdas())
    llaves = ["forma_gamma", "mezcla", "largo", "n_tray"]
    base = df.drop(columns=["id_celda", "dwell_medio_ms"])
    out = completa.merge(base, on=llaves, how="inner")
    out["escala_dwell_degenerada"] = out["dwell_medio_ms"] != DWELL_CANONICO
    return out


def construir_condiciones():
    """Condiciones de memoria. id_cond=0 es H0; el resto inyecta memoria."""
    conds = [dict(id_cond=0, modo="H0", delta=0.0, rho=0.0)]
    j = 1
    for d in DELTAS:
        conds.append(dict(id_cond=j, modo="limping", delta=d, rho=0.0)); j += 1
    for r in RHOS:
        conds.append(dict(id_cond=j, modo="ar1", delta=0.0, rho=r)); j += 1
    return conds


# --------------------------------------------------------------------------
# Ejecución
# --------------------------------------------------------------------------

DIR_BLOQUES = "bloques_a2"


def _tarea(arg):
    """Una (celda, condición): calcula todas las réplicas y persiste un bloque.

    El resultado se escribe en DIR_BLOQUES/<id_celda>_<id_cond>.npy mediante
    escritura atómica (fichero temporal + os.replace), de modo que el barrido
    es resumible: una tarea con bloque ya presente se omite. Esto hace el
    trabajo robusto frente a interrupciones sin alterar ningún resultado
    (cada instancia está determinada por su semilla derivada).
    """
    celda, cond, n_rep = arg
    ruta = os.path.join(DIR_BLOQUES, f"{celda['id_celda']:04d}_{cond['id_cond']:02d}.npy")
    if os.path.exists(ruta):
        return ruta

    pesos = np.array(MEZCLAS_PASO[celda["mezcla"]])
    filas = []
    for rep in range(n_rep):
        res = A.evaluar_instancia(
            n_tray=celda["n_tray"], largo=celda["largo"],
            dwell_medio_ms=celda["dwell_medio_ms"], forma_gamma=celda["forma_gamma"],
            pesos_paso=pesos, modo_memoria=cond["modo"],
            delta=cond["delta"], rho=cond["rho"],
            id_celda=celda["id_celda"], id_cond=cond["id_cond"], replica=rep)
        for r in res:
            filas.append((celda["id_celda"], cond["id_cond"], rep,
                          CANALES.index(r["canal"]), r["k"],
                          ESTIMADORES.index(r["estimador"]),
                          r["im_obs"], r["umbral_p95"], r["nulo_media"],
                          r["p_perm"], float(r["rechaza"])))
    arr = np.array(filas, dtype=np.float64)
    tmp = ruta + f".tmp{os.getpid()}"
    np.save(tmp, arr)
    os.replace(tmp + ".npy", ruta)
    return ruta


COLUMNAS = ["id_celda", "id_cond", "replica", "canal_id", "k", "estimador_id",
            "im_obs", "umbral_p95", "nulo_media", "p_perm", "rechaza"]


def ejecutar(n_procesos: int | None = None, salida: str = "."):
    celdas = construir_celdas_distintas()
    conds = construir_condiciones()

    tareas = []
    for c in celdas:
        for cond in conds:
            n_rep = N_REPLICAS_H0 if cond["id_cond"] == 0 else N_REPLICAS_MEM
            tareas.append((c, cond, n_rep))

    n_procesos = n_procesos or max(1, (os.cpu_count() or 2) - 1)
    print(f"[A2] {len(celdas)} celdas x {len(conds)} condiciones = {len(tareas)} tareas; "
          f"{n_procesos} procesos", flush=True)

    os.makedirs(DIR_BLOQUES, exist_ok=True)
    t0 = time.time()
    rutas = []
    hechas = 0
    # multiprocessing.Pool con contexto 'fork': ProcessPoolExecutor no está
    # disponible en este entorno (SC_SEM_NSEMS_MAX no consultable).
    with mp.get_context("fork").Pool(n_procesos) as ex:
        for ruta in ex.imap_unordered(_tarea, tareas, chunksize=1):
            rutas.append(ruta)
            hechas += 1
            if hechas % 100 == 0 or hechas == len(tareas):
                el = time.time() - t0
                print(f"[A2] {hechas}/{len(tareas)} tareas  {el/60:.1f} min  "
                      f"(ETA {el/hechas*(len(tareas)-hechas)/60:.1f} min)", flush=True)

    arr = np.concatenate([np.load(r) for r in sorted(rutas)], axis=0)
    df = pd.DataFrame(arr, columns=COLUMNAS)
    for c in ("id_celda", "id_cond", "replica", "k", "rechaza"):
        df[c] = df[c].astype(np.int32)
    df["canal"] = pd.Categorical.from_codes(df.pop("canal_id").astype(int), CANALES)
    df["estimador"] = pd.Categorical.from_codes(df.pop("estimador_id").astype(int), ESTIMADORES)
    df = df.merge(pd.DataFrame(celdas), on="id_celda").merge(pd.DataFrame(conds), on="id_cond")
    df = expandir_a_rejilla_completa(df)
    for c in ("im_obs", "umbral_p95", "nulo_media", "p_perm"):
        df[c] = df[c].astype(np.float32)

    ruta_inst = os.path.join(salida, "a2_por_instancia.csv.gz")
    df.to_csv(ruta_inst, index=False, compression="gzip", float_format="%.6g")
    print(f"[A2] por-instancia: {len(df):,} filas -> {ruta_inst}", flush=True)
    return df


# --------------------------------------------------------------------------
# Agregación
# --------------------------------------------------------------------------

CLAVES_CELDA = ["dwell_medio_ms", "forma_gamma", "mezcla", "largo", "n_tray"]


def agregar(df: pd.DataFrame, salida: str = "."):
    g = (df.groupby(["id_celda", *CLAVES_CELDA, "escala_dwell_degenerada",
                     "id_cond", "modo", "delta", "rho",
                     "canal", "k", "estimador"], observed=True)
           .agg(n_replicas=("rechaza", "size"),
                tasa_rechazo=("rechaza", "mean"),
                im_obs_media=("im_obs", "mean"),
                im_obs_sd=("im_obs", "std"),
                nulo_media=("nulo_media", "mean"),
                umbral_p95_medio=("umbral_p95", "mean"))
           .reset_index())
    # error estándar binomial de la tasa de rechazo
    p = g["tasa_rechazo"]
    g["ee_tasa"] = np.sqrt(p * (1 - p) / g["n_replicas"])
    g["metrica"] = np.where(g["modo"] == "H0", "FPR", "potencia")
    g["cumple_criterio"] = np.where(
        g["modo"] == "H0", g["tasa_rechazo"] <= A.FPR_MAX, g["tasa_rechazo"] >= A.POTENCIA_MIN)
    ruta = os.path.join(salida, "a2_por_celda.csv")
    g.to_csv(ruta, index=False, float_format="%.6g")
    print(f"[A2] por-celda: {len(g):,} filas -> {ruta}", flush=True)
    return g


def minima_detectable(g: pd.DataFrame, salida: str = "."):
    """Magnitud mínima (delta o rho) con potencia >= 0.80, por celda x canal x k
    x estimador. Sólo se declara detectable si la celda además PASA el control
    de FPR bajo H0 (criterio conjunto congelado)."""
    fpr = (g[g["modo"] == "H0"]
           .set_index(["id_celda", "canal", "k", "estimador"])["tasa_rechazo"]
           .rename("fpr_h0"))
    m = g[g["modo"] != "H0"].join(fpr, on=["id_celda", "canal", "k", "estimador"])
    m = m[m["fpr_h0"] <= A.FPR_MAX]
    m["magnitud"] = np.where(m["modo"] == "limping", m["delta"], m["rho"])

    filas = []
    llaves = ["id_celda", *CLAVES_CELDA, "canal", "k", "estimador", "modo"]
    for clave, sub in m.groupby(llaves, observed=True):
        ok = sub[sub["tasa_rechazo"] >= A.POTENCIA_MIN]["magnitud"]
        d = dict(zip(llaves, clave))
        d["fpr_h0"] = float(sub["fpr_h0"].iloc[0])
        d["magnitud_minima_detectable"] = float(ok.min()) if len(ok) else np.nan
        d["potencia_max_barrida"] = float(sub["tasa_rechazo"].max())
        filas.append(d)
    out = pd.DataFrame(filas)
    ruta = os.path.join(salida, "a2_minima_detectable.csv")
    out.to_csv(ruta, index=False, float_format="%.6g")
    print(f"[A2] mínima detectable: {len(out):,} filas -> {ruta}", flush=True)
    return out


if __name__ == "__main__":
    n_proc = int(sys.argv[1]) if len(sys.argv) > 1 else None
    df = ejecutar(n_proc)
    g = agregar(df)
    minima_detectable(g)
    print("[A2] listo.", flush=True)
