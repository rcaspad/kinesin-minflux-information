# -*- coding: utf-8 -*-
"""
Etapa C — Cargador de crudo MINFLUX + segmentador Kalafut-Visscher (BIC).

Congelación: C_CONGELACION.md. Independiente del pipeline de los autores:
posiciones FPGA crudas (no SCE), filtro mínimo declarado, KV sin parámetros.
"""
from __future__ import annotations

import glob
import os

import numpy as np

MIN_FOT, MAX_FOT = 7.0, 150.0     # constantes de los autores, sin suavizado


def cargar_fichero_crudo(ruta: str):
    """Trazas de un .txt MINFLUX: lista de dicts con t, x, y (filtradas)."""
    with open(ruta, "r", encoding="latin-1") as f:
        cabecera = f.readline()
    partes = cabecera.replace("data size", "").split()
    n_tray, n_loc = int(partes[0]), int(partes[1])
    datos = np.loadtxt(ruta, skiprows=1)
    if datos.shape[0] != n_tray * n_loc:
        raise ValueError(f"{ruta}: {datos.shape[0]} filas != {n_tray}x{n_loc}")
    trazas = []
    for k in range(n_tray):
        b = datos[k * n_loc:(k + 1) * n_loc]
        x, y = b[:, 0], b[:, 1]
        nx = b[:, 4:7].sum(axis=1)
        lx = b[:, 13]
        dt = b[:, 17]
        t = np.concatenate([[0.0], np.cumsum(dt[:-1])])
        vld = (lx > 0) & (nx >= MIN_FOT) & (nx <= MAX_FOT)
        if vld.sum() < 20:
            continue
        trazas.append(dict(t=t[vld], x=x[vld], y=y[vld]))
    return trazas


def proyectar_eje_principal(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Rotación al eje de avance (misma geometría que los autores)."""
    xc, yc = x - x.mean(), y - y.mean()
    pend = (xc * yc).sum() / (xc * xc).sum()
    a = -np.arctan(pend)
    return np.cos(a) * x - np.sin(a) * y


def kv_bic(z: np.ndarray, max_pasos: int = 200) -> np.ndarray:
    """Kalafut-Visscher: puntos de cambio por BIC gaussiano, sin parámetros.

    Devuelve los índices de cambio ordenados (frontera = primer índice del
    tramo nuevo). Voraz: añade el cambio que minimiza la SS residual y lo
    acepta solo si el BIC baja.
    """
    n = len(z)
    c1 = np.concatenate([[0.0], np.cumsum(z)])
    c2 = np.concatenate([[0.0], np.cumsum(z * z)])

    def ss(i, j):                     # suma de cuadrados residual de z[i:j]
        s, s2, m = c1[j] - c1[i], c2[j] - c2[i], j - i
        return s2 - s * s / m

    def mejor_corte(i, j):
        if j - i < 2:
            return None, np.inf
        m = np.arange(i + 1, j)
        s_izq = c1[m] - c1[i]
        q_izq = c2[m] - c2[i]
        s_der = c1[j] - c1[m]
        q_der = c2[j] - c2[m]
        ss_tot = (q_izq - s_izq ** 2 / (m - i)) + (q_der - s_der ** 2 / (j - m))
        b = int(np.argmin(ss_tot))
        return int(m[b]), float(ss_tot[b])

    cortes = []          # índices de cambio aceptados
    ss_actual = ss(0, n)
    bic_actual = n * np.log(max(ss_actual, 1e-12) / n) + 2 * np.log(n)
    for _ in range(max_pasos):
        bordes = [0] + sorted(cortes) + [n]
        cand, mejor = None, np.inf
        for i, j in zip(bordes[:-1], bordes[1:]):
            corte, ss_nueva = mejor_corte(i, j)
            if corte is None:
                continue
            delta = ss_nueva - ss(i, j)
            if delta < mejor:
                mejor, cand = delta, corte
        if cand is None:
            break
        ss_prop = ss_actual + mejor
        k = len(cortes) + 1
        bic_prop = n * np.log(max(ss_prop, 1e-12) / n) + (k + 2) * np.log(n)
        if bic_prop >= bic_actual:
            break
        cortes.append(cand)
        ss_actual, bic_actual = ss_prop, bic_prop
    return np.array(sorted(cortes), dtype=int)


def fusionar_tramos(bordes: np.ndarray, niveles: np.ndarray,
                    min_paso: float = 5.0):
    """Fusión iterativa de tramos con salto < min_paso (C_CONGELACION_v2).

    Funde en cada iteración el par adyacente con el salto más pequeño
    (nivel = media ponderada por longitud) hasta que todos los saltos
    sean >= min_paso.
    """
    b = list(bordes)
    nv = list(niveles)
    while len(nv) > 1:
        saltos = np.abs(np.diff(nv))
        i = int(np.argmin(saltos))
        if saltos[i] >= min_paso:
            break
        w1, w2 = b[i + 1] - b[i], b[i + 2] - b[i + 1]
        nv[i] = (nv[i] * w1 + nv[i + 1] * w2) / (w1 + w2)
        del nv[i + 1]
        del b[i + 1]
    return np.array(b), np.array(nv)


def eventos_de_traza(t: np.ndarray, z: np.ndarray):
    """Secuencia de eventos (dwell interior, |paso| al final del tramo)."""
    cortes = kv_bic(z)
    if len(cortes) < 3:
        return []
    bordes = np.concatenate([[0], cortes, [len(z)]])
    niveles = np.array([z[i:j].mean() for i, j in zip(bordes[:-1], bordes[1:])])
    bordes, niveles = fusionar_tramos(bordes, niveles)
    if len(niveles) < 4:
        return []
    eventos = []
    # tramos interiores: 1 .. len(niveles)-2 (primero y último censurados);
    # paso CON SIGNO (el nulo de artefacto lo necesita); clases usan |paso|
    for s in range(1, len(niveles) - 1):
        dwell = t[bordes[s + 1] - 1] - t[bordes[s]]
        paso = niveles[s + 1] - niveles[s]
        if dwell > 0:
            eventos.append((float(dwell), float(paso)))
    return eventos


def segmentar_celda(dir_celda: str):
    """Todas las trazas KV-segmentadas de una celda, en orden determinista."""
    salida = []
    for f in sorted(glob.glob(os.path.join(dir_celda, "*.txt"))):
        for tr in cargar_fichero_crudo(f):
            z = proyectar_eje_principal(tr["x"], tr["y"])
            ev = eventos_de_traza(tr["t"], z)
            salida.append(dict(t=tr["t"], z=z, eventos=ev))
    return salida
