# -*- coding: utf-8 -*-
"""
A2 — Figuras diagnósticas del barrido (FPR y curvas de potencia).

Entrada: a2_por_celda.csv (salida de a2_barrido.agregar).
Salida : a2_fig_fpr.png, a2_fig_potencia.png
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import a2_estimadores as A

ESTIMADORES = ["plugin", "mm", "pt", "nsb"]
ETIQUETA_EST = {"plugin": "plug-in", "mm": "Miller-Madow",
                "pt": "Panzeri-Treves", "nsb": "NSB"}
# Los dos estimadores del prereg §4 son los focales; plug-in y Miller-Madow
# son referencias de sesgo y se dibujan con menor peso visual (§4.2).
COLOR_EST = {"pt": "#1f77b4", "nsb": "#d1495b",
             "plugin": "#b8b8b8", "mm": "#8f8f8f"}
ANCHO_EST = {"pt": 1.8, "nsb": 1.8, "plugin": 1.0, "mm": 1.0}


def _estilo():
    """Aplica el estilo de figura si el plugin del skill está cargado."""
    try:
        apply_figure_style(frame="open", sizes=(9, 8, 7))   # noqa: F821
    except NameError:
        plt.rcParams.update({"font.size": 8, "axes.spines.top": False,
                             "axes.spines.right": False, "savefig.dpi": 300,
                             "legend.frameon": False})


def figura_fpr(g: pd.DataFrame, canal: str = "evento", ruta: str = "a2_fig_fpr.png"):
    """Distribución de la FPR sobre las celdas de la rejilla (§6.1).

    Se muestra una nube de puntos por celda en vez de la media, porque el
    criterio congelado se aplica a CADA celda: lo decisivo es la cola superior,
    no el promedio.
    """
    _estilo()
    h0 = g[(g["modo"] == "H0") & (g["canal"] == canal)]
    ntrays = sorted(h0["n_tray"].unique())
    n_celdas = h0["id_celda"].nunique()

    fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.4), sharey=True, squeeze=False)
    rng = np.random.default_rng(A.SEMILLA_RAIZ)

    for col, k in enumerate((4, 8)):
        ax = axes[0, col]
        for j, nt in enumerate(ntrays):
            for i, est in enumerate(ESTIMADORES):
                v = h0[(h0["k"] == k) & (h0["n_tray"] == nt) &
                       (h0["estimador"] == est)]["tasa_rechazo"].values
                if not len(v):
                    continue
                x0 = j + (i - 1.5) * 0.20
                ax.scatter(x0 + rng.uniform(-0.05, 0.05, len(v)), v, s=7,
                           color=COLOR_EST[est], alpha=0.75, lw=0,
                           zorder=3 if est in ("pt", "nsb") else 2)
                ax.plot([x0 - 0.075, x0 + 0.075], [np.median(v)] * 2,
                        color=COLOR_EST[est], lw=1.8, zorder=4,
                        solid_capstyle="butt")
        ax.axhline(A.FPR_MAX, color="k", ls="--", lw=1.0, zorder=1)
        ax.set_xticks(np.arange(len(ntrays)))
        ax.set_xticklabels(ntrays)
        ax.set_xlabel("N trayectorias")
        ax.set_title(f"k = {k}", loc="left")
        ax.margins(x=0.10)
    axes[0, 0].set_ylabel("tasa de falsos positivos")
    axes[0, 1].annotate("umbral congelado 0.05", xy=(1.005, A.FPR_MAX),
                        xycoords=("axes fraction", "data"), ha="left", va="center",
                        fontsize=7)

    # identidad de series por texto directo (§7.3), sin caja de leyenda
    for i, est in enumerate(ESTIMADORES):
        axes[0, 1].annotate(ETIQUETA_EST[est], xy=(1.005, 0.34 - 0.085 * i),
                            xycoords="axes fraction", ha="left", fontsize=7,
                            color=COLOR_EST[est],
                            weight="bold" if est in ("pt", "nsb") else "normal")

    fig.suptitle(f"Tasa de falsos positivos bajo H0 — una marca por celda de la "
                 f"rejilla ({n_celdas} celdas), canal «{canal}»",
                 x=0.01, ha="left", fontsize=9)
    fig.tight_layout(rect=(0, 0, 0.86, 0.92))
    fig.savefig(ruta)
    plt.close(fig)
    return ruta, n_celdas


def figura_potencia(g: pd.DataFrame, canal: str = "evento",
                    ruta: str = "a2_fig_potencia.png"):
    """Potencia frente a la magnitud de memoria inyectada, por N trayectorias."""
    _estilo()
    sub = g[(g["canal"] == canal) & (g["modo"] != "H0")]
    ntrays = sorted(sub["n_tray"].unique())
    n_celdas = sub["id_celda"].nunique()

    fig, axes = plt.subplots(2, len(ntrays), figsize=(2.5 * len(ntrays) + 0.9, 4.6),
                             sharey=True, sharex="row", squeeze=False)
    for fila, (modo, etiqueta) in enumerate((("limping", "limping (δ)"),
                                             ("ar1", "AR(1) (ρ)"))):
        eje = "delta" if modo == "limping" else "rho"
        for col, nt in enumerate(ntrays):
            ax = axes[fila, col]
            s = sub[(sub["modo"] == modo) & (sub["n_tray"] == nt)]
            for est in ESTIMADORES:
                for k, ls in ((4, "-"), (8, "--")):
                    c = s[(s["estimador"] == est) & (s["k"] == k)]
                    if not len(c):
                        continue
                    m = c.groupby(eje)["tasa_rechazo"].mean()
                    ax.plot(m.index, m.values, ls, color=COLOR_EST[est],
                            lw=ANCHO_EST[est], marker="o" if k == 4 else "s",
                            ms=2.8, zorder=3 if est in ("pt", "nsb") else 2)
            ax.axhline(A.POTENCIA_MIN, color="k", ls=":", lw=1.0, zorder=1)
            ax.set_ylim(-0.05, 1.05)
            ax.margins(x=0.06)
            if fila == 0:
                ax.set_title(f"N = {nt} trayectorias", loc="left")
            if fila == 1:
                ax.set_xlabel("magnitud de memoria inyectada")
            if col == 0:
                ax.set_ylabel(f"potencia — {etiqueta}")
    axes[0, 0].annotate("potencia 0.80 exigida", xy=(0.03, A.POTENCIA_MIN),
                        xycoords=("axes fraction", "data"),
                        xytext=(0, -10), textcoords="offset points", fontsize=7)
    for i, est in enumerate(("pt", "nsb", "plugin", "mm")):
        axes[1, len(ntrays) - 1].annotate(
            ETIQUETA_EST[est], xy=(0.97, 0.34 - 0.09 * i), xycoords="axes fraction",
            ha="right", fontsize=7, color=COLOR_EST[est],
            weight="bold" if est in ("pt", "nsb") else "normal")
    axes[1, len(ntrays) - 1].annotate(
        "línea continua: k = 4\nlínea discontinua: k = 8", xy=(0.97, 0.62),
        xycoords="axes fraction", ha="right", va="top", fontsize=7)
    fig.suptitle(f"Potencia frente a memoria inyectada — media sobre "
                 f"{n_celdas} celdas de la rejilla, canal «{canal}»",
                 x=0.01, ha="left", fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(ruta)
    plt.close(fig)
    return ruta, n_celdas


if __name__ == "__main__":
    g = pd.read_csv("a2_por_celda.csv")
    print(figura_fpr(g))
    print(figura_potencia(g))
