# -*- coding: utf-8 -*-
"""
A2 — Pruebas de verificación del código (ejecutar antes de confiar en un barrido).

    python a2_pruebas.py

Comprueba las propiedades de las que dependen los resultados:
  1. Determinismo bit a bit de una instancia.
  2. Reanudabilidad: un bloque ya presente se omite SIN reescribirse, y su
     recálculo tras borrarlo reproduce el contenido byte a byte.
  3. Marginales del generador: media y CV de la gamma; correlación de lag 1
     nula bajo H0, negativa bajo limping, positiva bajo AR(1).
  4. Los tamaños de paso son i.i.d. en TODA condición (base del control interno).
  5. Degeneración exacta del eje `dwell_medio_ms`.
  6. k = 4 es un engrosamiento exacto de k = 8.
  7. Estimadores de entropía contra un valor conocido (uniforme) y orden
     esperado del sesgo bajo submuestreo severo.
  8. Convergencia de la cuadratura NSB al número de nodos fijado.

Ninguna prueba toca datos reales.
"""

from __future__ import annotations

import hashlib
import os
import shutil
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import a2_estimadores as A
import a2_barrido as B

W = np.array([0.15, 0.70, 0.15])
FALLOS = []


def check(nombre, cond, detalle=""):
    estado = "OK  " if cond else "FALLA"
    print(f"[{estado}] {nombre}" + (f" — {detalle}" if detalle else ""))
    if not cond:
        FALLOS.append(nombre)


def p1_determinismo():
    a = A.evaluar_instancia(100, 100, 25., 2., W, "ar1", 0, .3, 1, 1, 1)
    b = A.evaluar_instancia(100, 100, 25., 2., W, "ar1", 0, .3, 1, 1, 1)
    check("1. determinismo bit a bit de evaluar_instancia", a == b)


def p2_reanudable():
    dir_orig = B.DIR_BLOQUES
    B.DIR_BLOQUES = "._prueba_reanuda"
    shutil.rmtree(B.DIR_BLOQUES, ignore_errors=True)
    os.makedirs(B.DIR_BLOQUES)
    try:
        tarea = (B.construir_celdas_distintas()[0], B.construir_condiciones()[1], 3)
        r1 = B._tarea(tarea)
        h1 = hashlib.sha256(open(r1, "rb").read()).hexdigest()
        mt1 = os.path.getmtime(r1)
        r2 = B._tarea(tarea)                       # debe omitir, no recalcular
        check("2a. bloque existente se omite sin reescribir",
              r1 == r2 and mt1 == os.path.getmtime(r2))
        os.remove(r1)
        h3 = hashlib.sha256(open(B._tarea(tarea), "rb").read()).hexdigest()
        check("2b. recálculo tras borrado es byte a byte idéntico", h1 == h3,
              f"sha256 {h1[:16]}")
    finally:
        shutil.rmtree(B.DIR_BLOQUES, ignore_errors=True)
        B.DIR_BLOQUES = dir_orig


def p3_generador():
    for modo, d, r, signo in (("H0", 0, 0, 0), ("limping", .4, 0, -1), ("ar1", 0, .5, +1)):
        dw, _ = A.generar_instancia(400, 200, 25., 2., W, modo, d, r, A.semilla_de(99))
        lag1 = float(np.mean([np.corrcoef(dw[i, :-1], dw[i, 1:])[0, 1]
                              for i in range(dw.shape[0])]))
        media_ok = abs(dw.mean() - 25.) < 0.5
        if signo == 0:
            sig_ok, esp = abs(lag1) < 0.02, "≈0"
        elif signo < 0:
            sig_ok, esp = lag1 < -0.05, "<0"
        else:
            sig_ok, esp = lag1 > 0.30, ">0"
        check(f"3. generador «{modo}»: media marginal y correlación de lag 1",
              media_ok and sig_ok, f"media={dw.mean():.2f} corr_lag1={lag1:+.4f} (esp. {esp})")


def p4_pasos_iid():
    """Los pasos deben ser i.i.d. en toda condición: la memoria vive en los dwells."""
    peor = 0.0
    for modo, d, r in (("H0", 0, 0), ("limping", .5, 0), ("ar1", 0, .5)):
        _, cl = A.generar_instancia(300, 200, 25., 2., W, modo, d, r, A.semilla_de(7))
        lag1 = float(np.mean([np.corrcoef(cl[i, :-1], cl[i, 1:])[0, 1]
                              for i in range(cl.shape[0])]))
        peor = max(peor, abs(lag1))
    check("4. tamaños de paso i.i.d. en toda condición (base del control interno)",
          peor < 0.02, f"max|corr_lag1| = {peor:.4f}")


def p5_degeneracion():
    peor, n = 0.0, 0
    for fg in (1., 2.):
        for modo, d, r in (("H0", 0, 0), ("limping", .2, 0), ("ar1", 0, .2)):
            ref = None
            for dm in B.DWELL_MEDIO_MS:
                f = A.evaluar_instancia(100, 100, dm, fg, W, modo, d, r, 11, 4, 2)
                v = np.array([x["im_obs"] for x in f] + [x["umbral_p95"] for x in f])
                if ref is None:
                    ref = v
                else:
                    peor = max(peor, float(np.abs(v - ref).max())); n += 1
    check("5. degeneración exacta del eje dwell_medio_ms", peor == 0.0,
          f"{n} comparaciones, max|dif| = {peor:.1e}")


def p6_engrosamiento():
    dw, cl = A.generar_instancia(200, 100, 25., 2., W, "H0", 0, 0, A.semilla_de(3))
    s8 = A.simbolizar(dw, cl, k=8)
    s4 = A.simbolizar(dw, cl, k=4)
    b8, b4 = s8 // A.N_CLASES_PASO, s4 // A.N_CLASES_PASO
    check("6. k=4 es engrosamiento exacto de k=8", np.array_equal(b8 // 2, b4))


def p7_entropias():
    rng = np.random.default_rng(0)
    K = 24
    c = rng.multinomial(50000, np.ones(K) / K)[None, :]
    n = c.sum(axis=1).astype(float)
    hp, hn = A._entropia_plugin(c, n)[0], A._entropia_nsb(c, n, K)[0]
    check("7a. entropía de la uniforme con N grande ≈ log2(K)",
          abs(hp - np.log2(K)) < 5e-3 and abs(hn - np.log2(K)) < 5e-3,
          f"plugin={hp:.4f} nsb={hn:.4f} exacta={np.log2(K):.4f}")
    c = rng.multinomial(30, np.ones(K) / K)[None, :]
    n = c.sum(axis=1).astype(float)
    hp, hn = A._entropia_plugin(c, n)[0], A._entropia_nsb(c, n, K)[0]
    r_pt, r_obs = A._r_bayes_pt(c, n, K)[0], int((c > 0).sum())
    check("7b. bajo submuestreo, plug-in subestima y NSB corrige al alza",
          hp < hn <= np.log2(K) + 1e-9, f"plugin={hp:.3f} < nsb={hn:.3f} ≤ {np.log2(K):.3f}")
    check("7c. Panzeri-Treves estima más bins que los observados",
          r_obs <= r_pt <= K, f"R_obs={r_obs} R_pt={r_pt:.0f} K={K}")


def p8_cuadratura_nsb():
    rng = np.random.default_rng(1)
    K, N = 576, 9900
    c = rng.multinomial(N, rng.dirichlet(np.ones(K) * .3))[None, :]
    n = c.sum(axis=1).astype(float)

    def con(q):
        A._CACHE_NSB.clear()
        orig = A._nodos_nsb
        A._nodos_nsb = lambda K_, n_nodos=q: orig(K_, q)
        try:
            return float(A._entropia_nsb(c, n, K)[0])
        finally:
            A._nodos_nsb = orig
            A._CACHE_NSB.clear()

    err = abs(con(256) - con(1200))
    check("8. cuadratura NSB convergida con 256 nodos (K=576)", err < 1e-9,
          f"|H(256) - H(1200)| = {err:.1e} bits")


if __name__ == "__main__":
    for f in (p1_determinismo, p2_reanudable, p3_generador, p4_pasos_iid,
              p5_degeneracion, p6_engrosamiento, p7_entropias, p8_cuadratura_nsb):
        f()
    print()
    if FALLOS:
        print(f"FALLARON {len(FALLOS)} prueba(s): {', '.join(FALLOS)}")
        sys.exit(1)
    print("Todas las pruebas pasan.")
