# -*- coding: utf-8 -*-
"""
A2 — Validación del estimador en datos sintéticos.
Prereg congelado v1.0, SHA-256 c73041ab8bdd81299f27cbdf039af6351b206171d4526bc73ace4a02cc6c71b2

Módulo determinista. SOLO datos sintéticos: este archivo no lee ni escribe
ninguna trayectoria MINFLUX real.

Contenido:
  1. Generador de un proceso de renovación semi-Markov de eventos de paso
     (pares dwell time / tamaño de paso), con nulo H0 (i.i.d., sin memoria)
     y dos variantes de memoria inyectada conocida:
       (a) limping: alternancia par/impar de dwells de magnitud delta
       (b) AR(1): correlación entre dwells consecutivos de coeficiente rho
  2. Discretización por cuantiles (k = 4 y k = 8) del dwell time; el tamaño
     de paso conserva sus 3 clases nativas (mezcla de 4/8/16 nm).
  3. Estimadores de información mutua I[evento_k ; evento_{k+1}]:
       - plug-in (referencia)
       - plug-in + corrección de sesgo Miller-Madow (referencia)
       - plug-in + corrección de sesgo Panzeri-Treves (prereg §4a)
       - NSB (Nemenman-Shafee-Bialek) (prereg §4b)
  4. Nulo por permutación dentro de trayectoria, 1000 permutaciones,
     percentil 95.

Semilla raíz: 20260812. Derivación determinista declarada en `semilla_de`.

Referencias de método (nombres estándar, Regla 0 — no se renombra nada):
  Panzeri S, Treves A (1996) Network 7:87-107.
  Panzeri S et al. (2007) J Neurophysiol 98:1064-1072.
  Nemenman I, Shafee F, Bialek W (2002) NIPS 14:471-478.
  Ince RAA et al. (2009) Front Neuroinform 3:4 (pyEntropy: conteo bayesiano PT).
  Wolff JO et al. (2023) Science 379:1004-1010 (contexto experimental, NO usado aquí).
"""

from __future__ import annotations

import numpy as np
from scipy.special import gammaln, psi, ndtr
from scipy.stats import gamma as gamma_dist
from scipy.optimize import brentq

# --------------------------------------------------------------------------
# 0. Constantes congeladas
# --------------------------------------------------------------------------

SEMILLA_RAIZ = 20260812
N_PERMUTACIONES = 1000          # prereg §4
PERCENTIL_NULO = 95.0           # prereg §4
FPR_MAX = 0.05                  # criterio congelado (i)
POTENCIA_MIN = 0.80             # criterio congelado (ii)

PASOS_NM = np.array([4.0, 8.0, 16.0])   # mezcla de tamaños de paso
N_CLASES_PASO = 3

LN2 = np.log(2.0)


def semilla_de(*etiquetas: int) -> np.random.Generator:
    """Derivación determinista de generadores a partir de la semilla raíz.

    El estado queda definido por la tupla (SEMILLA_RAIZ, *etiquetas) pasada
    como entropía a np.random.SeedSequence. Es reproducible entre máquinas y
    versiones de numpy (SeedSequence es estable por especificación).
    """
    return np.random.default_rng(np.random.SeedSequence([SEMILLA_RAIZ, *etiquetas]))


# --------------------------------------------------------------------------
# 1. Generador sintético
# --------------------------------------------------------------------------

def generar_instancia(n_tray: int,
                      largo: int,
                      dwell_medio_ms: float,
                      forma_gamma: float,
                      pesos_paso: np.ndarray,
                      modo_memoria: str = "H0",
                      delta: float = 0.0,
                      rho: float = 0.0,
                      rng: np.random.Generator | None = None):
    """Genera una instancia = n_tray trayectorias de `largo` eventos.

    Un evento es el par (dwell time [ms], tamaño de paso [nm]).

    modo_memoria:
      "H0"      : pares i.i.d., sin memoria entre eventos (nulo).
      "limping" : dwells modulados multiplicativamente por paridad,
                  factor (1+delta) / (1-delta); la fase de paridad se sortea
                  por trayectoria, de modo que la memoria es de lag 1 y no
                  un efecto de índice absoluto. La media marginal se conserva.
      "ar1"     : dwells con cópula gaussiana AR(1) de coeficiente rho;
                  el marginal gamma se conserva exactamente.

    Los tamaños de paso son SIEMPRE i.i.d. (multinomial con `pesos_paso`):
    la memoria se inyecta únicamente en la secuencia de dwells. Esto es
    deliberado y se declara en el informe: diluye la señal en el canal
    "evento" y da una cota conservadora de potencia.

    Devuelve (dwells, clases_paso) con forma (n_tray, largo).
    """
    if rng is None:
        rng = semilla_de(0)
    escala = float(dwell_medio_ms) / float(forma_gamma)

    if modo_memoria == "H0":
        dwells = rng.gamma(forma_gamma, escala, size=(n_tray, largo))

    elif modo_memoria == "limping":
        base = rng.gamma(forma_gamma, escala, size=(n_tray, largo))
        fase = rng.integers(0, 2, size=(n_tray, 1))
        paridad = (np.arange(largo)[None, :] + fase) % 2
        factor = np.where(paridad == 0, 1.0 + delta, 1.0 - delta)
        dwells = base * factor

    elif modo_memoria == "ar1":
        z = np.empty((n_tray, largo))
        z[:, 0] = rng.standard_normal(n_tray)
        s = np.sqrt(1.0 - rho ** 2)
        ruido = rng.standard_normal((n_tray, largo - 1))
        for t in range(1, largo):                     # AR(1) estacionario
            z[:, t] = rho * z[:, t - 1] + s * ruido[:, t - 1]
        u = np.clip(ndtr(z), 1e-12, 1.0 - 1e-12)   # cópula gaussiana
        dwells = gamma_dist.ppf(u, forma_gamma, scale=escala)

    else:
        raise ValueError(f"modo_memoria desconocido: {modo_memoria}")

    pesos = np.asarray(pesos_paso, dtype=float)
    pesos = pesos / pesos.sum()
    clases = rng.choice(N_CLASES_PASO, size=(n_tray, largo), p=pesos)
    return dwells, clases.astype(np.int8)


# --------------------------------------------------------------------------
# 2. Discretización por cuantiles
# --------------------------------------------------------------------------

def simbolizar(dwells: np.ndarray, clases_paso: np.ndarray, k: int = 8):
    """Símbolo de evento = bin_dwell * 3 + clase_paso, alfabeto A = 3k.

    Los bordes son los cuantiles empíricos de TODOS los dwells de la instancia
    (agrupados sobre trayectorias). Con k = 8 se usan octiles; los bordes de
    k = 4 (cuartiles) son un subconjunto exacto de los de k = 8, de modo que
    la partición k = 4 es un engrosamiento exacto de la k = 8.

    El tamaño de paso NO se discretiza por cuantiles: es una mezcla de 3
    átomos {4, 8, 16} nm y la discretización por cuantiles sobre una variable
    con 3 átomos es degenerada (bins vacíos o empatados). Se usan sus 3 clases
    nativas, que son la partición más fina posible de esa variable.
    """
    q = np.arange(1, k) / k
    bordes = np.quantile(dwells, q)
    bins = np.searchsorted(bordes, dwells, side="right").astype(np.int8)
    return (bins * N_CLASES_PASO + clases_paso).astype(np.int8)


# --------------------------------------------------------------------------
# 3. Conteos conjuntos de pares consecutivos + nulo por permutación
# --------------------------------------------------------------------------

def conteos_pares(simbolos: np.ndarray,
                  alfabeto: int,
                  n_perm: int,
                  rng: np.random.Generator,
                  tam_bloque: int = 64) -> np.ndarray:
    """Conteos conjuntos de pares (s_k, s_{k+1}) agrupados sobre trayectorias.

    Fila 0 = observado; filas 1..n_perm = permutaciones del orden de los
    eventos DENTRO de cada trayectoria (preserva los marginales de evento de
    cada trayectoria y destruye el orden temporal).

    Devuelve un array (1 + n_perm, alfabeto**2) de enteros.
    """
    n_tray, largo = simbolos.shape
    a2 = alfabeto * alfabeto
    salida = np.empty((1 + n_perm, a2), dtype=np.int64)

    codigos = simbolos[:, :-1].astype(np.int32) * alfabeto + simbolos[:, 1:]
    salida[0] = np.bincount(codigos.ravel(), minlength=a2)

    idx = 1
    while idx <= n_perm:
        p = min(tam_bloque, n_perm - idx + 1)
        rep = np.repeat(simbolos[None, :, :], p, axis=0).reshape(p * n_tray, largo)
        rep = rng.permuted(rep, axis=1)
        cod = rep[:, :-1].astype(np.int32) * alfabeto + rep[:, 1:]
        desplaz = (np.arange(p * n_tray, dtype=np.int32) // n_tray)[:, None] * a2
        plano = (cod + desplaz).ravel()
        salida[idx:idx + p] = np.bincount(plano, minlength=p * a2).reshape(p, a2)
        idx += p
    return salida


def derivar_canales(conteos_evento_k8: np.ndarray, k: int) -> dict:
    """De los conteos conjuntos del alfabeto de evento con k = 8 deriva, por
    marginalización exacta, los conteos de todos los canales y ambas k.

    Canales:
      "evento"     : I[(dwell,paso)_k ; (dwell,paso)_{k+1}]   (primario, prereg)
      "dwell"      : I[dwell_k ; dwell_{k+1}]
      "paso"       : I[paso_k ; paso_{k+1}]
      "dwell_paso" : I[dwell_k ; paso_{k+1}]
    """
    P = conteos_evento_k8.shape[0]
    J = conteos_evento_k8.reshape(P, 8, N_CLASES_PASO, 8, N_CLASES_PASO)
    if k == 4:                                   # engrosamiento exacto de octiles
        J = J.reshape(P, 4, 2, N_CLASES_PASO, 4, 2, N_CLASES_PASO).sum(axis=(2, 5))
    elif k != 8:
        raise ValueError("k debe ser 4 u 8")
    ev = J.reshape(P, k * N_CLASES_PASO, k * N_CLASES_PASO)
    return {
        "evento":     (ev.reshape(P, -1),                       k * N_CLASES_PASO, k * N_CLASES_PASO),
        "dwell":      (J.sum(axis=(2, 4)).reshape(P, -1),       k, k),
        "paso":       (J.sum(axis=(1, 3)).reshape(P, -1),       N_CLASES_PASO, N_CLASES_PASO),
        "dwell_paso": (J.sum(axis=(2, 3)).reshape(P, -1),       k, N_CLASES_PASO),
    }


# --------------------------------------------------------------------------
# 4. Estimadores de información mutua
# --------------------------------------------------------------------------

def _entropia_plugin(conteos: np.ndarray, n: np.ndarray) -> np.ndarray:
    """Entropía plug-in en bits, por filas. `conteos` (P, A), `n` (P,)."""
    with np.errstate(divide="ignore", invalid="ignore"):
        p = conteos / n[:, None]
        t = np.where(conteos > 0, p * np.log2(np.where(p > 0, p, 1.0)), 0.0)
    return -t.sum(axis=1)


def _multiplicidades(conteos: np.ndarray):
    """Representación por multiplicidades de los conteos NO nulos.

    Devuelve (u, M, R_obs): u = valores distintos de conteo no nulo (U,),
    M = (P, U) número de bins con cada valor, R_obs = (P,) bins ocupados.

    Se implementa con un histograma por filas (bincount desplazado) en lugar
    de np.unique: los conteos son enteros acotados por N y el histograma es
    O(P*A) frente al O(P*A log(P*A)) del ordenamiento.
    """
    P = conteos.shape[0]
    mx = int(conteos.max())
    if mx == 0:
        return np.array([1], dtype=np.int64), np.zeros((P, 1)), np.zeros(P)
    desplaz = (np.arange(P, dtype=np.int64) * (mx + 1))[:, None]
    H = np.bincount((conteos + desplaz).ravel(),
                    minlength=P * (mx + 1)).reshape(P, mx + 1)
    H[:, 0] = 0                                   # los bins vacíos no cuentan
    u = np.nonzero(H.any(axis=0))[0]
    M = H[:, u].astype(np.float64)
    return u.astype(np.int64), M, M.sum(axis=1)


# ---- NSB -----------------------------------------------------------------

_CACHE_NSB: dict = {}


def _nodos_nsb(K: int, n_nodos: int = 256):
    """Nodos y pesos de cuadratura para el prior NSB sobre beta.

    NSB fija un prior sobre la concentración beta de la Dirichlet tal que la
    entropía a priori xi(beta) = psi(K*beta+1) - psi(beta+1) sea uniforme en
    (0, log K). Se integra con Gauss-Legendre en xi y se resuelve beta(xi).

    n_nodos = 256 se fijó por estudio de convergencia: para el alfabeto
    conjunto más grande del barrido (K = 576) el error frente a una cuadratura
    de 1200 nodos es < 1e-14 bits, mientras que 48 nodos dejan ~3e-3 bits de
    error de cuadratura — comparable a las diferencias de IM que se están
    midiendo, y por tanto inaceptable aquí.
    """
    if K in _CACHE_NSB:
        return _CACHE_NSB[K]
    x, w = np.polynomial.legendre.leggauss(n_nodos)
    lim = np.log(K)
    xi = 0.5 * lim * (x + 1.0)
    pesos = 0.5 * lim * w

    def xi_de(lb, objetivo):
        b = np.exp(lb)
        return psi(K * b + 1.0) - psi(b + 1.0) - objetivo

    betas = np.empty(n_nodos)
    for i, obj in enumerate(xi):
        lo, hi = -30.0, 30.0
        betas[i] = np.exp(brentq(xi_de, lo, hi, args=(obj,), xtol=1e-12, rtol=1e-14))
    _CACHE_NSB[K] = (betas, pesos)
    return betas, pesos


def _entropia_nsb(conteos: np.ndarray, n: np.ndarray, K: int, mult=None) -> np.ndarray:
    """Entropía NSB en bits, por filas, vectorizada sobre filas y nodos.

    Nemenman, Shafee & Bialek (2002). K es el tamaño NOMINAL del alfabeto
    (incluye bins nunca observados).
    """
    betas, pesos = _nodos_nsb(K)
    u, M, R = _multiplicidades(conteos) if mult is None else mult
    N = n[0]
    Kb = K * betas                                        # (Q,)

    # log-evidencia (salvo constantes independientes de beta)
    g_u = gammaln(u[:, None] + betas[None, :]) - gammaln(betas)[None, :]   # (U,Q)
    logL = gammaln(Kb)[None, :] - gammaln(N + Kb)[None, :] + M @ g_u       # (P,Q)

    # entropía posterior esperada a beta fijo (nats)
    t_u = (u[:, None] + betas[None, :]) * psi(u[:, None] + betas[None, :] + 1.0)
    t_0 = betas * psi(betas + 1.0)                                         # bins vacíos
    suma = M @ t_u + (K - R)[:, None] * t_0[None, :]
    H_beta = psi(N + Kb + 1.0)[None, :] - suma / (N + Kb)[None, :]         # (P,Q)

    lw = logL + np.log(pesos)[None, :]
    lw -= lw.max(axis=1, keepdims=True)
    w = np.exp(lw)
    return (w * H_beta).sum(axis=1) / w.sum(axis=1) / LN2


# ---- Panzeri-Treves ------------------------------------------------------

def _r_bayes_pt(conteos: np.ndarray, n: np.ndarray, dim: int,
                max_iter: int = 600, mult=None) -> np.ndarray:
    """Conteo bayesiano de Panzeri-Treves del número de bins relevantes.

    Procedimiento iterativo de Panzeri & Treves (1996) en la variante
    implementada en pyEntropy (Ince et al. 2009), vectorizado sobre filas
    mediante la representación por multiplicidades.
    """
    P = conteos.shape[0]
    u, M, R = _multiplicidades(conteos) if mult is None else mult
    N = float(n[0])
    R = R.astype(np.float64)
    xtr = np.zeros(P)
    activo = R < dim
    if not activo.any():
        return R

    def esperado(xt, Rn):
        gam = xt * (1.0 - (N / (N + Rn)) ** (1.0 / N))
        pb = ((1.0 - gam) / (N + Rn))[:, None] * (u[None, :] + 1.0)
        ocup = (M * (1.0 - np.exp(N * np.log1p(-np.clip(pb, 0, 1 - 1e-15))))).sum(axis=1)
        pv = np.where(xt > 0, gam / np.maximum(xt, 1e-30), 0.0)
        vac = xt * (1.0 - np.exp(N * np.log1p(-np.clip(pv, 0, 1 - 1e-15))))
        return ocup + vac

    d_prev = np.full(P, float(dim))
    d_act = np.abs(R - esperado(np.zeros(P), R))
    for _ in range(max_iter):
        if not activo.any():
            break
        xtr_p = np.where(activo, xtr + 1.0, xtr)
        e = esperado(xtr_p, R)
        d_new = np.abs(R - e)
        sigue = activo & (d_new < d_act) & (R + xtr_p < dim)
        xtr = np.where(sigue, xtr_p, xtr)
        d_prev = np.where(sigue, d_act, d_prev)
        d_act = np.where(sigue, d_new, d_act)
        activo = sigue
    return np.minimum(R + xtr, float(dim))


def estimar_im(conteos_conj: np.ndarray, Ax: int, Ay: int) -> dict:
    """Información mutua en bits para cada fila de `conteos_conj` (P, Ax*Ay).

    Devuelve dict con las claves 'plugin', 'mm', 'pt', 'nsb'.
    """
    P = conteos_conj.shape[0]
    n = conteos_conj.sum(axis=1).astype(np.float64)
    J = conteos_conj.reshape(P, Ax, Ay)
    cx = J.sum(axis=2)
    cy = J.sum(axis=1)

    Hx = _entropia_plugin(cx, n)
    Hy = _entropia_plugin(cy, n)
    Hxy = _entropia_plugin(conteos_conj, n)
    plugin = Hx + Hy - Hxy

    Rx_n = (cx > 0).sum(axis=1).astype(float)
    Ry_n = (cy > 0).sum(axis=1).astype(float)
    Rxy_n = (conteos_conj > 0).sum(axis=1).astype(float)
    mm = plugin - (Rxy_n - Rx_n - Ry_n + 1.0) / (2.0 * n * LN2)

    mx, my, mxy = (_multiplicidades(cx), _multiplicidades(cy),
                   _multiplicidades(conteos_conj))

    Rx = _r_bayes_pt(cx, n, Ax, mult=mx)
    Ry = _r_bayes_pt(cy, n, Ay, mult=my)
    Rxy = _r_bayes_pt(conteos_conj, n, Ax * Ay, mult=mxy)
    pt = plugin - (Rxy - Rx - Ry + 1.0) / (2.0 * n * LN2)

    nsb = (_entropia_nsb(cx, n, Ax, mult=mx)
           + _entropia_nsb(cy, n, Ay, mult=my)
           - _entropia_nsb(conteos_conj, n, Ax * Ay, mult=mxy))

    return {"plugin": plugin, "mm": mm, "pt": pt, "nsb": nsb}


# --------------------------------------------------------------------------
# 5. Test de permutación por instancia
# --------------------------------------------------------------------------

def test_permutacion(valores: np.ndarray) -> dict:
    """valores[0] = observado, valores[1:] = nulo (1000 permutaciones)."""
    obs = float(valores[0])
    nulo = valores[1:]
    umbral = float(np.percentile(nulo, PERCENTIL_NULO))
    p = (1.0 + np.sum(nulo >= obs)) / (1.0 + nulo.size)
    return {"im_obs": obs, "umbral_p95": umbral,
            "nulo_media": float(nulo.mean()), "nulo_sd": float(nulo.std(ddof=1)),
            "p_perm": float(p), "rechaza": bool(obs > umbral)}


def evaluar_instancia(n_tray, largo, dwell_medio_ms, forma_gamma, pesos_paso,
                      modo_memoria, delta, rho, id_celda, id_cond, replica,
                      n_perm=N_PERMUTACIONES):
    """Pipeline completo de UNA instancia. Devuelve lista de filas por-instancia."""
    rng_dat = semilla_de(1, id_celda, id_cond, replica)
    rng_per = semilla_de(2, id_celda, id_cond, replica)

    dwells, clases = generar_instancia(n_tray, largo, dwell_medio_ms, forma_gamma,
                                       pesos_paso, modo_memoria, delta, rho, rng_dat)
    s8 = simbolizar(dwells, clases, k=8)
    C8 = conteos_pares(s8, 3 * 8, n_perm, rng_per)

    filas = []
    for k in (4, 8):
        for canal, (cc, Ax, Ay) in derivar_canales(C8, k).items():
            est = estimar_im(cc, Ax, Ay)
            for nombre, val in est.items():
                r = test_permutacion(val)
                r.update(canal=canal, k=k, estimador=nombre)
                filas.append(r)
    return filas
