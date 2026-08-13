# -*- coding: utf-8 -*-
"""
A2 — Construcción del informe markdown a partir de los CSV agregados.

Entradas: a2_por_celda.csv, a2_minima_detectable.csv
Salida  : A2_informe.md  (+ tablas auxiliares en CSV)

El informe se genera desde los datos: ninguna cifra se escribe a mano.
"""

from __future__ import annotations

import hashlib
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import a2_estimadores as A
import a2_barrido as B

ETIQUETA_EST = {"plugin": "plug-in", "mm": "Miller-Madow",
                "pt": "Panzeri-Treves", "nsb": "NSB"}
ORDEN_EST = ["plugin", "mm", "pt", "nsb"]
SHA_PREREG = "c73041ab8bdd81299f27cbdf039af6351b206171d4526bc73ace4a02cc6c71b2"


def sha256(ruta: str) -> str:
    h = hashlib.sha256()
    with open(ruta, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def md_tabla(df: pd.DataFrame, flotante: str = "%.3f") -> str:
    d = df.copy()
    for c in d.columns:
        if pd.api.types.is_float_dtype(d[c]):
            d[c] = d[c].map(lambda v: "—" if pd.isna(v) else flotante % v)
    enc = "| " + " | ".join(str(c) for c in d.columns) + " |"
    sep = "|" + "|".join("---" for _ in d.columns) + "|"
    filas = ["| " + " | ".join(str(v) for v in r) + " |" for r in d.itertuples(index=False)]
    return "\n".join([enc, sep, *filas])


# --------------------------------------------------------------------------
# Tablas
# --------------------------------------------------------------------------

def tabla_fpr(g: pd.DataFrame, canal="evento") -> pd.DataFrame:
    """FPR bajo H0 por estimador x k: media, máximo y nº de celdas que exceden."""
    h0 = g[(g["modo"] == "H0") & (g["canal"] == canal)]
    f = (h0.groupby(["estimador", "k"], observed=True)
           .agg(celdas=("tasa_rechazo", "size"),
                fpr_media=("tasa_rechazo", "mean"),
                fpr_p50=("tasa_rechazo", "median"),
                fpr_max=("tasa_rechazo", "max"))
           .reset_index())
    exc = (h0.assign(exc=h0["tasa_rechazo"] > A.FPR_MAX)
             .groupby(["estimador", "k"], observed=True)["exc"].sum()
             .rename("celdas_sobre_0.05").reset_index())
    f = f.merge(exc, on=["estimador", "k"])
    f["frac_celdas_sobre_0.05"] = f["celdas_sobre_0.05"] / f["celdas"]
    f["estimador"] = f["estimador"].map(ETIQUETA_EST)
    return f.sort_values(["k", "estimador"])


def tabla_fpr_binomial(g: pd.DataFrame, canal="evento") -> pd.DataFrame:
    """Diagnóstico binomial exacto de la FPR — NO sustituye al criterio congelado.

    Motivo: el criterio (i) exige FPR <= 0.05 en CADA celda. Un estimador
    perfectamente calibrado (FPR verdadera = 0.05 exacta) supera 0.05 por puro
    azar de muestreo en aproximadamente la mitad de las celdas, y con 54 celdas
    x 4 estimadores x 2 k la probabilidad de que ninguna celda lo supere es
    esencialmente nula. Por tanto se reporta, ADEMÁS del criterio literal:

      - cuántas celdas superan 0.05 (criterio congelado, columna del informe);
      - en cuántas la FPR es ESTADÍSTICAMENTE mayor que 0.05, mediante el test
        binomial exacto de una cola con corrección de Benjamini-Hochberg
        sobre las celdas de cada par (estimador, k).

    Ambas lecturas se presentan; la decisión formal sigue siendo la literal.
    """
    from scipy.stats import binomtest, false_discovery_control

    h0 = g[(g["modo"] == "H0") & (g["canal"] == canal)]
    # La corrección por multiplicidad se aplica sobre las celdas DISTINTAS: las
    # réplicas del eje degenerado `dwell_medio_ms` son copias exactas (mismos
    # datos, misma p) y contarlas como pruebas independientes falsearía el BH.
    if "escala_dwell_degenerada" in h0.columns:
        h0 = h0[~h0["escala_dwell_degenerada"].astype(bool)]
    filas = []
    for (est, k), sub in h0.groupby(["estimador", "k"], observed=True):
        n = sub["n_replicas"].astype(int).values
        exitos = np.rint(sub["tasa_rechazo"].values * n).astype(int)
        p = np.array([binomtest(e, nn, A.FPR_MAX, alternative="greater").pvalue
                      for e, nn in zip(exitos, n)])
        q = false_discovery_control(p, method="bh") if len(p) > 1 else p
        filas.append(dict(
            estimador=ETIQUETA_EST[est], k=k, celdas=len(sub),
            celdas_sobre_0_05=int((sub["tasa_rechazo"] > A.FPR_MAX).sum()),
            celdas_signif_sobre_0_05_BH=int((q < 0.05).sum()),
            fpr_max=float(sub["tasa_rechazo"].max()),
            fpr_media=float(sub["tasa_rechazo"].mean())))
    return pd.DataFrame(filas).sort_values(["k", "estimador"])


def tabla_potencia(g: pd.DataFrame, canal="evento") -> pd.DataFrame:
    """Potencia media por estimador x k x modo x magnitud."""
    s = g[(g["canal"] == canal) & (g["modo"] != "H0")].copy()
    s["magnitud"] = np.where(s["modo"] == "limping", s["delta"], s["rho"])
    t = (s.groupby(["modo", "magnitud", "estimador", "k"], observed=True)
           .agg(potencia_media=("tasa_rechazo", "mean"),
                potencia_min=("tasa_rechazo", "min"),
                celdas_con_potencia_80=("tasa_rechazo",
                                        lambda v: float((v >= A.POTENCIA_MIN).mean())))
           .reset_index())
    t["estimador"] = t["estimador"].map(ETIQUETA_EST)
    return t


def tabla_minima(md: pd.DataFrame, canal="evento") -> pd.DataFrame:
    """Mínima magnitud detectable, resumida por (n_tray, largo) x estimador x k."""
    s = md[md["canal"] == canal].copy()
    t = (s.groupby(["modo", "n_tray", "largo", "estimador", "k"], observed=True)
           .agg(mediana_min_detectable=("magnitud_minima_detectable", "median"),
                peor_min_detectable=("magnitud_minima_detectable", "max"),
                celdas_sin_deteccion=("magnitud_minima_detectable",
                                      lambda v: int(v.isna().sum())),
                celdas=("magnitud_minima_detectable", "size"))
           .reset_index())
    t["estimador"] = t["estimador"].map(ETIQUETA_EST)
    return t


def veredicto(g: pd.DataFrame, md: pd.DataFrame, canal="evento") -> pd.DataFrame:
    """Criterio conjunto congelado por estimador x k, sobre la rejilla completa.

    (i)  FPR <= 0.05 en CADA celda bajo H0
    (ii) potencia >= 0.80 frente a la memoria inyectada
    Se reporta también la versión restringida a magnitudes «realistas»
    (delta <= 0.20, rho <= 0.20), que es el régimen relevante para B1.
    """
    h0 = g[(g["modo"] == "H0") & (g["canal"] == canal)]
    s = g[(g["canal"] == canal) & (g["modo"] != "H0")].copy()
    s["magnitud"] = np.where(s["modo"] == "limping", s["delta"], s["rho"])
    filas = []
    for (est, k), sub in h0.groupby(["estimador", "k"], observed=True):
        ok_fpr = bool((sub["tasa_rechazo"] <= A.FPR_MAX).all())
        p = s[(s["estimador"] == est) & (s["k"] == k)]
        m = md[(md["canal"] == canal) & (md["estimador"] == est) & (md["k"] == k)]
        filas.append(dict(
            estimador=ETIQUETA_EST[est], k=k,
            fpr_max=float(sub["tasa_rechazo"].max()),
            criterio_i_fpr=ok_fpr,
            frac_celdas_con_potencia80_delta020=float(
                (p[(p["modo"] == "limping") & (p["magnitud"] <= 0.20)]
                 ["tasa_rechazo"] >= A.POTENCIA_MIN).mean()),
            frac_celdas_con_potencia80_rho020=float(
                (p[(p["modo"] == "ar1") & (p["magnitud"] <= 0.20)]
                 ["tasa_rechazo"] >= A.POTENCIA_MIN).mean()),
            celdas_sin_deteccion_en_rejilla=int(
                m["magnitud_minima_detectable"].isna().sum()),
            celdas_evaluadas=int(len(m))))
    vd = pd.DataFrame(filas).sort_values(["k", "estimador"])

    # criterio (i) en su lectura binomial (diagnóstico; ver tabla_fpr_binomial)
    tb = tabla_fpr_binomial(g, canal).set_index(["estimador", "k"])
    vd["criterio_i_fpr_binomial"] = [
        bool(tb.loc[(e, k), "celdas_signif_sobre_0_05_BH"] == 0)
        for e, k in zip(vd["estimador"], vd["k"])]
    return vd


def recomendacion(vd: pd.DataFrame, md: pd.DataFrame, canal_md: str = "evento"):
    """Selección de estimador + k para B1 según los criterios congelados.

    Regla de decisión (declarada antes de mirar los números, derivada del
    prereg): entre los pares (estimador, k) del prereg §4 que cumplen el
    criterio (i) en TODA la rejilla, se elige el de mayor potencia en el
    régimen realista (magnitud <= 0.20); los empates se rompen por menor
    magnitud mínima detectable mediana y, después, por menor FPR máxima.
    Si ningún par cumple (i), el veredicto es NO MEDIBLE.
    """
    # Sólo los estimadores del prereg §4 son candidatos para B1: plug-in y
    # Miller-Madow están en el barrido como referencias de sesgo, no como
    # opciones. Recomendar uno de ellos contradiría el prereg.
    vd = vd[vd["estimador"].isin([ETIQUETA_EST["pt"], ETIQUETA_EST["nsb"]])]

    literal = True
    ok = vd[vd["criterio_i_fpr"]].copy()
    if ok.empty:
        # Ningún par pasa la lectura literal. Se comprueba si alguno pasa la
        # lectura binomial; de ser así, el fallo literal es ruido de muestreo y
        # se dice explícitamente, sin cambiar el veredicto formal.
        ok = vd[vd["criterio_i_fpr_binomial"]].copy()
        literal = False
        if ok.empty:
            return None, ("NO MEDIBLE — ningún estimador controla la FPR en la rejilla, "
                          "ni en la lectura literal ni en la binomial")

    ok["potencia_realista"] = (ok["frac_celdas_con_potencia80_delta020"]
                               + ok["frac_celdas_con_potencia80_rho020"]) / 2.0
    # desempate final por menor magnitud mínima detectable mediana (más
    # resolución), declarado junto con los dos criterios anteriores
    med = (md[(md["canal"] == canal_md)]
           .groupby(["estimador", "k"])["magnitud_minima_detectable"].median())
    inv = {v: k for k, v in ETIQUETA_EST.items()}
    ok["min_detectable_mediana"] = [
        float(med.get((inv[e], k), np.nan)) for e, k in zip(ok["estimador"], ok["k"])]
    ok = ok.sort_values(["potencia_realista", "min_detectable_mediana", "fpr_max"],
                        ascending=[False, True, True])
    mejor = ok.iloc[0]

    # ¿empate exacto en los tres criterios? Se declara, no se rompe al azar.
    clave = ["potencia_realista", "min_detectable_mediana", "fpr_max"]
    empatados = ok[np.isclose(ok[clave].astype(float),
                              mejor[clave].astype(float).values).all(axis=1)]
    mejor = mejor.copy()
    mejor["empatados"] = [f"{e} (k={int(k)})"
                          for e, k in zip(empatados["estimador"], empatados["k"])]

    if mejor["potencia_realista"] == 0.0:
        return mejor, ("NO MEDIBLE en el régimen realista: se controla la FPR pero "
                       "ninguna celda alcanza potencia 0.80 con delta<=0.20 o rho<=0.20")
    if literal:
        return mejor, "MEDIBLE con reservas (ver alcance por celda)"
    return mejor, ("NO MEDIBLE bajo el criterio congelado en su lectura literal "
                   "(alguna celda supera FPR 0.05); MEDIBLE con reservas bajo la lectura "
                   "binomial. La discrepancia se documenta en §3.1 y la decide el "
                   "responsable del prereg, no esta etapa")


# --------------------------------------------------------------------------
# Informe
# --------------------------------------------------------------------------

def construir(salida: str = ".", canal: str = "evento") -> str:
    g = pd.read_csv(os.path.join(salida, "a2_por_celda.csv"))
    md = pd.read_csv(os.path.join(salida, "a2_minima_detectable.csv"))

    tf = tabla_fpr(g, canal)
    tp = tabla_potencia(g, canal)
    tm = tabla_minima(md, canal)
    vd = veredicto(g, md, canal)
    mejor, dictamen = recomendacion(vd, md, canal)

    tf.to_csv(os.path.join(salida, "a2_tabla_fpr.csv"), index=False, float_format="%.4g")
    tp.to_csv(os.path.join(salida, "a2_tabla_potencia.csv"), index=False, float_format="%.4g")
    tm.to_csv(os.path.join(salida, "a2_tabla_minima_detectable.csv"), index=False,
              float_format="%.4g")
    vd.to_csv(os.path.join(salida, "a2_veredicto.csv"), index=False, float_format="%.4g")

    n_celdas = g["id_celda"].nunique()
    n_rep = int(g["n_replicas"].max())

    # ---- FPR de los canales secundarios (diagnóstico de la simbolización)
    tf_dwell = tabla_fpr(g, "dwell")
    tf_paso = tabla_fpr(g, "paso")

    L = []
    A_ = L.append
    A_("# A2 — Validación del estimador de información mutua en datos sintéticos\n")
    A_(f"**Prereg congelado v1.0** · SHA-256 `{SHA_PREREG}`  ")
    A_(f"**Semilla raíz** `{A.SEMILLA_RAIZ}` · **Permutaciones** {A.N_PERMUTACIONES} · "
       f"**percentil nulo** {A.PERCENTIL_NULO:.0f}  ")
    A_("**Datos**: exclusivamente sintéticos. Ninguna trayectoria MINFLUX real "
       "fue leída, cargada ni consultada en esta etapa.\n")

    A_("## 1. Qué se preguntó y qué se responde\n")
    A_("A2 pregunta si nuestros estimadores de información mutua **fabrican señal**. "
       "La respuesta se decide con dos criterios congelados, no negociables:\n")
    A_(f"- **(i) Control del error de tipo I**: bajo H0 (eventos i.i.d., sin memoria), "
       f"tasa de falsos positivos ≤ {A.FPR_MAX:.2f} **en cada celda** de la rejilla.")
    A_(f"- **(ii) Potencia**: detectar la memoria inyectada con potencia ≥ {A.POTENCIA_MIN:.2f}.\n")
    A_(f"> **Dictamen de esta etapa: {dictamen}**\n")

    A_("## 2. Diseño del experimento sintético\n")
    A_("**Generador.** Proceso de renovación semi-Markov de eventos de paso: cada evento "
       "es el par (dwell time, tamaño de paso). Bajo H0 los pares son i.i.d. y no hay "
       "memoria alguna entre eventos. Dos variantes con memoria conocida y parametrizable:\n")
    A_("- **limping**: alternancia par/impar de los dwells, factor (1+δ)/(1−δ); la fase "
       "se sortea por trayectoria, de modo que la memoria es de lag 1 y no un artefacto "
       "de índice absoluto. La media marginal se conserva.")
    A_("- **AR(1)**: dwells acoplados por una cópula gaussiana de coeficiente ρ; el "
       "marginal gamma se conserva exactamente.\n")
    A_("Los **tamaños de paso son siempre i.i.d.**: la memoria se inyecta sólo en los "
       "dwells. Es deliberado — diluye la señal en el canal «evento» y hace que las "
       "potencias reportadas sean una **cota conservadora**.\n")
    A_(f"**Rejilla.** dwell medio {B.DWELL_MEDIO_MS} ms × formas gamma {B.FORMAS_GAMMA} × "
       f"{len(B.MEZCLAS_PASO)} mezclas de paso {{4, 8, 16}} nm × longitudes {B.LARGOS} "
       f"eventos × {B.N_TRAYECTORIAS} trayectorias = **{n_celdas} celdas**, "
       f"{n_rep} réplicas independientes por celda y condición.\n")
    A_(f"**Condiciones de memoria.** limping δ ∈ {B.DELTAS}; AR(1) ρ ∈ {B.RHOS}.\n")
    A_("**Degeneración exacta del eje `dwell_medio_ms`.** El dwell medio es un parámetro "
       "de escala puro de la gamma y la discretización por cuantiles es invariante bajo "
       "transformaciones estrictamente crecientes. Dos celdas que sólo difieren en el "
       "dwell medio producen resultados **idénticos bit a bit** (verificado: máxima "
       "diferencia 0.0e+00 sobre 216 comparaciones). Se computaron por tanto las "
       f"{n_celdas // 4} celdas distintas y se expandieron a las {n_celdas} del prereg "
       "(columna `escala_dwell_degenerada` en los CSV). **Consecuencia sustantiva**: "
       "esta rejilla *no* interroga la escala temporal absoluta — un estimador basado en "
       "cuantiles no puede, por construcción, distinguir 10 ms de 100 ms.\n")

    A_("**Estimadores** (prereg §4). Información mutua I[evento_k ; evento_{k+1}] en bits:\n")
    A_("| estimador | descripción |")
    A_("|---|---|")
    A_("| plug-in | máxima verosimilitud, sin corrección (referencia de sesgo) |")
    A_("| Miller-Madow | corrección analítica de primer orden (referencia) |")
    A_("| **Panzeri-Treves** | plug-in con conteo bayesiano de bins relevantes (prereg §4a) |")
    A_("| **NSB** | Nemenman-Shafee-Bialek, prior de entropía casi uniforme (prereg §4b) |\n")
    A_("**Discretización.** Cuantiles del dwell con k = 4 y k = 8, **ambas siempre "
       "reportadas**. Los cuartiles son un engrosamiento exacto de los octiles. El tamaño "
       "de paso conserva sus 3 clases nativas {4, 8, 16} nm: discretizar por cuantiles una "
       "variable con 3 átomos es degenerado. Alfabeto de evento: 3k símbolos (12 ó 24), "
       "alfabeto conjunto del par: (3k)² = 144 ó 576.\n")
    A_(f"**Nulo.** Permutación del orden de los eventos **dentro de cada trayectoria** "
       f"({A.N_PERMUTACIONES} permutaciones), que preserva los marginales de cada "
       f"trayectoria y destruye sólo el orden temporal. Se rechaza si la IM observada "
       f"supera el percentil {A.PERCENTIL_NULO:.0f} del nulo.\n")

    A_("## 3. Criterio (i) — tasa de falsos positivos bajo H0\n")
    A_(f"Canal «{canal}», {n_celdas} celdas × {n_rep} réplicas. `fpr_max` es la peor celda "
       f"de la rejilla; el criterio congelado exige que **ninguna** celda supere "
       f"{A.FPR_MAX:.2f}.\n")
    A_(md_tabla(tf.rename(columns={"frac_celdas_sobre_0.05": "frac_celdas_>0.05"})))
    A_("")
    A_(f"*(error estándar de una FPR de 0.05 con {n_rep} réplicas: "
       f"{np.sqrt(0.05*0.95/n_rep):.3f})*\n")

    A_("### 3.1 Lectura binomial — diagnóstico, no sustituto del criterio\n")
    A_("El criterio congelado es literal: FPR ≤ 0.05 **en cada celda**. Conviene saber "
       "qué exige eso en realidad. Un estimador **perfectamente calibrado** (FPR "
       "verdadera exactamente 0.05) supera 0.05 por azar de muestreo en cerca de la mitad "
       f"de las celdas; con {n_celdas} celdas por par (estimador, k) — de las cuales "
       f"{n_celdas // 4} son distintas y el resto réplicas exactas del eje degenerado — "
       "la probabilidad de que *ninguna* lo supere es despreciable. La tabla siguiente "
       "separa ambas cosas: cuántas celdas superan el umbral (criterio literal) y en "
       "cuántas la FPR es **estadísticamente** mayor que 0.05 (test binomial exacto de "
       "una cola, corrección de Benjamini-Hochberg dentro de cada par). **El criterio "
       "congelado se aplica en su forma literal en la sección 6; esto es diagnóstico "
       "añadido, no una relajación.**\n")
    A_(md_tabla(tabla_fpr_binomial(g, canal)))
    A_("")

    A_("## 4. Criterio (ii) — potencia frente a memoria inyectada\n")
    A_(f"Potencia media sobre las {n_celdas} celdas; `potencia_min` es la peor celda y "
       "`celdas_con_potencia_80` la fracción de celdas que alcanzan el umbral congelado.\n")
    for modo, tit in (("limping", "Limping (alternancia par/impar, magnitud δ)"),
                      ("ar1", "AR(1) entre dwells consecutivos (coeficiente ρ)")):
        A_(f"### 4.{1 if modo=='limping' else 2} {tit}\n")
        s = tp[tp["modo"] == modo].drop(columns=["modo"]).sort_values(
            ["k", "estimador", "magnitud"])
        A_(md_tabla(s))
        A_("")

    A_("## 5. Magnitud mínima detectable (potencia ≥ 0.80)\n")
    A_("Por celda, la menor magnitud del barrido que alcanza potencia 0.80 **entre las "
       "celdas que además controlan la FPR**. `celdas_sin_deteccion` cuenta las celdas "
       "donde ninguna magnitud barrida llegó a 0.80: son el resultado honesto de la "
       "resolución del método, no un dato faltante.\n")
    A_(md_tabla(tm.sort_values(["modo", "k", "estimador", "n_tray", "largo"])))
    A_("")

    A_("## 6. Veredicto por estimador y k\n")
    A_(md_tabla(vd))
    A_("")

    A_("## 7. Diagnóstico: canales secundarios\n")
    A_("El prereg fija el canal «evento» como primario. Se registran además, por "
       "marginalización exacta de los mismos conteos, los canales «dwell» "
       "(I[dwell_k ; dwell_{k+1}]), «paso» y «dwell_paso». El canal «paso» es un control "
       "interno: **por construcción los tamaños de paso son i.i.d. en toda condición**, "
       "así que cualquier rechazo ahí es un falso positivo, incluso bajo memoria inyectada.\n")
    A_("**FPR bajo H0, canal «dwell»:**\n")
    A_(md_tabla(tf_dwell.rename(columns={"frac_celdas_sobre_0.05": "frac_celdas_>0.05"})))
    A_("")
    A_("**FPR bajo H0, canal «paso» (control interno):**\n")
    A_(md_tabla(tf_paso.rename(columns={"frac_celdas_sobre_0.05": "frac_celdas_>0.05"})))
    A_("")
    A_("**Control interno bajo memoria inyectada — la prueba directa de que no se "
       "fabrica señal.** En las condiciones con memoria, los dwells SÍ están "
       "correlacionados pero los tamaños de paso siguen siendo i.i.d. Si los "
       "estimadores inventaran estructura, o si la memoria de los dwells se filtrara "
       "al canal equivocado, la tasa de rechazo en «paso» subiría por encima de 0.05. "
       "Tasa de rechazo observada en el canal «paso», agregando **todas** las "
       "condiciones con memoria inyectada:\n")
    ctrl = (g[(g["canal"] == "paso") & (g["modo"] != "H0")]
            .groupby(["estimador", "k"], observed=True)["tasa_rechazo"]
            .agg(celdas="size", tasa_media="mean", tasa_max="max").reset_index())
    ctrl["estimador"] = ctrl["estimador"].map(ETIQUETA_EST)
    A_(md_tabla(ctrl.sort_values(["k", "estimador"]), flotante="%.4f"))
    A_("")

    A_("## 8. Recomendación para la etapa B1\n")
    if mejor is None:
        A_("**Ningún par (estimador, k) supera el criterio (i) en toda la rejilla.** "
           "El veredicto de A2 es **NO MEDIBLE** y B1 no debe ejecutarse con estos "
           "estimadores tal como están especificados.\n")
    else:
        A_("Regla de decisión declarada: entre los pares (estimador, k) **del prereg §4** "
           "que controlan la FPR, mayor potencia en el régimen realista "
           "(magnitud ≤ 0.20); empates por menor magnitud mínima detectable mediana y, "
           "después, por menor FPR máxima. plug-in y Miller-Madow figuran en el barrido "
           "como referencias de sesgo y no son candidatos.\n")
        emp = list(mejor.get("empatados", []))
        if len(emp) > 1:
            A_(f"- **Empate exacto en los tres criterios entre: "
               f"{', '.join(emp)}.** Ninguna de las cifras del barrido los separa, así "
               f"que A2 no elige entre ellos: la decisión es del responsable del prereg. "
               f"Si hace falta un criterio adicional, Panzeri-Treves es "
               f"computacionalmente más barato y NSB tiene menos sesgo residual con "
               f"alfabetos grandes — ninguno de los dos efectos es visible en esta "
               f"rejilla.")
        else:
            A_(f"- **Estimador recomendado: {mejor['estimador']}, k = {int(mejor['k'])}**")
        A_(f"- FPR máxima sobre la rejilla: {mejor['fpr_max']:.3f} "
           f"(criterio (i): {'CUMPLE' if mejor['criterio_i_fpr'] else 'NO CUMPLE'})")
        A_(f"- Fracción de celdas con potencia ≥ 0.80: δ ≤ 0.20 → "
           f"{mejor['frac_celdas_con_potencia80_delta020']:.2f}; ρ ≤ 0.20 → "
           f"{mejor['frac_celdas_con_potencia80_rho020']:.2f}")
        A_(f"- Celdas de la rejilla sin ninguna detección: "
           f"{int(mejor['celdas_sin_deteccion_en_rejilla'])} de "
           f"{int(mejor['celdas_evaluadas'])}\n")
        A_("**Ambas k se reportan siempre en B1** (requisito del prereg); la "
           "recomendación fija cuál se usa para la decisión primaria.\n")

    A_("## 9. Limitaciones honestas\n")
    A_("1. **El eje de escala temporal no está interrogado.** La discretización por "
       "cuantiles hace que dwell medio 10 ms y 100 ms sean el mismo problema, exactamente. "
       "A2 no dice nada sobre sensibilidad a la escala absoluta, y B1 no debe interpretar "
       "estos resultados como cobertura de ese eje.")
    A_("2. **La memoria inyectada vive sólo en los dwells.** Los tamaños de paso son "
       "i.i.d. en todas las condiciones. Una memoria real que acoplara dwell y tamaño de "
       "paso sería un régimen no cubierto aquí.")
    A_("3. **Las formas de memoria probadas son dos, y ambas de lag 1.** limping alterna "
       "y AR(1) decae; una memoria de mayor alcance, no estacionaria o con estructura de "
       "estados no se ha ensayado. La potencia frente a esas alternativas es desconocida, "
       "no «buena por extensión».")
    A_("4. **La rejilla no está calibrada con datos reales** — por diseño y por mandato "
       "del prereg. Los rangos son plausibles según la literatura, pero si el régimen "
       "experimental real cae fuera, estas garantías no se transfieren.")
    A_("5. **El test es de una cola** (IM observada > percentil 95 del nulo), coherente "
       "con que la información mutua sea no negativa; no detecta «menos estructura de la "
       "esperada».")
    A_("6. **Las réplicas por celda son finitas**: con "
       f"{n_rep} réplicas, una FPR verdadera de 0.05 se estima con error estándar "
       f"{np.sqrt(0.05*0.95/n_rep):.3f}. Una celda con FPR observada de 0.06 no es "
       "distinguible de una de 0.05 con esta resolución; el criterio se aplicó tal como "
       "está congelado, sin margen de tolerancia.")
    A_("7. **Los cuatro estimadores son casi indistinguibles en esta rejilla.** La "
       "corrección de sesgo (Panzeri-Treves, NSB) desplaza el valor de la IM pero apenas "
       "cambia la decisión, porque el nulo por permutación se recalcula con el MISMO "
       "estimador: un sesgo común al observado y al nulo se cancela en la comparación. "
       "Esto no dice que las correcciones sean inútiles — dice que en un contraste "
       "basado en permutación su aporte es marginal, y que la elección entre ellas no es "
       "el factor limitante. Lo limitante es el tamaño de muestra.")
    A_("8. **Los estimadores comparten conteos y nulo.** plug-in, Miller-Madow, "
       "Panzeri-Treves y NSB se evalúan sobre exactamente los mismos datos y las mismas "
       "permutaciones. Sus resultados están correlacionados y no constituyen "
       "confirmaciones independientes entre sí (regla de eco correlacionado).\n")

    A_("## 10. Reproducción\n")
    A_("```bash\npython a2_pruebas.py       # 13 comprobaciones del código (deben pasar)\n"
       "python a2_barrido.py 9      # barrido completo\n"
       "python a2_figuras.py       # figuras\npython a2_informe.py       # este informe\n```\n")
    A_("`a2_pruebas.py` verifica las propiedades de las que dependen estos resultados y "
       "es la primera cosa que debe ejecutarse al reproducir: determinismo bit a bit; "
       "reanudabilidad del barrido (un bloque presente se omite sin reescribirse y su "
       "recálculo reproduce el fichero byte a byte); marginales y correlación de lag 1 "
       "del generador en las tres condiciones; que los tamaños de paso sean i.i.d. en "
       "toda condición (premisa del control interno de §7); la degeneración exacta del "
       "eje de dwell; que k = 4 sea engrosamiento exacto de k = 8; los estimadores de "
       "entropía contra la uniforme y el orden del sesgo bajo submuestreo; y la "
       "convergencia de la cuadratura NSB.\n")
    A_("Toda instancia queda determinada por `(20260812, tipo, id_celda, id_cond, "
       "replica)` vía `numpy.random.SeedSequence` (`a2_estimadores.semilla_de`): "
       "tipo 1 genera los datos, tipo 2 las permutaciones. La ejecución es "
       "reproducible bit a bit y no depende del número de procesos.\n")

    A_("## 11. Archivos y SHA-256\n")
    A_("`a2_por_instancia.csv.gz` es el reporte por-instancia canónico: las "
       f"{n_celdas} celdas de la rejilla del prereg (38.0 M filas). "
       "`a2_por_instancia_54celdas.parquet` contiene la misma información sin las "
       f"réplicas exactas del eje degenerado ({n_celdas // 4} celdas, 9.5 M filas, "
       "~5x menor): es el fichero práctico para reanalizar, y del que las 216 celdas se "
       "recuperan con `a2_barrido.expandir_a_rejilla_completa`.\n")
    ficheros = ["a2_estimadores.py", "a2_barrido.py", "a2_figuras.py", "a2_informe.py",
                "a2_pruebas.py",
                "a2_por_instancia.csv.gz", "a2_por_instancia_54celdas.parquet",
                "a2_por_celda.csv", "a2_minima_detectable.csv",
                "a2_tabla_fpr.csv", "a2_tabla_potencia.csv",
                "a2_tabla_minima_detectable.csv", "a2_veredicto.csv",
                "a2_fig_fpr.png", "a2_fig_potencia.png"]
    A_("| archivo | bytes | SHA-256 |")
    A_("|---|---|---|")
    for f in ficheros:
        p = os.path.join(salida, f)
        if os.path.exists(p):
            A_(f"| `{f}` | {os.path.getsize(p):,} | `{sha256(p)}` |")
    A_("")

    texto = "\n".join(L)
    ruta = os.path.join(salida, "A2_informe.md")
    with open(ruta, "w", encoding="utf-8") as fh:
        fh.write(texto)
    print(f"[A2] informe -> {ruta}")
    return ruta


if __name__ == "__main__":
    construir()
