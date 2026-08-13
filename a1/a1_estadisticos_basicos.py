#!/usr/bin/env python3
"""A1 — Gate de reproducibilidad (prereg v1.0, SHA c73041ab…cc6c71b2).

Lee las tablas de pasos de los autores (allsteps_reeval.xls) y computa los
estadísticos básicos por constructo × [ATP]. Determinista: sin aleatoriedad.
Semántica de columnas decodificada de functions/process_MF_data.m línea 71:
  [stepx, stepy, tau, sigma_x, sigma_y, phot_x, phot_y, end_flag, transitions]
- end_flag==1 marca fila terminadora de traza (stepx=0): separa trazas.
- tau==0 en filas de padding; dwells válidos: tau>0.
Salida: a1_resultados.csv (por celda) — el desglose por-instancia se genera
en B1; aquí solo el gate de consistencia.
"""
import pandas as pd, numpy as np
from pathlib import Path

BASE = Path(__file__).parent / (
    "extracted/Repository for MINFLUX dissects the unimpeded walking of "
    "kinesin-1/Data repository/KinesinDataFiles")
COLS = ["stepx","stepy","tau","sig_x","sig_y","phot_x","phot_y",
        "end_flag","transitions"]

def main():
    rows = []
    for constructo in ["E215C","K28C","T324C"]:
        for cond in ["10uM","100uM","1mM"]:
            f = BASE / constructo / cond / "allsteps_reeval.xls"
            if not f.exists():
                continue
            df = pd.read_excel(f, header=0)
            df.columns = COLS
            n_traces = int(df["end_flag"].sum())
            steps = df[df["end_flag"] == 0]
            dwells = steps.loc[steps["tau"] > 0, "tau"]
            sx = np.abs(steps["stepx"])
            rows.append({
                "constructo": constructo, "ATP": cond,
                "n_trazas": n_traces, "n_pasos": len(steps),
                "n_dwells_validos": len(dwells),
                "step_mediana_nm": round(float(np.median(sx)), 2),
                "step_media_nm": round(float(np.mean(sx)), 2),
                "dwell_mediana_s": round(float(np.median(dwells)), 4),
                "dwell_media_s": round(float(np.mean(dwells)), 4),
                "sigma_x_mediana_nm": round(
                    float(np.median(steps.loc[steps["sig_x"] > 0, "sig_x"])), 2),
            })
    out = pd.DataFrame(rows)
    out.to_csv(Path(__file__).parent / "a1_resultados.csv", index=False)
    print(out.to_string(index=False))

if __name__ == "__main__":
    main()
