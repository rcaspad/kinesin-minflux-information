# Information-theoretic structure of kinesin-1 stepping (MINFLUX, preregistered)

A preregistered, fully deterministic analysis of the public MINFLUX
kinesin-1 dataset of **Wolff et al., *Science* 379, 1004 (2023)**
([Zenodo 7565676](https://doi.org/10.5281/zenodo.7565676), CC-BY-4.0),
asking whether the walk of a molecular motor carries temporal information
beyond its known stepping cycle — and attributing every bit found.

**Author / correspondence:** Raúl Casado Padilla (rcaspad@gmail.com),
independent researcher. Analyses executed with AI assistance (documented
in the audit trail); every result re-executes deterministically from the
code and seeds in this repository.

## Results in one paragraph

Sequential structure is strong in all nine construct × [ATP] conditions
(mutual information between consecutive (dwell, step-size) events:
**0.06–0.51 bits/event**, permutation p at floor), but it is entirely
attributable to the **known 16-nm step / 8-nm substep cycle**, dominated
by labeled-head dynamics. A **calibrated artifact null** (renewal
surrogates with per-trace measured noise, passed through the full
localization + step-detection chain) shows that step detection alone
manufactures adjacent-dwell anticorrelation of the observed order
(lag-1 ≈ −0.08 to −0.11) — killing that naive signature as evidence of
memory. A residual dwell→dwell mutual information that initially survived
in 3 of 7 conditions was then attributed by two further preregistered
tests: it **vanishes in center-of-mass-labeled traces** (N356C; p =
0.30–1.0) and a **cross-motor test on 19 individual microtubules** (2,566
same-site pairs) shows **no shared slow sites** between motors (Spearman
ρ = −0.02, p = 0.77). Conclusion: **no motor memory beyond the stepping
cycle, and no measurable substrate heterogeneity, in this assay** — with
the attribution chain, not by assumption.

## Why the documents are in Spanish

The preregistrations and stage documents were frozen with SHA-256 hashes
**before** data access, and several hashes are recorded in signed
documents; translating them would break hash verification. They are
published verbatim. The code is language-neutral; this README summarizes
everything needed to navigate.

## Repository map (chronological chain)

| dir | stage | key documents |
|---|---|---|
| `prereg/` | Preregistrations | v1 (frozen SHA `c73041ab…`), errata v1.1/v1.2, v2 (post-signature SHA `aa461a6e…`), dataset verification gate |
| `a1/` | Reproducibility gate vs published values | PASS |
| `a2/` | Estimator validation on synthetics (216-cell grid, 500 replicas) | calibration clean; PT/NSB tie |
| `a2_prima/` | Re-validation with the *real* trace-length multisets | PASS + sensitivity map |
| `b1/` | First B1 attempt | INCONCLUSIVE by its own frozen eligibility rule (honest death) |
| `b1v2/` | Confirmatory test on real trajectories (pooled statistic) | CONFIRMS 9/9; survives parity (limping) guard |
| `c/` | Independent Kalafut-Visscher segmentation + **calibrated artifact null** | anticorrelation signature = segmentation-made; residual in 3/7 cells |
| `d/` | Attribution tests (revisits, spatial decay, run lengths, stationary control) | UNRESOLVED by frozen rule |
| `e/` | Tiebreakers: center-of-mass control + cross-motor same-microtubule test | **RESOLVED: labeled-head dynamics** |

## Reproduction

```bash
pip install numpy scipy pandas xlrd
```

Download the dataset ZIP from Zenodo 7565676 (MD5
`473e28f9444427c55540362f84df350e`) and extract it; scripts expect the
authors' directory layout (`KinesinDataFiles/<construct>/<ATP>/`). Every
analysis is seeded from root seed `20260812` via
`numpy.random.SeedSequence` (see `a2/a2_estimadores.py::semilla_de`);
results reproduce **bit-for-bit** across machines (verified across
Python 3.13/3.14, NumPy 2.4.4/2.4.6, SciPy 1.17/1.18 — see
`a2/RECEPCION_A2.md`). Run `a2/a2_pruebas.py` (13 self-checks) before
trusting any reproduction.

## Integrity

- Preregistrations frozen (SHA-256) before opening any trajectory file;
  thresholds never adjusted after data.
- Three public errata issued against our own documents — all caught by
  our own gates.
- Large per-instance result files (549 MB CSV, 114 MB Parquet) are not in
  the repository; their SHA-256 hashes are listed in `a2/A2_informe.md`
  §11 and they are regenerable from seed.

## License

Code: MIT. Documents and figures: CC-BY-4.0. Underlying dataset ©
Wolff & Scheiderer (Zenodo 7565676, CC-BY-4.0) — cite Wolff et al.,
*Science* 379, 1004–1010 (2023), DOI 10.1126/science.ade2650.
