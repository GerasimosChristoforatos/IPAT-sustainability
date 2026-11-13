# IPAT-sustainability

A small, reproducible study that quantifies an extended IPAT identity,

```
I = P × ( Aᵉˢˢ + Aᵐᵃᵗ / Aᵈᵉᵐ ) × T ,   I ≤ Iᵖˡ
```

for CO₂. It supports a Perspective article. Read `REPORT.md` for results and limitations.

## Reproduce

```bash
pip install -r requirements.txt
```

```bash
make all
```

`make all` downloads the pinned OWID file, verifies SHA-256 hashes of all raw data, rebuilds
`data/processed/`, `results/` and `figures/`, and runs the tests.

## Layout

| Path | Content |
|---|---|
| `config.yaml` | every modelling choice and assumption |
| `data/sources.yaml` | every external input: value, unit, year, source, where in the source, retrieval date, conversion |
| `data/raw/` | unmodified data files (`SHA256SUMS`); `papers/` is git-ignored |
| `src/load.py` | raw to processed tables |
| `src/model.py` | identity and closed-form solutions |
| `src/floor.py` | RQ2 floor triangulation, RQ3 emissions above the floor |
| `src/frontier.py` | RQ6 wellbeing frontier and observed Aᵈᵉᵐ |
| `src/scenarios.py` | RQ1 burn rate, RQ4 cost of existing, RQ5 required Aᵈᵉᵐ |
| `src/montecarlo.py` | one uncertainty propagation for required Aᵈᵉᵐ |
| `src/figures.py` | the four figures |
| `tests/` | validation identity, closed forms, estimator checks |
| `results/` | one CSV per research question, each with a `type` column (data / derived / stretched) |
| `PLAN.md` | the approved Phase 0 plan |

## Notes

- Article PDFs are not redistributed. Values taken from them are recorded with section or
  table references in `data/sources.yaml`.
- The UN and World Happiness Report files were supplied by hand because their hosts were
  not reachable from the build environment.
- The three T decline rates are assumptions fixed in `config.yaml`; see the PwC entry in
  `data/sources.yaml` for what was and was not verified.
