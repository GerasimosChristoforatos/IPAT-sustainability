# IPAT-sustainability

**The idea.** Classic IPAT (impact = population × affluence × technology) treats affluence as
one block. We split it into a sufficiency floor Aᵉˢˢ and consumption above it, Aᵐᵃᵗ, and add
a fourth lever, dematerialism Aᵈᵉᵐ, which divides only the part above the floor:

```
I = P × ( Aᵉˢˢ + Aᵐᵃᵗ / Aᵈᵉᵐ ) × T
```

Using public data, we estimate four things:

- the CO₂ cost of a decent life for all;
- how much of today's CO₂ is consumption above the floor;
- how fast carbon intensity T would have to fall on its own;
- how much dematerialism is already observed across countries, and whether the two levers
  together fit a carbon budget.

The study supports a Perspective article. Results are in `REPORT.md`, quotable numbers in
`results/key_numbers.md`.

## Run

```bash
pip install -r requirements.txt
```

```bash
make all
```

`make all` does the following:

1. Downloads the pinned OWID file.
2. Checks SHA-256 hashes of all raw data.
3. Rebuilds `data/processed/`, `results/`, `figures/` and `results/key_numbers.md`.
4. Runs the tests.

## Repository map

| Path | Content |
|---|---|
| `config.yaml` | every modelling choice: years, scenarios, frames, horizon, Monte Carlo, colour palette |
| `data/sources.yaml` | every external input with value, unit, year, source, location in the source, retrieval date and conversion |
| `data/raw/` | unmodified data files (`SHA256SUMS`); `papers/` is git-ignored and never committed |
| `src/load.py` | raw data to processed tables |
| `src/model.py` | identity and closed forms (static and path frames) |
| `src/floor.py` | floor triangulation (RQ2) and emissions above the floor (RQ3) |
| `src/frontier.py` | wellbeing frontier and observed Aᵈᵉᵐ (RQ6) |
| `src/scenarios.py` | burn rate and required rates (RQ1), cost of existing (RQ4), required Aᵈᵉᵐ* (RQ5) |
| `src/montecarlo.py` | the single uncertainty propagation, both frames |
| `src/figures.py` | five main figures (`figures/`) and three supplementary figures (`figures/si/`) |
| `src/keynumbers.py` | `results/key_numbers.md` |
| `tests/` | identity, closed forms (incl. path frame vs quadrature), estimator checks |
| `results/` | one CSV per research question, each with a `type` column (data / derived / stretched) |
| `figures/captions.md` | captions for all figures |
| `PLAN.md` | the approved Phase 0 plan |

Article PDFs are not redistributed; values taken from them carry section or table references in
`data/sources.yaml`.
