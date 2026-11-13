# PLAN.md — Phase 0

Status: **awaiting approval.** Nothing has been downloaded, modelled or computed yet.
Every numeric value quoted below is from the brief or from memory and is **unverified**
until it is recorded in `data/sources.yaml` with a table/page reference.

## 0. What Phase 0 established

| Check | Result |
|---|---|
| Repo | Empty (LICENSE, README stub). Branch `claude/ipat-dematerialism-research-jfscga`. |
| Python libs | pandas, numpy, statsmodels, matplotlib, pyyaml present. **scipy missing** — will `pip install` and pin all versions in `requirements.txt`. |
| Tests | `unittest` (stdlib), because pytest is outside the allowed package list. |
| Network | `raw.githubusercontent.com` reachable. **Blocked (proxy 403):** worldhappiness.report, population.un.org, essd.copernicus.org, zenodo.org, wid.world, doi.org, pwc.co.uk, api.worldbank.org. |
| Upload 1 | Vélez-Henao, Streeck, Kikstra & Pauliuk, *Ecol. Econ.* 240 (2026) 108819 — floor route, usable as is. |
| Upload 2 | WHR 2026 **Ch. 1 (executive summary)**. Narrative only: it does not contain the country ladder table, so it cannot supply RQ6 data. |

Consequence: apart from OWID, all data must arrive via `data/raw/` (list in §7).

## 1. Model and units

`I = P × (Aᵉˢˢ + Aᵐᵃᵗ/Aᵈᵉᵐ) × T`, constraint `I ≤ Iᵖˡ`.

| Quantity | Unit | Type | Definition |
|---|---|---|---|
| `P` | persons | data | OWID (historical), UN WPP 2024 medium variant (peak) |
| `A` | 2011 int-$ / person / yr | data | GDP per capita, Maddison via OWID, 2022 |
| `I` | kg CO₂ / yr (reported Gt) | data | Consumption-based fossil CO₂, Global Carbon Project via OWID, 2022 |
| `T` | kg CO₂ / 2011 int-$ | derived | `I / (P × A)` |
| `I_pl` | Gt CO₂ / yr | stretched | Remaining carbon budget ÷ horizon |
| `A_ess` | 2011 int-$ / person / yr | stretched | Sufficiency floor (RQ2) |
| `A_mat` | 2011 int-$ / person / yr | derived | `A − A_ess` |
| `A_dem` | dimensionless | stretched | Efficiency factor on `A_mat` |

Points a reviewer will raise, stated up front:

1. **The validation test is an identity.** Because `T` is defined as `I/(P·A)`, `P×A×T = I`
   holds by construction. The test therefore verifies data plumbing (EU27 aggregation,
   units, year alignment), not the model. It will be described as such.
2. **Reference lifestyles are counterfactuals:** `I = P_world × A_ref × T_scenario`,
   i.e. everyone at US (EU, world) GDP per capita with *world* carbon intensity.
3. **Fossil CO₂ only.** GCP consumption accounts exclude land-use change (~4 Gt/yr, to be
   verified), while the remaining carbon budget covers all anthropogenic CO₂. `I` is
   therefore understated relative to the budget, which is conservative *against* the thesis.
   Stated as a limitation; no adjustment in the main specification.
4. **Year mismatch.** Activity data 2022; budget counted from its own base year
   (2025 or 2026, see §3). Emissions in between are not deducted. Stated, not patched.
5. **GDP is a proxy for consumption, not for wanting.** Repeated wherever A is interpreted.

## 2. T scenarios (RQ1, RQ4, RQ5)

- `T_2050 = T_2022 × (1 − r)^(2050 − 2022)` for r = 1.5 %, 6.9 %, 11.5 % per year.
  Base and end year in `config.yaml`.
- PwC rates are **production-based energy CO₂ per GDP (market/PPP as PwC defines it)**;
  ours is consumption-based per Maddison GDP. Rates are borrowed as scenario parameters,
  not as measurements of our T. Said in every table caption that uses them.
- **Real-world anchor**, "world at the best current national T": the unrestricted minimum
  is expected to be a low-income country whose low T reflects energy poverty, not
  decarbonisation. Proposed rule (in config): lowest 2022 consumption-based T among
  countries with `A ≥ A_world` and population ≥ 1 million. The unrestricted minimum is
  reported next to it for transparency.

**Proposed addition to RQ1 (one extra column, needs your approval).** The brief's design
is comparative-static: the 2050 endpoint T is applied to a budget counted from today.
That ignores emissions along the path. I propose to report both:

- *Static* (as specified): `years = budget / (P × A × T_2050)`.
- *Path-consistent*: cumulative emissions under exponential decline from 2022,
  `∫ P·A·T₀·e^(−rt) dt`, closed form. Under sustained decline the cumulative total is
  bounded by `I₀/r`, so the budget may never be exhausted.

Reason: the static frame holds A constant, whereas PwC's 6.9 % is the rate required
*given GDP growth*. The path version shows how much of thesis claim 1 depends on the
frame. This is the place the data are most likely to limit the claim, so it should be
visible rather than buried.

## 3. Planetary allowance Iᵖˡ

- Remaining carbon budget at 1.5 / 1.7 / 2.0 °C, 50 %; quantiles 33 % and 67 % for the
  Monte Carlo. Horizon 50 yr; 25 and 100 as sensitivity.
- **Discrepancy to resolve before use.** The brief gives 160 / 530 / 1080 Gt from the
  start of 2025 (IGCC 2024 corrigendum, Table 8). My recollection of the IGCC 2024 paper
  (Forster et al. 2025, ESSD 17, 2641) is 130 / 490 / 1050. I cannot check either from
  here. I will use whatever Table 8 of the latest valid version says, record the version
  and corrigendum status, and if the 2026 edition (ESSD 18, 3889) supersedes it, use that
  with its own base year.
- First attempt: the IGCC public data repository on GitHub (`ClimateIndicator/data`),
  pinned to a commit. Fallback: the paper PDF from you.
- A flat annual allowance is a **stretched concept**: a carbon budget is a cumulative
  quantity and implies net-zero CO₂ after the horizon. Labelled accordingly.

## 4. RQ2 — floor triangulation

**Common unit problem.** Routes come in GJ, tonnes of materials, shares and tCO₂e, while
the model needs `A_ess` in $ per person. Method, applied uniformly:

1. Express each route as a dimensionless ratio `φ = floor / current world mean`, both in
   the route's native unit and year.
2. `A_ess = φ × A_world(2022)`.

Assumption this requires: at today's average technology, the native quantity scales
proportionally with consumption. That is false in detail (decent-living studies assume
needs-targeted, often more efficient provisioning), and is the main reason `A_ess` is
typed *stretched*. Results are reported both in $/person and as `A_ess × T` in
tCO₂/person so the reader can judge.

| # | Route | Native floor | Denominator | Note |
|---|---|---|---|---|
| R1 | Decent living energy, Kikstra et al. 2021 | GJ/cap final energy | World final energy/cap, same year | |
| R2 | DLE, Millward-Hopkins, Hickel & Nag 2025 | GJ/cap | as reported in paper | Same paper as R3 |
| R3 | Decent living materials, same paper | t/cap | t/cap current | Different physical dimension, not independent authorship |
| R4 | Needs-based share, Hickel & Sullivan 2024 | 30 % of current use for 8.5 bn | — | `φ = 0.30 × P_ref / 8.5 bn`, P_ref = population of their reference year |
| R5 | DLS GHG, Vélez-Henao et al. 2026 | tCO₂e/cap | World GHG/cap, same year | **Use the recent-technology (2015) scenario**, not the efficient one. The efficient scenario embeds renewable energy, i.e. a change in T; using it with today's T would double-count. Efficient-scenario value shown separately as a T-adjusted reference. |
| R6 | Global median footprint, Chancel 2022 | tCO₂e/cap | World mean | Descriptive (what the median person uses), not needs-based. Flagged as different in kind. |
| R7 | Monetary poverty line (cross-check) | $/day, World Bank | World mean | Direct monetary route; also the yardstick for thesis claim 2. Needs two WDI series (§7). Optional. |

- Central value: **median of route central values**. Range: min–max of route centrals.
  Each route's own range is shown in the forest plot (Fig. 2).
- No weighting of routes; there is no defensible basis for weights with n ≈ 6.
- Independence is limited (shared authors, shared decent-living framework). Stated.

## 5. RQ3 — price of materialism

- Floor in footprint units: `f = φ_central × mean footprint` (consistent with §4).
- Above-floor emissions: `Σ max(eᵢ − f, 0)`; share applied to `I_world` in CO₂.
  Chancel's figures are CO₂e and include investment-attributed emissions. Stated.
- **Bound.** With group means only, `max(ē_g − f, 0)` understates the true above-floor
  amount, because `max(· − f, 0)` is convex (Jensen). The group-mean estimate is therefore
  a **lower bound**, exact and assumption-free.
- Point estimate: preferably from WID g-percentile data (127 fractiles; optional file,
  §7), which tightens the bound. Fallback: a parametric Lorenz curve (lognormal and
  Singh–Maddala) fitted to the cumulative group shares, fit error reported.
- Also reported: the **deficit** `Σ max(f − eᵢ, 0)`, because `P × A_ess × T` (RQ4) exceeds
  what people below the floor emit today. Fig. 3 shows floor met, floor unmet, above floor.
- Internal check (test): Chancel group means reproduce the stated world mean.

## 6. RQ6 — wellbeing frontier

**Data.** WHR ladder (three-year average) vs `log A` (GDP per capita, 2022, same A as the
model). EU27 ladder = population-weighted mean of members. All countries with both
variables; no exclusions; unweighted (each country is one observed combination).

**Method chosen: linear quantile regression of ladder on log A at τ = 0.9**
(Koenker & Bassett 1978), with pairs bootstrap over countries (2,000 replications, fixed
seed) for confidence intervals.

Why not a DEA / FDH envelope as the primary estimator:
- Deterministic envelopes attribute all deviation to inefficiency. National ladder scores
  carry sampling error (roughly ±0.1 at ~1,000 respondents/yr), so one noisy small country
  would define the frontier.
- Envelopes give no standard inference with scipy/statsmodels; quantile regression does.
- Log-linear form matches the established income–wellbeing literature, so the slope can
  be benchmarked against published values (~0.4 individual, ~0.6 cross-country).

Why FDH is still reported, once, as a secondary result: "already observed somewhere on
Earth" is literally a free-disposal-hull statement (Deprins, Simar & Tulkens 1984). For
each reference, the FDH peer is the lowest-A country with ladder ≥ the reference's ladder.
This names the actual country, and makes no functional-form assumption.

**Observed Aᵈᵉᵐ.** `A_f = exp((W_ref − α̂)/β̂)` from the τ = 0.9 line;
`A_dem_obs = (A_ref − A_ess) / (A_f − A_ess)`. If `A_f ≤ A_ess` it is reported as
"frontier at or below floor", not as a number. FDH version uses the peer's A.

**Checks:** τ = 0.95 as one robustness row; a quadratic term in log A tested for
satiation (reported, not adopted unless significant).

**What it does not show.** Cross-sectional association only. It shows a combination of
wellbeing and consumption that exists, not that a reference country could reach it, nor
by what intervention. Frontier countries differ in climate, institutions, inequality and
T. Counter-evidence cited in the report: Moser & Kleinhückelkotten 2018; Brockway et al.
2021 (rebound).

## 7. Monte Carlo (one propagation)

10,000 draws, fixed seed, for each temperature target, US and EU references.

| Input | Distribution | Basis |
|---|---|---|
| `A_ess` | Triangular(min, central, max) from RQ2 | Range from literature, mode at the median route |
| Budget | Likelihood u ~ U(0.33, 0.67), linear interpolation of published quantiles | Brief |
| Horizon | Triangular(25, 50, 100) yr | Brief |
| T decline rate | See decision D3 | |

- Outcome is on the extended real line: "none needed" is `A_dem* ≤ 1`, "floor fails" is
  `+∞`. Median and 5–95 % interval remain well defined. Probabilities of both categories
  are reported alongside.
- **Caveat stated in the report:** horizon and T rate are choices, not uncertain
  measurements. Mixing them with epistemic uncertainty yields a scenario-weighted
  distribution, not a probability forecast.

## 8. Data I need from you (one batch)

Blocked hosts; I will not work around them. Please place in `data/raw/`. Where I am not
certain of the exact published file name I say so; the name you save under is what I
will record.

**Required**

| # | File | Where |
|---|---|---|
| 1 | WHR 2026 country ladder data, "Data for Figure 2.1" (`.xlsx`; exact name unconfirmed, likely `WHR26_Data_Figure_2.1.xlsx`) | https://www.worldhappiness.report/data-sharing/ |
| 2 | UN WPP 2024, `WPP2024_GEN_F01_DEMOGRAPHIC_INDICATORS_COMPACT.xlsx` (Estimates + Medium variant) | https://population.un.org/wpp/downloads |
| 3 | IGCC carbon-budget table: PDF of the latest Forster et al. paper **and** its corrigendum if one exists (only if the GitHub route in §3 fails; I will tell you) | https://essd.copernicus.org/articles/17/2641/2025/ and the 2026 edition |
| 4 | PwC Net Zero Economy Index 2024, report PDF | https://www.pwc.co.uk/services/sustainability-climate-change/insights/net-zero-economy-index.html |
| 5 | Kikstra et al. 2021, *Environ. Res. Lett.* (PDF) | doi:10.1088/1748-9326/ac1c27 |
| 6 | Hickel & Sullivan 2024, *World Development Perspectives* (PDF) | doi:10.1016/j.wdp.2024.100612 |
| 7 | Millward-Hopkins, Hickel & Nag 2025 (PDF; I do not have a DOI I trust) | — |
| 8 | Chancel 2022, *Nature Sustainability* (PDF) | doi:10.1038/s41893-022-00955-z |
| 9 | IZA DP 18662 (PDF) | https://www.iza.org/publications/dp |

**Optional (each improves one result)**

| # | File | Improves |
|---|---|---|
| 10 | WID world per-capita carbon footprint by g-percentile, 2019 (CSV export) from https://wid.world/data/ | RQ3 point estimate |
| 11 | WDI: `NY.GDP.PCAP.PP.KD` and `NE.CON.PRVT.ZS`, World, 2022 (CSV) from https://data.worldbank.org | Route R7 and the poverty yardstick for thesis claim 2 |

DOIs in rows 5, 6, 8 are from memory; please correct me if they resolve elsewhere.
Paper PDFs will be read for the specific table/page only and **not committed** (copyright);
only extracted values with page references go into `sources.yaml`. Data files are
committed with SHA-256 hashes if their licence allows (to be checked per file).

Fetched by me: OWID CO₂ dataset from GitHub, pinned to a commit hash.

## 9. Decisions I need from you

| # | Decision | My recommendation |
|---|---|---|
| D1 | Add the path-consistent column to RQ1 (§2)? | Yes. One column, closed form, and it is the fairest test of claim 1. |
| D2 | Verifying literature constants: may I use the connected literature tools (SciSpace / Consensus / web fetch of open-access pages) to read published values, or only files you place? | Files you place for anything that enters `config.yaml`; tools only to confirm bibliographic details. Most traceable. |
| D3 | Monte Carlo treatment of T | Run conditional on each T scenario (rate fixed at 1.5 / 6.9 / 11.5 %), plus one pooled row with r ~ U(1.5 %, 11.5 %). A pooled-only result would hide that the T prior drives the answer, and a uniform prior has a median (6.5 %) far above anything sustained historically. |
| D4 | "Best current national T" rule (§2) | Restricted minimum as proposed, unrestricted alongside. |
| D5 | Poverty yardstick for claim 2 | World Bank $2.15, $3.65, $6.85 per day (2017 PPP), converted via row 11. If you prefer no extra data, claim 2 is reported relative to current national GDP per capita levels only. |
| D6 | Alternative to §8 | You could instead allow the blocked hosts in this environment's network settings; I would then fetch the data files myself. Papers behind paywalls would still need to come from you. |

## 10. Phases after approval

1. **Data** — fetch/ingest, `sources.yaml`, `data/processed/`. Commit.
2. **Model + tests** — `model.py`, closed forms, identity test, Chancel consistency. Commit.
3. **RQ1–RQ6 + Monte Carlo** — one CSV per RQ with a `type` column. Commit.
4. **Figures (4) + REPORT.md + README + Makefile.** Commit, 10-line summary.

Stop rule: two failed attempts at anything, then I report what was tried.
