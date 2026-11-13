# Two levers, one budget: a quantitative check of dematerialism (v2.0)

Every number below is generated from data by `make all` and listed, with its source and
uncertainty, in `results/key_numbers.md`. Figures are in `figures/`, captions in
`figures/captions.md`.

## The idea

Classic IPAT with affluence opened up:

```
I = P × ( Aᵉˢˢ + Aᵐᵃᵗ / Aᵈᵉᵐ ) × T          constraint: cumulative I ≤ remaining budget
```

| Term | Meaning |
|---|---|
| Aᵉˢˢ | The sufficiency floor |
| Aᵐᵃᵗ = A − Aᵉˢˢ | Consumption above the floor |
| Aᵈᵉᵐ | Dematerialism: a divisor on Aᵐᵃᵗ only |
| T | Carbon intensity of consumption |

As Aᵈᵉᵐ → ∞, emissions fall to P × Aᵉˢˢ × T, the cost of existing, which only T can reduce.
We estimate whether T alone can meet the budget, how much dematerialism would be needed
otherwise, and whether that much is already observed somewhere.

**Two frames.** The *path frame* (main) lets T fall at a constant annual rate r from its 2026
level and compares cumulative emissions over the horizon H with the remaining budget:
P(Aᵉˢˢ + Aᵐᵃᵗ/Aᵈᵉᵐ)T₀(1 − e^(−kH))/k, with k = −ln(1 − r). This compares a cumulative budget
with a cumulative path. The *static frame* (v1, robustness) compares annual emissions at the
2050 endpoint T with a flat allowance Iᵖˡ = budget ÷ H.

**Counterfactuals.** All scenarios put the whole population at a reference GDP per capita (US,
EU27 or world, 2022) from 2026. They are not forecasts.

## Data and methods in brief

**Data**

| Input | Source |
|---|---|
| Population, GDP (Maddison, 2011 int-$), consumption-based fossil CO₂ for 2022 | OWID, pinned commit |
| Remaining budget from 2026 | Forster et al. 2026, Table 8 (IGCC data repository) |
| Population peak | UN WPP 2024 |
| Global footprint distribution | Chancel 2022 |
| National life evaluation | WHR 2026 data file, 2021–2023 window |

**Floor.** Six published decent-living or distributional routes. Each is converted to a share φ
of the current world mean in its own unit; the floor is the median of the six.

**Frontier.** Quantile regression of life evaluation on log GDP per capita at τ = 0.9, with a
pairs bootstrap.

**Monte Carlo.** One propagation (10,000 draws) over the floor, budget likelihood (33–67 %) and
horizon (25–100 years).

**Tests.** 25 tests pass. They include the P × A × T identity for US, EU27 and World, and the
path closed form checked against numerical quadrature.

## 1. The allowance

The 1.7 °C budget (50 %) is 499 Gt CO₂ from the start of 2026 (391–650 Gt at 67–33 %). Spread
over 50 years, that is an allowance Iᵖˡ of 10 Gt per year. The world emitted 37 Gt of
consumption-based CO₂ in 2022, 3.8 times Iᵖˡ (Fig. 1). At 1.5 °C the budget is 134 Gt, under
four years of current emissions.

## 2. The cost of existing

The six routes put the floor at **Aᵉˢˢ ≈ $6,400 per person per year** (range $4,500–8,700),
39 % of today's world mean. Everyone at the floor at today's T emits **15 Gt per year,
1.5 times Iᵖˡ** (1.0–2.0 across the floor range). At today's carbon intensity, even a decent
life for all does not fit; T must improve (Fig. 1).

The path frame adds a qualification the data impose. If T keeps falling at the historical
rate, the floor alone uses 97 % of the 1.7 °C budget over 50 years. The margin is thin and
fragile: 0.68–1.3 across the floor range, and 1.25 at peak population (10.3 bn in 2084). At
1.5 °C the floor fits only if T falls by at least 9.7 % a year. The statement "T must
improve" therefore holds against today's T. Against the historical trend, it holds for
1.5 °C and for anything above the floor.

## 3. The price of materialism

**At least 69 % of today's CO₂, about 26 Gt per year, is consumption above the floor**
(62–76 % across the floor range).

- **It is a lower bound.** The estimate from group means is a lower bound by Jensen's
  inequality; a Singh–Maddala Lorenz fit to Chancel's groups gives the same 69 %.
- **It is concentrated.** The richest 10 % account for 64 % of the above-floor CO₂.
- **The poorest fall short.** The poorest half lack 3.1 Gt per year of floor consumption
  (Fig. 2).

## 4. Technology alone

With no dematerialism, the budget lasts 4.0 years for a US lifestyle for all at today's T.
To keep a US lifestyle for all within the 1.7 °C budget over 50 years, **T must fall by 22 %
every year, sustained** (Fig. 3).

| Lifestyle | Required decline of T, 1.7 °C | Required decline of T, 2.0 °C |
|---|---|---|
| US | 22 %/yr | 11 %/yr |
| EU27 | 15 %/yr | 7.4 %/yr |
| World average | 6.6 %/yr | 2.3 %/yr |

For comparison:

- The historical scenario rate is 1.5 %. Our consumption-based world T fell 1.76 % per year
  over 2000–2022.
- PwC's 2 °C-consistent rate is 6.9 %.
- No G20 country has exceeded 11.5 % in a single year since 2000.

The US requirement is about 15 times the historical rate and twice the best single year of
any G20 country.

## 5. Dematerialism exists

On the τ = 0.9 frontier, the US life evaluation (6.72) is reached at $21,000 per person, not
$57,000. The frontier slope is 0.71 ladder points per log unit (95 % CI 0.47–0.86); the
published crude cross-country slope is about 0.6. The frontier country uses 29 % of US
above-floor consumption, so **observed Aᵈᵉᵐ = 3.4 for the US (95 % CI 2.4–5.7) and 2.9 for
EU27 (2.0–5.3)** (Fig. 4).

**Named peers.** The free-disposal hull names actual countries: Costa Rica matches the US life
evaluation (Aᵈᵉᵐ 5.4) and Mexico matches EU27's (3.4).

**Sensitivity.**

| Specification | Observed Aᵈᵉᵐ, US |
|---|---|
| τ = 0.9, 2021–2023 ladder (main) | 3.4 |
| τ = 0.95 | 6.3 |
| 2023–2025 ladder | 3.7 |

A quadratic term is not significant (p = 0.53): we find no satiation within the sample.

## 6. Both levers together

Fig. 5a shows the two levers on one map for a US lifestyle at 1.7 °C.

- **Technology alone** needs 22 % a year.
- **Infinite dematerialism** works only once T falls faster than 1.4 % a year, because below
  that the floor itself overshoots.
- **Together, with the dematerialism already observed (3.4), T needs to fall 8.7 % a year**
  (6.2–11 % across the observed Aᵈᵉᵐ interval).

That is 40 % of the rate technology alone would need.

The data do not fully support the step as phrased in the brief, "an ambitious-but-not-extreme
T with an already-observed Aᵈᵉᵐ". The rate needed with observed Aᵈᵉᵐ, 8.7 %, lies between
PwC's Ambitious 6.9 % and Extreme 11.5 %.

**At exactly the Ambitious rate (US, 1.7 °C):**

- Required Aᵈᵉᵐ* is 4.8 (Monte Carlo 90 % interval 3.3–7.5).
- That exceeds the observed point estimate of 3.4, though it lies inside its CI.
- Observed Aᵈᵉᵐ suffices in 8 % of Monte Carlo draws.

**Where Ambitious T and observed Aᵈᵉᵐ do close the gap:**

- US at 2.0 °C: required 1.7, in 100 % of draws.
- EU27 at 1.7 °C: required 3.0 against 2.9 observed, which is the margin; 48 % of draws.

**Where neither lever at any observed value suffices:**

- 1.5 °C: the floor fails at the historical and Ambitious rates. At the Extreme rate the US
  needs Aᵈᵉᵐ ≈ 42.
- The historical rate at 1.7 °C: the US needs 260, and the floor fails in 56 % of draws.

**Static frame (robustness).** The v1 frame is far more lenient: it requires only 2.0 for the
US at Ambitious, 1.7 °C. It counts emissions only at the 2050 endpoint and ignores those on
the way there (S3). Its conclusions on the need for both levers are qualitatively the same,
but it overstates what Ambitious T achieves.

## The result that most challenges the thesis

**Path frame.** The need for dematerialism is a property of rich lifestyles and tight targets,
not of human life as such.

- A world-average lifestyle (about $16,000 per person, close to Costa Rica's GDP per capita)
  fits the 1.7 °C budget with no dematerialism at the Ambitious rate (required 0.94).
- Equal shares at that rate are $17,000 per person, about the Dominican Republic's level.
- Claim 2 of the thesis ("equal shares land at poverty") therefore holds only at today's T
  ($4,300, Ghana) or the historical trend ($6,500, Bolivia).
- At 2.0 °C, the Extreme rate alone just suffices for a US lifestyle (0.96).

**Static frame.** The challenge is stronger: the Ambitious rate alone fits a US lifestyle at
2.0 °C (0.83), as reported in v1.

**What survives in both frames.** At every rate actually sustained so far, T alone fails for
any rich lifestyle at 1.7 °C or below. There, dematerialism of the observed size is what
brings the required decarbonisation down from unprecedented (22 %) to demanding (8.7 %).

## What this study does not show

- **The frontier is cross-sectional, not causal.** It shows that combinations of high life
  evaluation and modest GDP exist. It does not show that a rich country could move to them,
  or by what intervention. Frontier countries differ in climate, institutions, inequality,
  prices and T. The bootstrap covers country sampling, not survey error.
- **GDP is a proxy for consumption, not for wanting.** "Observed Aᵈᵉᵐ" is a ratio of
  consumption levels at equal life evaluation, not a measured capacity to want less. Income,
  not green identity, predicts footprints (Moser & Kleinhückelkotten 2018). Efficiency gains
  rebound (Brockway et al. 2021), and so may money not spent.
- **Climate only.** No other boundary is assessed. Fossil CO₂ only: land-use CO₂ (4.3 Gt in
  2022) counts against the budget but is not in I, so burn rates are understated.
- **Population is held at its 2022 value** except in the rows marked "peak", which hold it at
  the 2084 UN peak.
- **Flat allowances and budget shares are stretched concepts.** A remaining budget is
  cumulative and implies net-zero CO₂ afterwards. Iᵖˡ = budget ÷ 50 years is a device, and
  the path frame does not allocate the budget between countries.
- **Other assumptions.**
  - The floor conversion assumes physical requirements scale with GDP at today's technology.
  - Activity data are for 2022, carried to 2026 at the historical rate.
  - Budget likelihoods cover TCRE uncertainty only.
  - The three T rates are scenario assumptions. They come from a production-based metric
    (PwC), and the PwC press release announced for v2.0 was not available, so they are not
    verified in a primary document.

## What the numbers support

**Support:**

- Today's emissions are about four times the 1.7 °C allowance.
- Even the floor exceeds the allowance at today's T.
- About two-thirds of emissions are consumption above the floor, mostly by the richest tenth.
- T alone would need to fall at rates never sustained for any rich lifestyle.
- The dematerialism already visible across countries would cut the required decarbonisation
  by more than half.

**Do not support:**

- That observed dematerialism combined with PwC's 2 °C rate is enough for a US lifestyle at
  1.7 °C.
- That anything observed is enough for 1.5 °C.
- That the equal-share outcome is poverty once T falls at an ambitious rate.
