# Figure captions

Values quoted here are in `results/key_numbers.md`. Colours are fixed per concept in every
figure: allowance Iᵖˡ black, floor Aᵉˢˢ blue, materialism Aᵐᵃᵗ orange, technology T green,
dematerialism Aᵈᵉᵐ purple.

## Main figures

**Figure 1 | The gap.** Annual CO₂ in Gt per year.

- *Top bar:* world consumption-based fossil CO₂ in 2022 (Global Carbon Budget via OWID), split
  into consumption at or below the floor and above it. The split uses the 2019 global
  distribution of individual footprints (Chancel 2022) and the central floor.
- *Bottom bar:* the cost of existing, P × Aᵉˢˢ × T. This is everyone at the floor
  (Aᵉˢˢ = $6,400 per person per year, 2011 int-$) at the 2022 world T. The hatched part is
  floor consumption not yet met today.
- *Black line:* the allowance Iᵖˡ, the remaining budget for 1.7 °C at 50 % likelihood from the
  start of 2026 (499 Gt; Forster et al. 2026, Table 8) divided by a 50-year horizon.
- *Caveats:* Iᵖˡ is a flat share of a cumulative budget (a stretched concept). Land-use CO₂
  (4.3 Gt in 2022) is excluded from the bars but counts against the budget. With T falling at
  the historical rate, the floor uses 97 % of the 1.7 °C budget over 50 years (Fig. 3
  framework; `results/rq4_cost_of_existing.csv`).

**Figure 2 | Where materialism sits.** The same CO₂ as Fig. 1 by group of the world population,
ranked by individual footprint.

- *Method:* group means from Chancel (2022, Table 1), with a Singh–Maddala Lorenz curve fitted
  to them (RMSE of cumulative shares 0.003). The floor is φ × Chancel's world mean of 6.0 t
  CO₂e, i.e. 2.3 t CO₂e per person, where φ = 0.39 is the floor's share of the current world
  mean. In CO₂-only accounts the same floor is Aᵉˢˢ × T = 1.8 t CO₂ per person; only the share
  φ is carried across.
- *Above floor:* the sum over individuals of footprint minus floor. With group means alone it
  is a lower bound (Jensen's inequality); the fitted curve gives the same 69 %.
- *Deficit (hatched):* floor consumption the poorest groups do not reach.
- *Caveats:* Chancel's footprints are all-GHG CO₂e for 2019 and include emissions attributed to
  investments. Only the shares are applied to 2022 CO₂.

**Figure 3 | Technology alone.**

- *x:* constant annual rate of decline of T (consumption CO₂ per 2011 int-$) from 2026.
- *y:* years until the 1.7 °C budget (50 %) is used up. The whole 2022 population lives at
  world-average, EU27 or US GDP per capita.
- *Method:* cumulative emissions P·A·T₀·(1 − e^(−kH))/k with k = −ln(1 − r). T₀ is 2022 world
  T carried to 2026 at the historical rate. Dots mark the rate at which the budget lasts the
  full 50-year horizon.
- *Rates:* historical 1.5 %/yr is a scenario value; consumption-based world T in our data fell
  1.76 %/yr over 2000–2022. Ambitious (6.9 %, 2 °C-consistent) and Extreme (11.5 %, faster
  than any G20 country's best single year since 2000) follow PwC's Net Zero Economy Index
  2024. PwC's metric is production-based energy CO₂ per GDP, and these rates were not verified
  in a primary PwC document.
- *Caveat:* counterfactual, not a forecast; GDP held at 2022 levels.

**Figure 4 | Dematerialism already exists.**

- *Data:* national life evaluation (Cantril ladder, three-year average 2021–2023, World
  Happiness Report 2026 data file) against GDP per person in 2022 (Maddison via OWID);
  142 countries.
- *Frontier (purple line):* linear quantile regression on log A at τ = 0.9, solved exactly by
  linear programming. Slope 0.71 ladder points per log unit (pairs bootstrap 95 % CI
  0.47–0.86, 2,000 replications).
- *Arrows:* from each reference to the frontier at the same life evaluation. Observed
  Aᵈᵉᵐ = (A_ref − Aᵉˢˢ)/(A_frontier − Aᵉˢˢ). EU27 is the population-weighted mean of its
  members.
- *Corner note:* sensitivity to τ = 0.95 and to the 2023–2025 window. The free-disposal-hull
  peers are Costa Rica (US, 5.4) and Mexico (EU27, 3.4).
- *Caveats:* cross-sectional association, not a causal effect. GDP per capita is a proxy for
  consumption, not for wanting less.

**Figure 5 | Both levers together.**

- *Panel a:* cumulative CO₂ over 50 years divided by the 1.7 °C budget, for a US lifestyle for
  all, over the decline rate of T (x) and Aᵈᵉᵐ (y, log). White region: within budget. Black
  contour: cumulative CO₂ = budget. Green line: historical rate.
- *Panel a, purple:* the observed US Aᵈᵉᵐ (line) and its bootstrap 95 % CI (band).
- *Panel a, dots:* the rate needed with technology alone (Aᵈᵉᵐ = 1) and with the observed
  Aᵈᵉᵐ.
- *Panel b:* required Aᵈᵉᵐ* = Aᵐᵃᵗ / (B/(P·T₀·F) − Aᵉˢˢ) for each temperature and T scenario,
  US (filled) and EU27 (open), path frame. "Floor fails" means the floor alone exceeds the
  budget, so no Aᵈᵉᵐ suffices. "None needed" means Aᵈᵉᵐ* ≤ 1.
- *Caveats:* Aᵉˢˢ, Aᵈᵉᵐ and the budget allocation are stretched concepts, and population is held
  at 2022. The static frame (S3) gives lower requirements under fast decline.

## Supplementary figures

**Figure S1 | Floor triangulation.** Aᵉˢˢ from six published routes, each expressed as
φ = floor ÷ current world mean in its native unit (final energy, materials, GHG, share of
resource use, median footprint) and multiplied by world GDP per capita in 2022.

- *Bars:* each route's own spread (regional or scenario), not a confidence interval.
- *Band:* min–max of route central values. The central value is their median.
- *R5:* its low end is the efficient scenario, which embeds a change in T.
- *Caveat:* the routes share lineages, so they are not independent.

**Figure S2 | Monte Carlo distribution of required Aᵈᵉᵐ\*.** Empirical CDFs from 10,000 draws
(seed 20261113), path frame, US and EU27, by temperature and T scenario.

- *Draws:* Aᵉˢˢ ~ triangular over the RQ2 range; budget likelihood ~ U(33 %, 67 %), with
  budgets interpolated in the published quantiles; horizon ~ triangular(25, 50, 100) years.
- *Rates:* fixed per scenario, or ~ U(1.5 %, 11.5 %) in the pooled curve.
- *Reading the curves:* a curve ending below 1 means the floor fails in the remaining draws.
- *Caveat:* horizon and rate are choices, not measurements, so this is a scenario-weighted
  distribution, not a probability forecast.

**Figure S3 | Static versus path frame.** Required Aᵈᵉᵐ* (central floor, 50 % budget, 50-year
horizon) in the path frame (filled) and the static frame (open).

- *Static frame:* annual impact at the 2050 endpoint T against budget ÷ horizon, as in v1.
- *Path frame:* cumulative emissions along the declining path against the budget.
- *Why they differ:* under fast decline the static frame ignores emissions before 2050, so it
  requires much less dematerialism.
