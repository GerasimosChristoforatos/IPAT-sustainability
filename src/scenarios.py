"""RQ1 (burn rate, technology alone), RQ4 (cost of existing, equal shares) and
RQ5 (required dematerialism), in two frames.

path (main)        cumulative emissions over the horizon H, with T falling at a constant
                   annual rate r from T0 (its level at the budget base year), against the
                   remaining budget. Physically consistent: a cumulative budget against a
                   cumulative path.
static (robustness) annual impact at the 2050 endpoint T against a flat allowance
                   I_pl = budget / H, as in v1.

All scenarios are counterfactuals: the whole population (2022, or the UN peak) lives at
the reference A (US, EU27 or World GDP per capita of 2022) from the budget base year.
They are not forecasts. T rates are scenario assumptions (see sources.yaml, pwc_nzei_2024).
"""
import math
from pathlib import Path

import numpy as np
import pandas as pd

import load
import model

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
A_COL, T_COL = "A_int2011usd_per_cap", "T_kgCO2_per_int2011usd"
LEVELS = ("central", "low", "high")


# ----------------------------------------------------------------------------
# Inputs shared by both frames
# ----------------------------------------------------------------------------
def rates(cfg: dict) -> dict:
    """Scenario name -> annual decline rate of T. 'Today' = no further decline."""
    return {"Today": 0.0, **{k: v["rate"] for k, v in cfg["t_scenarios"].items()}}


def t0_path(cfg: dict, units: pd.DataFrame) -> float:
    """T at the budget base year [kg CO2 / 2011 int-$]: the activity-year world T, carried
    forward over the gap years at the BAU (historical) rate (config: frames.path_t0_rule)."""
    yrs = cfg["years"]
    T = units.set_index("unit").loc["World", T_COL]
    return model.t_endpoint(T, cfg["t_scenarios"]["BAU"]["rate"], yrs["budget_base"] - yrs["activity"])


def t_scenarios(cfg: dict, units: pd.DataFrame, countries: pd.DataFrame) -> pd.DataFrame:
    """Static-frame T levels: today, three 2050 endpoints, and national anchors."""
    w = units.set_index("unit").loc["World"]
    T_now, span = w[T_COL], cfg["years"]["t_scenario_end"] - cfg["years"]["activity"]
    rows = [dict(scenario="Today", kind="observed", rate=0.0, T=T_now,
                 note=f"world T, {cfg['years']['activity']}", type="derived")]
    for name, spec in cfg["t_scenarios"].items():
        rows.append(dict(scenario=name, kind="rate endpoint", rate=spec["rate"],
                         T=model.t_endpoint(T_now, spec["rate"], span),
                         note=f"{spec['label']}: {spec['rate']:.1%}/yr to {cfg['years']['t_scenario_end']}",
                         type="derived (rate is an assumption)"))
    rule = cfg["best_national_t"]
    elig = countries[(countries.P_persons >= rule["min_population"]) &
                     (countries[A_COL] >= rule["min_a_relative_to_world"] * w[A_COL])].sort_values(T_COL)
    for rank, label in ((0, "Best national"), (1, "Second-best national")):
        c = elig.iloc[rank]
        rows.append(dict(scenario=label, kind="national anchor", rate=np.nan, T=c[T_COL],
                         note=f"{c.country} (A >= world mean, population >= {rule['min_population']:,})",
                         type="derived"))
    if rule["report_unrestricted"]:
        c = countries.sort_values(T_COL).iloc[0]
        rows.append(dict(scenario="Lowest national, unrestricted", kind="national anchor", rate=np.nan,
                         T=c[T_COL], note=c.country, type="derived"))
    return pd.DataFrame(rows)


def budgets(cfg: dict, budget: pd.DataFrame) -> dict:
    """Remaining budget [kg CO2] by target at the central likelihood."""
    b = budget[np.isclose(budget.likelihood, cfg["carbon_budget"]["likelihood_central"])]
    return {float(r.target_degC): r.budget_GtCO2 * model.KG_PER_GT for _, r in b.iterrows()}


def _refs(cfg, units):
    u = units.set_index("unit")
    return [(ref, u.loc[ref, A_COL]) for ref in cfg["units"]]


# ----------------------------------------------------------------------------
# RQ1
# ----------------------------------------------------------------------------
def rq1(cfg, units, tsc, pop, bud) -> pd.DataFrame:
    """Burn rate. Path frame: cumulative emissions over H / budget and years to exhaustion
    (unbounded horizon). Static frame: annual impact at the endpoint T / I_pl."""
    H, sens = cfg["carbon_budget"]["horizon_years"], cfg["carbon_budget"]["horizon_sensitivity_years"]
    T0 = t0_path(cfg, units)
    rows = []
    for ref, A in _refs(cfg, units):
        for _, p in pop.iterrows():
            for target, B in bud.items():
                base = dict(reference=ref, A_ref=A, P_label=p.label, P_persons=p.P_persons,
                            target_degC=target, budget_GtCO2=B / model.KG_PER_GT)
                I0 = p.P_persons * A * T0
                for name, r in rates(cfg).items():
                    y = model.years_to_exhaustion_path(B, I0, r) if r > 0 else B / I0
                    row = dict(base, frame="path", T_scenario=name, rate=r, T_start=T0,
                               I_start_GtCO2_per_yr=I0 / model.KG_PER_GT,
                               cumulative_over_budget=I0 * model.path_factor(r, H) / B,
                               years_to_exhaustion=y,
                               exhaustion_year=(math.inf if math.isinf(y) else cfg["years"]["budget_base"] + y))
                    for h in sens:
                        row[f"cumulative_over_budget_h{h}"] = I0 * model.path_factor(r, h) / B
                    rows.append(row)
                for _, t in tsc.iterrows():
                    I = p.P_persons * A * t["T"]
                    row = dict(base, frame="static", T_scenario=t.scenario, rate=t.rate, T_start=t["T"],
                               I_start_GtCO2_per_yr=I / model.KG_PER_GT,
                               cumulative_over_budget=I / model.allowance(B, H),   # = I / I_pl
                               years_to_exhaustion=B / I, exhaustion_year=cfg["years"]["budget_base"] + B / I)
                    for h in sens:
                        row[f"cumulative_over_budget_h{h}"] = I / model.allowance(B, h)
                    rows.append(row)
    out = pd.DataFrame(rows)
    out["horizon_yr"] = H
    out["note"] = np.where(out.frame == "path", "cumulative over H / budget",
                           "annual impact at endpoint T / I_pl (= budget / H)")
    out["type"] = "derived (counterfactual; I_pl is stretched)"
    return out


def required_rates(cfg, units, pop, bud, floor, observed: pd.DataFrame) -> pd.DataFrame:
    """Technology alone and both levers together (path frame): the smallest constant annual
    decline rate of T that keeps cumulative emissions over H within the budget, for
    A_dem = 1 (technology alone), A_dem = observed (central and 95 % CI), and the floor only."""
    H = cfg["carbon_budget"]["horizon_years"]
    T0 = t0_path(cfg, units)
    qr = observed[observed.method.str.startswith("quantile") & (observed.A_ess_level == "central")].set_index("reference")
    rows = []
    for ref, A in _refs(cfg, units):
        cases = [("technology alone (A_dem = 1)", 1.0)]
        if ref in qr.index:
            o = qr.loc[ref]
            cases += [("both levers (observed A_dem, central)", o.A_dem_observed),
                      ("both levers (observed A_dem, CI low)", o.A_dem_observed_ci_low),
                      ("both levers (observed A_dem, CI high)", o.A_dem_observed_ci_high)]
        cases += [("floor only (A_dem -> inf)", math.inf)]
        for _, p in pop.iterrows():
            for target, B in bud.items():
                for case, d in cases:
                    I0 = model.impact(p.P_persons, floor["central"], A - floor["central"], d, T0)
                    rows.append(dict(frame="path", reference=ref, P_label=p.label, target_degC=target,
                                     horizon_yr=H, case=case, A_dem=d, I_start_GtCO2_per_yr=I0 / model.KG_PER_GT,
                                     required_rate=model.required_rate_path(I0, B, H),
                                     required_rate_unbounded_horizon=min(1.0, 1.0 - math.exp(-I0 / B)),
                                     type="derived (counterfactual)"))
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------
# RQ4
# ----------------------------------------------------------------------------
def _nearest_country(A_eq: float, gdp_countries: pd.DataFrame) -> str:
    """Country with the closest GDP per capita (log distance); flagged if out of range."""
    if A_eq > gdp_countries[A_COL].max():
        return "above every country in the sample"
    if A_eq < gdp_countries[A_COL].min():
        return "below every country in the sample"
    i = (np.log(gdp_countries[A_COL]) - math.log(A_eq)).abs().idxmin()
    return gdp_countries.loc[i, "country"]


def rq4(cfg, units, tsc, pop, bud, floor, gdp_countries) -> pd.DataFrame:
    """Cost of existing (everyone at A_ess) against the budget, and the equal-share A."""
    H = cfg["carbon_budget"]["horizon_years"]
    T0 = t0_path(cfg, units)
    rows = []
    for _, p in pop.iterrows():
        for target, B in bud.items():
            base = dict(P_label=p.label, P_persons=p.P_persons, target_degC=target, horizon_yr=H,
                        budget_GtCO2=B / model.KG_PER_GT, I_pl_GtCO2_per_yr=model.allowance(B, H) / model.KG_PER_GT)
            for name, r in rates(cfg).items():
                F = model.path_factor(r, H)
                row = dict(base, frame="path", T_scenario=name, rate=r, T_start=T0)
                for lv in LEVELS:
                    c = model.cost_of_existing(p.P_persons, floor[lv], T0)
                    row[f"cost_of_existing_GtCO2_per_yr_{lv}"] = c / model.KG_PER_GT
                    row[f"cost_over_budget_{lv}"] = c * F / B
                    row[f"floor_fits_{lv}"] = bool(c * F <= B)
                A_eq = model.equal_share_A_path(B, p.P_persons, T0, r, H)
                row.update(equal_share_A_int2011usd_per_cap=A_eq, equal_share_A_over_A_ess_central=A_eq / floor["central"],
                           equal_share_lives_like=_nearest_country(A_eq, gdp_countries))
                rows.append(row)
            for _, t in tsc.iterrows():
                I_pl = model.allowance(B, H)
                row = dict(base, frame="static", T_scenario=t.scenario, rate=t.rate, T_start=t["T"])
                for lv in LEVELS:
                    c = model.cost_of_existing(p.P_persons, floor[lv], t["T"])
                    row[f"cost_of_existing_GtCO2_per_yr_{lv}"] = c / model.KG_PER_GT
                    row[f"cost_over_budget_{lv}"] = c / I_pl
                    row[f"floor_fits_{lv}"] = bool(c <= I_pl)
                A_eq = model.equal_share_A(I_pl, p.P_persons, t["T"])
                row.update(equal_share_A_int2011usd_per_cap=A_eq, equal_share_A_over_A_ess_central=A_eq / floor["central"],
                           equal_share_lives_like=_nearest_country(A_eq, gdp_countries))
                rows.append(row)
    out = pd.DataFrame(rows)
    out["type"] = "derived (A_ess, I_pl are stretched)"
    return out


# ----------------------------------------------------------------------------
# RQ5
# ----------------------------------------------------------------------------
def rq5(cfg, units, tsc, pop, bud, floor, observed: pd.DataFrame) -> pd.DataFrame:
    H = cfg["carbon_budget"]["horizon_years"]
    T0 = t0_path(cfg, units)
    obs = observed.set_index(["method", "reference", "A_ess_level"])
    qr = next(m for m in observed.method.unique() if m.startswith("quantile"))
    rows = []
    for ref, A in _refs(cfg, units):
        for _, p in pop.iterrows():
            for target, B in bud.items():
                base = dict(reference=ref, A_ref=A, P_label=p.label, P_persons=p.P_persons,
                            target_degC=target, horizon_yr=H)
                cases = [dict(base, frame="path", T_scenario=n, rate=r, T_start=T0,
                              _solve=lambda lv, r=r: model.required_dem_path(A, floor[lv], B, p.P_persons, T0, r, H))
                         for n, r in rates(cfg).items()]
                cases += [dict(base, frame="static", T_scenario=t.scenario, rate=t.rate, T_start=t["T"],
                               _solve=lambda lv, T=t["T"]: model.required_dem(A, floor[lv], model.allowance(B, H), p.P_persons, T))
                          for _, t in tsc.iterrows()]
                for row in cases:
                    solve = row.pop("_solve")
                    for lv in LEVELS:
                        row[f"A_dem_required_{lv}"], row[f"status_{lv}"] = solve(lv)
                    if ref in cfg["frontier"]["references"]:
                        o, f = obs.loc[(qr, ref, "central")], obs.loc[("FDH peer", ref, "central")]
                        req = row["A_dem_required_central"]
                        row.update(A_dem_observed_frontier=o.A_dem_observed, A_dem_observed_ci_low=o.A_dem_observed_ci_low,
                                   A_dem_observed_ci_high=o.A_dem_observed_ci_high,
                                   A_dem_observed_fdh_peer=f.A_dem_observed, fdh_peer=f.peer,
                                   observed_frontier_meets_required=bool(o.A_dem_observed >= req),
                                   observed_ci_low_meets_required=bool(o.A_dem_observed_ci_low >= req),
                                   observed_ci_high_meets_required=bool(o.A_dem_observed_ci_high >= req))
                    rows.append(row)
    out = pd.DataFrame(rows)
    out["type"] = "derived (A_ess, A_dem, I_pl are stretched)"
    return out


def main():
    import floor as fl
    import frontier
    cfg = load.load_config()
    y = cfg["years"]["activity"]
    units = pd.read_csv(load.PROC / f"units_{y}.csv")
    countries = pd.read_csv(load.PROC / f"countries_{y}.csv")
    pop = pd.read_csv(load.PROC / "population.csv")
    bud = budgets(cfg, pd.read_csv(load.PROC / "carbon_budget.csv"))
    floor = fl.floor_values(pd.read_csv(RES / "rq2_floor.csv"))
    rq6 = pd.read_csv(RES / "rq6_frontier.csv")
    observed = frontier.observed_lookup(rq6)
    ladder = pd.read_csv(load.PROC / "ladder.csv")
    gdp_countries = ladder[ladder.ladder_year == cfg["frontier"]["ladder_year_main"]].reset_index(drop=True)
    tsc = t_scenarios(cfg, units, countries)
    tsc.to_csv(load.PROC / "t_scenarios.csv", index=False)
    rq1(cfg, units, tsc, pop, bud).to_csv(RES / "rq1_burn_rate.csv", index=False)
    obs_main = rq6[(rq6.spec == "main")]
    required_rates(cfg, units, pop, bud, floor, obs_main).to_csv(RES / "rq1_required_rate.csv", index=False)
    rq4(cfg, units, tsc, pop, bud, floor, gdp_countries).to_csv(RES / "rq4_cost_of_existing.csv", index=False)
    rq5(cfg, units, tsc, pop, bud, floor, observed).to_csv(RES / "rq5_required_dematerialism.csv", index=False)


if __name__ == "__main__":
    main()
