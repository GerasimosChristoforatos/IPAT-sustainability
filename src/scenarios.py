"""RQ1 (burn rate), RQ4 (cost of existing, equal shares) and RQ5 (required dematerialism).

All scenarios are counterfactuals: the whole world population lives at the reference
A (US, EU27 or World GDP per capita of the activity year) with a given T. They are not
forecasts. T rates are borrowed from a production-based metric (see sources.yaml).
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


def t_scenarios(cfg: dict, units: pd.DataFrame, countries: pd.DataFrame) -> pd.DataFrame:
    """T levels [kg CO2 / 2011 int-$]: today, three 2050 endpoints, and national anchors."""
    w = units.set_index("unit").loc["World"]
    T0, span = w[T_COL], cfg["years"]["t_scenario_end"] - cfg["years"]["activity"]
    rows = [dict(scenario="Today", kind="observed", rate=0.0, T=T0, note=f"world T, {cfg['years']['activity']}", type="derived")]
    for name, spec in cfg["t_scenarios"].items():
        rows.append(dict(scenario=name, kind="rate endpoint", rate=spec["rate"],
                         T=model.t_endpoint(T0, spec["rate"], span),
                         note=f"{spec['rate']:.1%}/yr to {cfg['years']['t_scenario_end']}", type="derived (rate is an assumption)"))
    rule = cfg["best_national_t"]
    elig = countries[(countries.P_persons >= rule["min_population"]) &
                     (countries[A_COL] >= rule["min_a_relative_to_world"] * w[A_COL])].sort_values(T_COL)
    for rank, label in ((0, "Best national"), (1, "Second-best national")):
        c = elig.iloc[rank]
        rows.append(dict(scenario=label, kind="national anchor", rate=np.nan, T=c[T_COL],
                         note=f"{c.country} (A >= world mean, population >= {rule['min_population']:,})", type="derived"))
    if rule["report_unrestricted"]:
        c = countries.sort_values(T_COL).iloc[0]
        rows.append(dict(scenario="Lowest national, unrestricted", kind="national anchor", rate=np.nan,
                         T=c[T_COL], note=c.country, type="derived"))
    return pd.DataFrame(rows)


def budgets(cfg: dict, budget: pd.DataFrame) -> dict:
    """Remaining budget [kg CO2] by target at the central likelihood."""
    b = budget[np.isclose(budget.likelihood, cfg["carbon_budget"]["likelihood_central"])]
    return {float(r.target_degC): r.budget_GtCO2 * model.KG_PER_GT for _, r in b.iterrows()}


def _grid(cfg, units, tsc, pop, bud):
    u = units.set_index("unit")
    for ref in cfg["units"]:
        for _, t in tsc.iterrows():
            for _, p in pop.iterrows():
                for target, B in bud.items():
                    yield ref, u.loc[ref, A_COL], t, p, target, B


def rq1(cfg, units, tsc, pop, bud) -> pd.DataFrame:
    cb, yrs = cfg["carbon_budget"], cfg["years"]
    T0 = units.set_index("unit").loc["World", T_COL]
    rows = []
    for ref, A, t, p, target, B in _grid(cfg, units, tsc, pop, bud):
        I = p.P_persons * A * t["T"]
        row = dict(reference=ref, A_ref=A, T_scenario=t.scenario, T=t["T"], P_label=p.label, P_persons=p.P_persons,
                   target_degC=target, budget_GtCO2=B / model.KG_PER_GT, I_GtCO2_per_yr=I / model.KG_PER_GT)
        for h in [cb["horizon_years"], *cb["horizon_sensitivity_years"]]:
            row[f"I_over_Ipl_h{h}"] = I / model.allowance(B, h)
        row["years_to_exhaustion_static"] = model.years_to_exhaustion_static(B, I)
        # Path-consistent counterfactual: everyone at A_ref from now, T declining from today's
        # observed level at the scenario rate; emissions counted from the budget base year.
        if t.kind == "rate endpoint":
            I_start = p.P_persons * A * model.t_endpoint(T0, t.rate, yrs["budget_base"] - yrs["activity"])
            y = model.years_to_exhaustion_path(B, I_start, t.rate)
            row["years_to_exhaustion_path"] = y
            row["cumulative_bound_GtCO2_path"] = I_start / model.decay_constant(t.rate) / model.KG_PER_GT
        elif t.kind == "observed":
            row["years_to_exhaustion_path"] = row["years_to_exhaustion_static"]
            row["cumulative_bound_GtCO2_path"] = math.inf
        else:
            row["years_to_exhaustion_path"] = np.nan
            row["cumulative_bound_GtCO2_path"] = np.nan
        y = row["years_to_exhaustion_path"]
        row["path_outcome"] = ("not applicable" if isinstance(y, float) and math.isnan(y)
                               else "budget never exhausted" if math.isinf(y)
                               else f"exhausted in {yrs['budget_base'] + y:.0f}")
        row["type"] = "derived (counterfactual; I_pl is stretched)"
        rows.append(row)
    return pd.DataFrame(rows)


def _nearest_country(A_eq: float, gdp_countries: pd.DataFrame) -> str:
    """Country with the closest GDP per capita (log distance); flagged if out of range."""
    if A_eq > gdp_countries[A_COL].max():
        return "above every country in the sample"
    if A_eq < gdp_countries[A_COL].min():
        return "below every country in the sample"
    i = (np.log(gdp_countries[A_COL]) - math.log(A_eq)).abs().idxmin()
    return gdp_countries.loc[i, "country"]


def rq4(cfg, units, tsc, pop, bud, floor, gdp_countries) -> pd.DataFrame:
    """Cost of existing P*A_ess*T against I_pl; equal-share A and the country it resembles."""
    h = cfg["carbon_budget"]["horizon_years"]
    seen, rows = set(), []
    for _, _, t, p, target, B in _grid(cfg, units, tsc, pop, bud):
        key = (t.scenario, p.label, target)
        if key in seen:
            continue
        seen.add(key)
        I_pl = model.allowance(B, h)
        A_eq = model.equal_share_A(I_pl, p.P_persons, t["T"])
        row = dict(T_scenario=t.scenario, T=t["T"], P_label=p.label, P_persons=p.P_persons, target_degC=target,
                   horizon_yr=h, I_pl_GtCO2_per_yr=I_pl / model.KG_PER_GT)
        for level in ("central", "low", "high"):
            c = model.cost_of_existing(p.P_persons, floor[level], t["T"])
            row[f"cost_of_existing_GtCO2_per_yr_{level}"] = c / model.KG_PER_GT
            row[f"cost_over_Ipl_{level}"] = c / I_pl
            row[f"floor_within_allowance_{level}"] = bool(c <= I_pl)
        row["equal_share_A_int2011usd_per_cap"] = A_eq
        row["equal_share_A_over_A_ess_central"] = A_eq / floor["central"]
        row["equal_share_lives_like"] = _nearest_country(A_eq, gdp_countries)
        row["type"] = "derived (A_ess, I_pl are stretched)"
        rows.append(row)
    return pd.DataFrame(rows)


def rq5(cfg, units, tsc, pop, bud, floor, observed: pd.DataFrame) -> pd.DataFrame:
    h = cfg["carbon_budget"]["horizon_years"]
    obs = observed.set_index(["method", "reference", "A_ess_level"])
    qr = next(m for m in observed.method.unique() if m.startswith("quantile"))
    rows = []
    for ref, A, t, p, target, B in _grid(cfg, units, tsc, pop, bud):
        I_pl = model.allowance(B, h)
        row = dict(reference=ref, A_ref=A, T_scenario=t.scenario, T=t["T"], P_label=p.label, P_persons=p.P_persons,
                   target_degC=target, horizon_yr=h)
        for level in ("central", "low", "high"):
            v, status = model.required_dem(A, floor[level], I_pl, p.P_persons, t["T"])
            row[f"A_dem_required_{level}"] = v
            row[f"status_{level}"] = status
        if ref in cfg["frontier"]["references"]:
            o = obs.loc[(qr, ref, "central")]
            row["A_dem_observed_frontier"] = o.A_dem_observed
            row["A_dem_observed_ci_low"] = o.A_dem_observed_ci_low
            row["A_dem_observed_ci_high"] = o.A_dem_observed_ci_high
            f = obs.loc[("FDH peer", ref, "central")]
            row["A_dem_observed_fdh_peer"] = f.A_dem_observed
            row["fdh_peer"] = f.peer
            req = row["A_dem_required_central"]
            row["observed_frontier_meets_required"] = bool(o.A_dem_observed >= req)
            row["observed_ci_low_meets_required"] = bool(o.A_dem_observed_ci_low >= req)
        row["type"] = "derived (A_ess, A_dem, I_pl are stretched)"
        rows.append(row)
    return pd.DataFrame(rows)


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
    observed = frontier.observed_lookup(pd.read_csv(RES / "rq6_frontier.csv"))
    ladder = pd.read_csv(load.PROC / "ladder.csv")
    gdp_countries = ladder[ladder.ladder_year == cfg["frontier"]["ladder_year_main"]].reset_index(drop=True)
    tsc = t_scenarios(cfg, units, countries)
    tsc.to_csv(load.PROC / "t_scenarios.csv", index=False)
    rq1(cfg, units, tsc, pop, bud).to_csv(RES / "rq1_burn_rate.csv", index=False)
    rq4(cfg, units, tsc, pop, bud, floor, gdp_countries).to_csv(RES / "rq4_cost_of_existing.csv", index=False)
    rq5(cfg, units, tsc, pop, bud, floor, observed).to_csv(RES / "rq5_required_dematerialism.csv", index=False)


if __name__ == "__main__":
    main()
