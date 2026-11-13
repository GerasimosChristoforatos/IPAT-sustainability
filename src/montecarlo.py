"""One Monte Carlo propagation for the required dematerialism A_dem* (US, EU27 references).

Draws (common random numbers across scenarios, fixed seed):
    A_ess    triangular(min, central, max) of the RQ2 floor
    budget   likelihood u ~ U(0.33, 0.67); budget by linear interpolation of the
             published IGCC quantiles (higher likelihood = smaller budget)
    horizon  triangular(25, 50, 100) years
    T        fixed at each scenario's 2050 endpoint, plus one pooled row with the
             decline rate ~ U(BAU rate, Extreme rate)

Horizon and T rate are choices, not measurements: the result is a scenario-weighted
distribution, not a probability forecast. A_dem* lives on the extended real line:
'floor fails' = +inf, 'none needed' = A_dem* <= 1, so quantiles are well defined.
"""
from pathlib import Path

import numpy as np
import pandas as pd

import load
import model

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
A_COL, T_COL = "A_int2011usd_per_cap", "T_kgCO2_per_int2011usd"


def required_dem_vec(A_ref, A_ess, I_pl, P, T):
    """Vectorised model.required_dem (value only; +inf where the floor fails)."""
    headroom = I_pl / (P * T) - A_ess
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(headroom > 0, (A_ref - A_ess) / headroom, np.inf)


def run(cfg, units, budget, pop, floor, observed: dict) -> pd.DataFrame:
    mc, yrs = cfg["montecarlo"], cfg["years"]
    n = mc["draws"]
    rng = np.random.default_rng(mc["seed"])
    A_ess = rng.triangular(floor["low"], floor["central"], floor["high"], n)
    u = rng.uniform(*cfg["carbon_budget"]["likelihood_range"], n)
    lo, mode, hi = mc["horizon_triangular_years"]
    horizon = rng.triangular(lo, mode, hi, n)
    rates = {k: v["rate"] for k, v in cfg["t_scenarios"].items()}
    pooled_rate = rng.uniform(min(rates.values()), max(rates.values()), n)

    uw = units.set_index("unit")
    T0, span = uw.loc["World", T_COL], yrs["t_scenario_end"] - yrs["activity"]
    P = float(pop.set_index("label").loc["now", "P_persons"])
    T_draws = {name: np.full(n, model.t_endpoint(T0, r, span)) for name, r in rates.items()}
    if mc["pooled_rate_uniform"]:
        T_draws["Pooled (rate uniform)"] = T0 * (1.0 - pooled_rate) ** span

    q_lo, q_hi = mc["interval"]
    rows = []
    for target in cfg["carbon_budget"]["targets_degC"]:
        b = budget[budget.target_degC == target].sort_values("likelihood")
        B = np.interp(u, b.likelihood.values, b.budget_GtCO2.values) * model.KG_PER_GT
        I_pl = B / horizon
        for ref in cfg["frontier"]["references"]:
            for name, T in T_draws.items():
                d = required_dem_vec(uw.loc[ref, A_COL], A_ess, I_pl, P, T)
                rows.append(dict(
                    reference=ref, target_degC=target, T_scenario=name, draws=n, P_label="now",
                    # empirical-CDF quantiles (no interpolation): well defined with +inf draws
                    A_dem_required_median=np.quantile(d, 0.5, method="inverted_cdf"),
                    A_dem_required_p05=np.quantile(d, q_lo, method="inverted_cdf"),
                    A_dem_required_p95=np.quantile(d, q_hi, method="inverted_cdf"),
                    prob_required_le_observed=float(np.mean(d <= observed[ref])),
                    prob_floor_fails=float(np.mean(np.isinf(d))),
                    prob_none_needed=float(np.mean(d <= 1)),
                    type="derived (scenario-weighted; A_ess, A_dem, I_pl are stretched)"))
    return pd.DataFrame(rows)


def main():
    import floor as fl
    cfg = load.load_config()
    units = pd.read_csv(load.PROC / f"units_{cfg['years']['activity']}.csv")
    budget = pd.read_csv(load.PROC / "carbon_budget.csv")
    pop = pd.read_csv(load.PROC / "population.csv")
    floor = fl.floor_values(pd.read_csv(RES / "rq2_floor.csv"))
    obs = pd.read_csv(RES / "rq6_frontier.csv")
    obs = obs[(obs.spec == "main") & (obs.A_ess_level == "central") & obs.method.str.startswith("quantile")]
    # observed A_dem is held at its central point estimate (its own uncertainty is in RQ6)
    out = run(cfg, units, budget, pop, floor, dict(zip(obs.reference, obs.A_dem_observed)))
    out = out.merge(obs[["reference", "A_dem_observed", "A_dem_observed_ci_low", "A_dem_observed_ci_high"]],
                    on="reference", how="left")
    out.to_csv(RES / "montecarlo_required_dematerialism.csv", index=False)


if __name__ == "__main__":
    main()
