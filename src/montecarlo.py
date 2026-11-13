"""One Monte Carlo propagation for the required dematerialism A_dem* (US, EU27), in the
path frame (main) and the static frame (robustness), with common random numbers.

Draws (fixed seed):
    A_ess    triangular(min, central, max) of the RQ2 floor
    budget   likelihood u ~ U(0.33, 0.67); budget by linear interpolation of the
             published IGCC quantiles (higher likelihood = smaller budget)
    horizon  triangular(25, 50, 100) years
    T        each scenario's rate held fixed, plus one pooled row with the rate
             ~ U(BAU rate, Extreme rate)

Path frame:   A_dem* = A_mat / ( B / (P T0 F(r, H)) - A_ess ),  F = (1 - e^(-kH))/k
Static frame: A_dem* = A_mat / ( (B/H) / (P T_2050(r)) - A_ess )

Horizon and T rate are choices, not measurements: the result is a scenario-weighted
distribution, not a probability forecast. A_dem* lives on the extended real line
('floor fails' = +inf), so empirical-CDF quantiles are used.
"""
from pathlib import Path

import numpy as np
import pandas as pd

import load
import model
import scenarios

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
A_COL, T_COL = "A_int2011usd_per_cap", "T_kgCO2_per_int2011usd"


def required_dem_vec(A_ref, A_ess, I_pl, P, T):
    """Vectorised model.required_dem (value only; +inf where the floor fails)."""
    headroom = I_pl / (P * T) - A_ess
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(headroom > 0, (A_ref - A_ess) / headroom, np.inf)


def path_factor_vec(rate, H):
    """Vectorised model.path_factor."""
    rate, H = np.broadcast_arrays(np.asarray(rate, float), np.asarray(H, float))
    k = -np.log1p(-rate)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(k > 0, -np.expm1(-k * H) / np.where(k > 0, k, 1), H)


def draws(cfg, units, budget, pop, floor) -> dict:
    """All Monte Carlo draws: {(frame, reference, target, scenario): array of A_dem*}."""
    mc, yrs = cfg["montecarlo"], cfg["years"]
    n = mc["draws"]
    rng = np.random.default_rng(mc["seed"])
    A_ess = rng.triangular(floor["low"], floor["central"], floor["high"], n)
    u = rng.uniform(*cfg["carbon_budget"]["likelihood_range"], n)
    lo, mode, hi = mc["horizon_triangular_years"]
    H = rng.triangular(lo, mode, hi, n)
    rts = {k: v["rate"] for k, v in cfg["t_scenarios"].items()}
    pooled = rng.uniform(min(rts.values()), max(rts.values()), n)
    rate_draws = {name: np.full(n, r) for name, r in rts.items()}
    if mc["pooled_rate_uniform"]:
        rate_draws["Pooled (rate uniform)"] = pooled

    uw = units.set_index("unit")
    T_now, span = uw.loc["World", T_COL], yrs["t_scenario_end"] - yrs["activity"]
    T0 = scenarios.t0_path(cfg, units)
    P = float(pop.set_index("label").loc["now", "P_persons"])
    out = {}
    for target in cfg["carbon_budget"]["targets_degC"]:
        b = budget[budget.target_degC == target].sort_values("likelihood")
        B = np.interp(u, b.likelihood.values, b.budget_GtCO2.values) * model.KG_PER_GT
        for ref in cfg["frontier"]["references"]:
            A = uw.loc[ref, A_COL]
            for name, r in rate_draws.items():
                out[("path", ref, target, name)] = required_dem_vec(A, A_ess, B / path_factor_vec(r, H), P, T0)
                out[("static", ref, target, name)] = required_dem_vec(A, A_ess, B / H, P, T_now * (1 - r) ** span)
    return out


def summarise(d: dict, cfg: dict, observed: dict) -> pd.DataFrame:
    q_lo, q_hi = cfg["montecarlo"]["interval"]
    rows = []
    for (frame, ref, target, name), x in d.items():
        rows.append(dict(
            frame=frame, reference=ref, target_degC=target, T_scenario=name, draws=len(x), P_label="now",
            A_dem_required_median=np.quantile(x, 0.5, method="inverted_cdf"),
            A_dem_required_p05=np.quantile(x, q_lo, method="inverted_cdf"),
            A_dem_required_p95=np.quantile(x, q_hi, method="inverted_cdf"),
            prob_required_le_observed=float(np.mean(x <= observed[ref])),
            prob_floor_fails=float(np.mean(np.isinf(x))),
            prob_none_needed=float(np.mean(x <= 1)),
            A_dem_observed=observed[ref],
            type="derived (scenario-weighted; A_ess, A_dem, I_pl are stretched)"))
    out = pd.DataFrame(rows)
    order = {"path": 0, "static": 1}
    return out.sort_values(["frame", "reference", "target_degC"], key=lambda s: s.map(order) if s.name == "frame" else s,
                           kind="stable").reset_index(drop=True)


def observed_central(rq6: pd.DataFrame) -> dict:
    """Observed A_dem held at its central point estimate (its own uncertainty is in RQ6)."""
    o = rq6[(rq6.spec == "main") & (rq6.A_ess_level == "central") & rq6.method.str.startswith("quantile")]
    return dict(zip(o.reference, o.A_dem_observed))


def main():
    import floor as fl
    cfg = load.load_config()
    units = pd.read_csv(load.PROC / f"units_{cfg['years']['activity']}.csv")
    budget = pd.read_csv(load.PROC / "carbon_budget.csv")
    pop = pd.read_csv(load.PROC / "population.csv")
    floor = fl.floor_values(pd.read_csv(RES / "rq2_floor.csv"))
    obs = observed_central(pd.read_csv(RES / "rq6_frontier.csv"))
    summarise(draws(cfg, units, budget, pop, floor), cfg, obs).to_csv(
        RES / "montecarlo_required_dematerialism.csv", index=False)


if __name__ == "__main__":
    main()
