"""RQ6: empirical wellbeing frontier and the dematerialism already observed.

Primary estimator: linear quantile regression of the national Cantril ladder on
log GDP per capita at tau = 0.9 (Koenker & Bassett 1978), pairs bootstrap over countries.
Secondary: free-disposal-hull peer (Deprins, Simar & Tulkens 1984): the lowest-A country
whose ladder is at least the reference's ladder.

Cross-sectional association only: it shows a combination that exists, not what an
intervention would cause. GDP per capita is a proxy for consumption, not for wanting.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.optimize import linprog

import load
import model

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
A_COL = "A_int2011usd_per_cap"


def sample(ladder: pd.DataFrame, year: int) -> pd.DataFrame:
    return ladder[ladder.ladder_year == year].reset_index(drop=True)


def references(cfg: dict, units: pd.DataFrame, s: pd.DataFrame) -> pd.DataFrame:
    """Ladder and A for the reference units. EU27 ladder = population-weighted member mean."""
    u = units.set_index("unit")
    usa = s[s.iso_code == "USA"].iloc[0]
    eu = s[s.iso_code.isin(cfg["eu27_iso"])]
    if len(eu) != len(cfg["eu27_iso"]):
        raise ValueError(f"EU27 ladder incomplete: {len(eu)} members")
    rows = [("US", usa.ladder, u.loc["US", A_COL]),
            ("EU27", float(np.average(eu.ladder, weights=eu.P_persons)), u.loc["EU27", A_COL])]
    return pd.DataFrame(rows, columns=["reference", "ladder", A_COL])


def fit_qr(s: pd.DataFrame, tau: float, quadratic: bool = False):
    x = np.log(s[A_COL].values)
    X = np.column_stack([np.ones_like(x), x] + ([x ** 2] if quadratic else []))
    return sm.QuantReg(s.ladder.values, X).fit(q=tau)


def qr_params(x_log: np.ndarray, y: np.ndarray, tau: float) -> np.ndarray:
    """Exact linear quantile regression (intercept, slope) by linear programming:
    minimise tau*sum(u+) + (1-tau)*sum(u-)  s.t.  X b + u+ - u- = y,  u+, u- >= 0.
    Used for point estimates and the bootstrap; avoids the iteration limit of the
    iteratively reweighted solver on resamples with ties."""
    n = len(y)
    X = np.column_stack([np.ones(n), x_log])
    c = np.concatenate([np.zeros(2), tau * np.ones(n), (1 - tau) * np.ones(n)])
    A_eq = np.hstack([X, np.eye(n), -np.eye(n)])
    res = linprog(c, A_eq=A_eq, b_eq=y, bounds=[(None, None)] * 2 + [(0, None)] * (2 * n), method="highs")
    if not res.success:
        raise RuntimeError(res.message)
    return res.x[:2]


def frontier_A(alpha: float, beta: float, ladder_ref: float) -> float:
    """A at which the fitted frontier reaches ladder_ref. inf if the slope is not positive."""
    return float(np.exp((ladder_ref - alpha) / beta)) if beta > 0 else np.inf


def estimate(s, refs, tau, A_ess_levels, reps, seed, spec):
    """Point estimates and bootstrap percentile intervals for one specification."""
    a, b = qr_params(np.log(s[A_COL].values), s.ladder.values, tau)
    rng = np.random.default_rng(seed)
    boot = np.empty((reps, 2))
    n = len(s)
    for i in range(reps):
        bs = s.iloc[rng.integers(0, n, n)]
        boot[i] = qr_params(np.log(bs[A_COL].values), bs.ladder.values, tau)
    rows = []
    for _, r in refs.iterrows():
        A_f = frontier_A(a, b, r.ladder)
        A_f_boot = np.array([frontier_A(x, y, r.ladder) for x, y in boot])
        for level, A_ess in A_ess_levels.items():
            obs_boot = np.array([model.observed_dem(r[A_COL], v, A_ess) for v in A_f_boot])
            rows.append(dict(
                spec=spec, method=f"quantile regression tau={tau}", reference=r.reference,
                ladder_year=int(s.ladder_year.iloc[0]), n_countries=n,
                ladder_ref=r.ladder, A_ref=r[A_COL], A_ess_level=level, A_ess=A_ess,
                intercept=a, slope=b, slope_ci_low=np.percentile(boot[:, 1], 2.5),
                slope_ci_high=np.percentile(boot[:, 1], 97.5),
                A_frontier=A_f, A_frontier_ci_low=np.percentile(A_f_boot, 2.5),
                A_frontier_ci_high=np.percentile(A_f_boot, 97.5),
                A_dem_observed=model.observed_dem(r[A_COL], A_f, A_ess),
                A_dem_observed_ci_low=np.percentile(obs_boot, 2.5),
                A_dem_observed_ci_high=np.percentile(obs_boot, 97.5),
                share_boot_frontier_below_floor=float(np.mean(~np.isfinite(obs_boot))),
                peer="", peer_ladder=np.nan, peer_ladder_lower=np.nan))
    return rows


def fdh(s, refs, A_ess_levels, spec):
    """Free-disposal-hull peer: lowest-A country with ladder >= reference ladder."""
    rows = []
    for _, r in refs.iterrows():
        cand = s[s.ladder >= r.ladder]
        p = cand.loc[cand[A_COL].idxmin()]
        for level, A_ess in A_ess_levels.items():
            rows.append(dict(
                spec=spec, method="FDH peer", reference=r.reference, ladder_year=int(s.ladder_year.iloc[0]),
                n_countries=len(s), ladder_ref=r.ladder, A_ref=r[A_COL], A_ess_level=level, A_ess=A_ess,
                A_frontier=p[A_COL], A_dem_observed=model.observed_dem(r[A_COL], p[A_COL], A_ess),
                peer=p.country, peer_ladder=p.ladder, peer_ladder_lower=p.ladder_lower))
    return rows


def rq6(cfg: dict, units: pd.DataFrame, ladder: pd.DataFrame, floor: dict) -> pd.DataFrame:
    f = cfg["frontier"]
    levels = {k: floor[k] for k in ("central", "low", "high")}
    main = sample(ladder, f["ladder_year_main"])
    rob = sample(ladder, f["ladder_year_robustness"])
    refs_main, refs_rob = references(cfg, units, main), references(cfg, units, rob)
    rows = estimate(main, refs_main, f["tau"], levels, f["bootstrap_reps"], cfg["seed"], "main")
    rows += fdh(main, refs_main, levels, "main")
    rows += estimate(main, refs_main, f["tau_robustness"], levels, f["bootstrap_reps"], cfg["seed"], "robustness: tau")
    rows += estimate(rob, refs_rob, f["tau"], levels, f["bootstrap_reps"], cfg["seed"], "robustness: ladder window")
    out = pd.DataFrame(rows)
    quad = fit_qr(main, f["tau"], quadratic=True)
    out["quadratic_term_coef"] = quad.params[2]
    out["quadratic_term_pvalue"] = quad.pvalues[2]
    ols = sm.OLS(main.ladder.values, sm.add_constant(np.log(main[A_COL].values))).fit()
    out["ols_slope_for_benchmark"] = ols.params[1]
    out["type"] = "derived (A_dem is stretched)"
    return out


def observed_lookup(rq6_table: pd.DataFrame) -> pd.DataFrame:
    m = rq6_table[rq6_table.spec == "main"]
    return m[["method", "reference", "A_ess_level", "A_frontier", "A_dem_observed",
              "A_dem_observed_ci_low", "A_dem_observed_ci_high", "peer"]]


def main():
    cfg = load.load_config()
    import floor as fl
    units = pd.read_csv(load.PROC / f"units_{cfg['years']['activity']}.csv")
    ladder = pd.read_csv(load.PROC / "ladder.csv")
    floor = fl.floor_values(pd.read_csv(RES / "rq2_floor.csv"))
    rq6(cfg, units, ladder, floor).to_csv(RES / "rq6_frontier.csv", index=False)


if __name__ == "__main__":
    main()
