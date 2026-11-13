"""RQ2 (floor triangulation) and RQ3 (emissions above the floor).

Conversion used for every route (config: floor.conversion):
    phi   = route floor / current world mean, in the route's native unit and year
    A_ess = phi * A_world                      [2011 int-$ per person per year]
This assumes the native quantity scales proportionally with consumption at today's
average technology; A_ess is therefore typed 'stretched'.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import ndtri

import load

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"


def _v(src: dict, name: str) -> float:
    return float(src["values"][name]["value"])


# ----------------------------------------------------------------------------
# Global distribution of individual footprints (Chancel 2022)
# ----------------------------------------------------------------------------
def chancel_groups(sources: dict) -> pd.DataFrame:
    """Five exhaustive groups: population share and mean footprint [t CO2e/person/yr]."""
    v = sources["chancel_2022"]["values"]
    top10, top1 = v["top10"], v["top1"]
    next9_mean = (top10["pop_share"] * top10["value"] - top1["pop_share"] * top1["value"]) / (
        top10["pop_share"] - top1["pop_share"])
    rows = [("bottom20", v["bottom20"]["pop_share"], v["bottom20"]["value"]),
            ("next30", v["next30"]["pop_share"], v["next30"]["value"]),
            ("middle40", v["middle40"]["pop_share"], v["middle40"]["value"]),
            ("next9", top10["pop_share"] - top1["pop_share"], next9_mean),
            ("top1", top1["pop_share"], top1["value"])]
    g = pd.DataFrame(rows, columns=["group", "pop_share", "mean_tCO2e"])
    g["p_upper"] = g.pop_share.cumsum()
    g["p_lower"] = g.p_upper - g.pop_share
    return g


def _quantiles(family: str, params, u):
    """Quantile function, normalised to mean 1 on the grid u."""
    if family == "lognormal":
        q = np.exp(params[0] * ndtri(u))
    elif family == "singh_maddala":          # Burr XII: Q(u) = ((1-u)^(-1/d) - 1)^(1/c)
        c, d = params
        q = ((1.0 - u) ** (-1.0 / d) - 1.0) ** (1.0 / c)
    else:
        raise ValueError(family)
    return q / q.mean()


def fit_lorenz(groups: pd.DataFrame, family: str, n_grid: int) -> dict:
    """Least-squares fit of a parametric Lorenz curve to the observed cumulative shares.
    Returns the fitted quantile array in t CO2e (mean equal to the group-implied mean)."""
    mean = float((groups.pop_share * groups.mean_tCO2e).sum())
    p_obs = groups.p_upper.values[:-1]
    L_obs = (groups.pop_share * groups.mean_tCO2e).cumsum().values[:-1] / mean
    u = (np.arange(n_grid) + 0.5) / n_grid
    idx = np.round(p_obs * n_grid).astype(int) - 1

    def lorenz(logp):
        q = _quantiles(family, np.exp(logp), u)
        return np.cumsum(q)[idx] / n_grid

    def loss(logp):
        L = lorenz(logp)
        return np.inf if not np.all(np.isfinite(L)) else float(((L - L_obs) ** 2).sum())

    x0 = np.log([1.2] if family == "lognormal" else [1.5, 1.5])
    res = minimize(loss, x0, method="Nelder-Mead", options={"xatol": 1e-8, "fatol": 1e-14, "maxiter": 4000})
    q = _quantiles(family, np.exp(res.x), u) * mean
    return {"family": family, "params": np.exp(res.x), "rmse_lorenz": float(np.sqrt(res.fun / len(L_obs))),
            "u": u, "q_tCO2e": q, "mean_tCO2e": mean, "median_tCO2e": float(np.median(q))}


def fits(cfg: dict, sources: dict) -> list:
    g = chancel_groups(sources)
    n = cfg["floor"]["lorenz_grid"]
    return sorted((fit_lorenz(g, f, n) for f in ("lognormal", "singh_maddala")), key=lambda d: d["rmse_lorenz"])


# ----------------------------------------------------------------------------
# RQ2
# ----------------------------------------------------------------------------
def routes(cfg: dict, sources: dict) -> pd.DataFrame:
    kik, mh, hs = sources["kikstra_2021"], sources["millward_hopkins_2025"], sources["hickel_sullivan_2024"]
    vh = sources["velez_henao_2026"]
    e_now = _v(mh, "current_final_energy")
    p_ref = load.world_population(int(_v(hs, "reference_year")))
    k_hs = p_ref / _v(hs, "population")        # total share for 8.5 bn -> per-person ratio
    f = fits(cfg, sources)
    med = [d["median_tCO2e"] / d["mean_tCO2e"] for d in f]
    rows = [
        dict(route="R1", label="Decent living energy (Kikstra et al. 2021)", lineage="IIASA DLE model",
             native="17 GJ/cap vs 55 GJ/cap final energy",
             phi=_v(kik, "dle_per_cap") / e_now,
             phi_low=_v(kik, "dle_region_low") / e_now, phi_high=_v(kik, "dle_region_high") / e_now,
             range_meaning="regional DLE (South Asia to North America)",
             conversion="DLE per person / current world final energy per person (denominator from Millward-Hopkins et al. 2025 table, 2015)"),
        dict(route="R2", label="Decent living energy (Millward-Hopkins et al. 2025)", lineage="Millward-Hopkins 2022 DLE model",
             native="22 GJ/cap vs 55 GJ/cap final energy",
             phi=_v(mh, "dle_global_avg") / e_now, phi_low=_v(mh, "dle_low") / e_now, phi_high=_v(mh, "dle_high") / e_now,
             range_meaning="national DLE thresholds (119 countries)",
             conversion="global-average DLE / current global average, both from the paper's table"),
        dict(route="R3", label="Decent living materials (Millward-Hopkins et al. 2025)", lineage="Velez-Henao & Pauliuk LCA model",
             native="4.6 t/cap vs 12 t/cap material footprint",
             phi=_v(mh, "dlm") / _v(mh, "current_materials"), phi_low=np.nan, phi_high=np.nan,
             range_meaning="",
             conversion="DLM / current global average material footprint, both from the paper's table"),
        dict(route="R4", label="Needs-based share (Hickel & Sullivan 2024)", lineage="synthesis of the R2 and R3 lineages",
             native="30 % of current resource and energy use for 8.5 bn",
             phi=_v(hs, "share_headline") * k_hs,
             phi_low=min(_v(hs, "share_energy_low"), _v(hs, "share_materials_low")) * k_hs,
             phi_high=max(_v(hs, "share_energy_high"), _v(hs, "share_materials_high")) * k_hs,
             range_meaning="29-44 % across energy and materials, efficient to current technology",
             conversion=f"share x P_2019 / 8.5 bn (P_2019 = {p_ref:.4g} from OWID)"),
        dict(route="R5", label="Decent living GHG, recent technology (Velez-Henao et al. 2026)", lineage="Velez-Henao & Pauliuk LCA model",
             native="29 Gt CO2e vs 54 Gt CO2e (2015)",
             phi=_v(vh, "dls_ghg_recent_tech") / _v(vh, "current_ghg"),
             phi_low=_v(vh, "dls_ghg_efficient") / _v(vh, "current_ghg"), phi_high=np.nan,
             range_meaning="low end = efficient scenario (embeds a change in T; reference only)",
             conversion="global DLS GHG / current global GHG, same year and population"),
        dict(route="R6", label="Global median footprint (Chancel 2022)", lineage="empirical distribution (descriptive, not needs-based)",
             native="median vs mean individual footprint, 2019",
             phi=med[0], phi_low=min(med), phi_high=max(med),
             range_meaning="two fitted Lorenz families",
             conversion=f"median / mean of the fitted distribution (best fit: {f[0]['family']})"),
    ]
    return pd.DataFrame(rows)


def rq2(cfg: dict, sources: dict, units: pd.DataFrame) -> pd.DataFrame:
    w = units.set_index("unit").loc["World"]
    A_w, T_w = w.A_int2011usd_per_cap, w.T_kgCO2_per_int2011usd
    r = routes(cfg, sources)
    r["type"] = "stretched"
    central = pd.DataFrame([dict(
        route="CENTRAL", label="Chosen floor: median of route central values", lineage="",
        native="", phi=r.phi.median(), phi_low=r.phi.min(), phi_high=r.phi.max(),
        range_meaning="min-max of route central values", conversion="config: floor.central_rule, floor.range_rule",
        type="stretched")])
    out = pd.concat([r, central], ignore_index=True)
    for s in ("", "_low", "_high"):
        out[f"A_ess{s}_int2011usd_per_cap"] = out[f"phi{s}"] * A_w
        out[f"floor{s}_tCO2_per_cap_at_T{cfg['years']['activity']}"] = out[f"phi{s}"] * A_w * T_w / 1e3
    return out


def floor_values(rq2_table: pd.DataFrame) -> dict:
    c = rq2_table.set_index("route").loc["CENTRAL"]
    return {"central": c.A_ess_int2011usd_per_cap, "low": c.A_ess_low_int2011usd_per_cap,
            "high": c.A_ess_high_int2011usd_per_cap,
            "phi": {"central": c.phi, "low": c.phi_low, "high": c.phi_high}}


# ----------------------------------------------------------------------------
# RQ3
# ----------------------------------------------------------------------------
def _split(e, f):
    """Per-person means of (floor met, above floor, unmet floor) for footprints e and floor f."""
    return np.minimum(e, f).mean(), np.maximum(e - f, 0).mean(), np.maximum(f - e, 0).mean()


def rq3(cfg: dict, sources: dict, units: pd.DataFrame, floor: dict) -> pd.DataFrame:
    """Share and amount of global CO2 attributable to consumption above the floor.

    Floor in footprint units: f = phi * mean footprint. Shares are applied to observed
    world consumption CO2 of the activity year (the distribution is CO2e, 2019)."""
    g = chancel_groups(sources)
    I_world_Gt = units.set_index("unit").loc["World"].I_kgCO2 / load.KG_PER_GT
    rows = []

    def add(method, bound, level, grp, pop_share, met, above, unmet, mean):
        rows.append(dict(method=method, bound=bound, A_ess_level=level, group=grp, pop_share=pop_share,
                         floor_tCO2e_per_cap=f, mean_tCO2e_per_cap=mean,
                         floor_met_tCO2e_per_cap=met, above_floor_tCO2e_per_cap=above,
                         unmet_floor_tCO2e_per_cap=unmet))

    fitted = fits(cfg, sources)
    for level, phi in floor["phi"].items():
        mean_all = float((g.pop_share * g.mean_tCO2e).sum())
        f = phi * mean_all
        # (a) group means only: exact lower bound on 'above' and on 'unmet' (Jensen)
        parts = [_split(np.array([m]), f) for m in g.mean_tCO2e]
        for (_, gr), (met, above, unmet) in zip(g.iterrows(), parts):
            add("group means", "lower bound on above-floor", level, gr.group, gr.pop_share, met, above, unmet, gr.mean_tCO2e)
        tot = [sum(w * p[i] for w, p in zip(g.pop_share, parts)) for i in range(3)]
        add("group means", "lower bound on above-floor", level, "World", 1.0, *tot, mean_all)
        # (b) fitted continuous distributions: point estimates
        for d in fitted:
            name = f"Lorenz fit: {d['family']}"
            for _, gr in g.iterrows():
                sel = (d["u"] > gr.p_lower) & (d["u"] <= gr.p_upper)
                add(name, "point estimate", level, gr.group, gr.pop_share, *_split(d["q_tCO2e"][sel], f),
                    float(d["q_tCO2e"][sel].mean()))
            add(name, "point estimate", level, "World", 1.0, *_split(d["q_tCO2e"], f), d["mean_tCO2e"])
    out = pd.DataFrame(rows)
    world_mean = out.groupby(["method", "A_ess_level"]).mean_tCO2e_per_cap.transform("last")
    contrib = out.pop_share / world_mean          # converts a group's per-person t into a share of the world total
    out["share_of_world_above_floor"] = out.above_floor_tCO2e_per_cap * contrib
    out["share_of_world_floor_met"] = out.floor_met_tCO2e_per_cap * contrib
    out["unmet_floor_as_share_of_world"] = out.unmet_floor_tCO2e_per_cap * contrib
    out["above_floor_GtCO2"] = out.share_of_world_above_floor * I_world_Gt
    out["floor_met_GtCO2"] = out.share_of_world_floor_met * I_world_Gt
    out["unmet_floor_GtCO2"] = out.unmet_floor_as_share_of_world * I_world_Gt
    rm = {f"Lorenz fit: {d['family']}": d["rmse_lorenz"] for d in fitted}
    out["lorenz_rmse"] = out.method.map(rm)
    out["type"] = "derived (floor is stretched)"
    return out


def main():
    cfg, sources = load.load_config(), load.load_sources()
    units = pd.read_csv(load.PROC / f"units_{cfg['years']['activity']}.csv")
    RES.mkdir(exist_ok=True)
    t2 = rq2(cfg, sources, units)
    t2.to_csv(RES / "rq2_floor.csv", index=False)
    rq3(cfg, sources, units, floor_values(t2)).to_csv(RES / "rq3_materialism.csv", index=False)


if __name__ == "__main__":
    main()
