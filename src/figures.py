"""The four figures (PNG + SVG). Reads results/ and data/processed/ only."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import frontier
import load
import model

ROOT = Path(__file__).resolve().parents[1]
RES, FIG = ROOT / "results", ROOT / "figures"
A_COL, T_COL = "A_int2011usd_per_cap", "T_kgCO2_per_int2011usd"

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"   # categorical slots 1-3, fixed order
INK, INK2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"

plt.rcParams.update({
    "font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
    "axes.facecolor": SURFACE, "figure.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
    "legend.frameon": False, "svg.fonttype": "none",
})


def _save(fig, name):
    FIG.mkdir(exist_ok=True)
    for ext in ("png", "svg"):
        fig.savefig(FIG / f"{name}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def fig1_feasibility(cfg, units, tsc, pop, budget, floor, rq6):
    u = units.set_index("unit")
    P = float(pop.set_index("label").loc["now", "P_persons"])
    h = cfg["carbon_budget"]["horizon_years"]
    b50 = budget[np.isclose(budget.likelihood, 0.5)].set_index("target_degC").budget_GtCO2
    target = 1.7
    I_pl = model.allowance(b50[target] * model.KG_PER_GT, h)
    A_ess = floor["central"]
    T = np.logspace(np.log10(0.006), np.log10(0.45), 600)

    def boundary(A_ref, I_allow):
        head = I_allow / (P * T) - A_ess
        return np.where(head > 0, (A_ref - A_ess) / np.where(head > 0, head, 1), np.nan)

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    y_top = 200
    us = boundary(u.loc["US", A_COL], I_pl)
    T_fail = I_pl / (P * A_ess)
    ax.fill_between(T, np.clip(np.nan_to_num(us, nan=y_top), 0.5, y_top), y_top, where=T < T_fail,
                    color=BLUE, alpha=0.14, linewidth=0)
    for ref, color in (("US", BLUE), ("EU27", ORANGE), ("World", AQUA)):
        ax.plot(T, boundary(u.loc[ref, A_COL], I_pl), color=color, lw=2, label=f"{ref} lifestyle")
    for tg, ls in ((1.5, ":"), (2.0, "--")):
        ax.plot(T, boundary(u.loc["US", A_COL], model.allowance(b50[tg] * model.KG_PER_GT, h)),
                color=BLUE, lw=1, ls=ls)
    ax.axvline(T_fail, color=INK2, lw=0.8, ls="-.")
    ax.text(T_fail * 1.04, 110, "floor alone\nexceeds Iᵖˡ", fontsize=7.5, color=INK2, va="top")
    ax.axhline(1, color=INK2, lw=0.8)
    ax.text(0.0063, 0.86, "below this line: no dematerialism needed", fontsize=7.5, color=INK2, va="top")

    short = {"Today": "Today", "BAU": "BAU 2050", "Ambitious": "Ambitious 2050", "Extreme": "Extreme 2050",
             "Best national": "Norway", "Second-best national": "Sweden"}
    for _, t in tsc[tsc.scenario.isin(short)].iterrows():
        ax.axvline(t["T"], color=GRID, lw=1.2, zorder=0)
        ax.text(t["T"], y_top * 1.08, short[t.scenario], rotation=40, fontsize=7.5, ha="left", va="bottom", color=INK2)
        v, status = model.required_dem(u.loc["US", A_COL], A_ess, I_pl, P, t["T"])
        if status == model.REQUIRED and v < y_top:
            ax.plot(t["T"], v, "o", ms=6, color=BLUE, mec=SURFACE, mew=1.2, zorder=5)
            ax.annotate(f"{v:.1f}", (t["T"], v), textcoords="offset points", xytext=(6, 4), fontsize=7.5)

    o = rq6[(rq6.spec == "main") & (rq6.A_ess_level == "central") & rq6.method.str.startswith("quantile")].set_index("reference")
    for ref, color in (("US", BLUE), ("EU27", ORANGE)):
        r = o.loc[ref]
        ax.axhspan(r.A_dem_observed_ci_low, r.A_dem_observed_ci_high, color=color, alpha=0.07, linewidth=0)
        ax.axhline(r.A_dem_observed, color=color, lw=1.2, ls=(0, (5, 2)))
        ax.text(0.44, r.A_dem_observed * (1.07 if ref == "US" else 0.80), f"observed Aᵈᵉᵐ, {ref}: {r.A_dem_observed:.1f}",
                fontsize=7.5, ha="right")
    ax.text(0.0075, 30, "feasible for the\nUS lifestyle\n(I ≤ Iᵖˡ)", fontsize=8.5, color=INK2)
    ax.set(xscale="log", yscale="log", xlim=(0.006, 0.45), ylim=(0.5, y_top),
           xlabel="T, carbon intensity of consumption (kg CO₂ per 2011 int-$, log scale)",
           ylabel="Aᵈᵉᵐ, dematerialism (log scale)")
    ax.legend(loc="lower right", fontsize=8, title=f"Required Aᵈᵉᵐ*, {target} °C (50 %), {h}-yr horizon",
              title_fontsize=7.5, frameon=True, facecolor=SURFACE, edgecolor="none", framealpha=0.92)
    fig.text(0.125, -0.02, "Thin blue lines: US lifestyle at 1.5 °C (dotted) and 2.0 °C (dashed). Shaded bands: 95 % bootstrap "
             "interval of observed Aᵈᵉᵐ.\nCounterfactual: world population of 2022 at the reference GDP per capita. "
             "Aᵉˢˢ, Aᵈᵉᵐ and Iᵖˡ are stretched concepts.", fontsize=7, color=INK2, va="top")
    _save(fig, "fig1_feasibility_map")


def fig2_floor(rq2):
    r = rq2[rq2.route != "CENTRAL"].iloc[::-1].reset_index(drop=True)
    c = rq2[rq2.route == "CENTRAL"].iloc[0]
    k = "A_ess{}_int2011usd_per_cap"
    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    ax.axvspan(c[k.format("_low")], c[k.format("_high")], color=BLUE, alpha=0.10, linewidth=0)
    ax.axvline(c[k.format("")], color=BLUE, lw=2)
    for i, row in r.iterrows():
        lo, hi = row[k.format("_low")], row[k.format("_high")]
        lo = row[k.format("")] if np.isnan(lo) else lo
        hi = row[k.format("")] if np.isnan(hi) else hi
        ax.plot([lo, hi], [i, i], color=INK2, lw=1.5, solid_capstyle="round")
        ax.plot(row[k.format("")], i, "o", ms=7, color=BLUE, mec=SURFACE, mew=1.2, zorder=3)
    ax.set_yticks(range(len(r)))
    ax.set_yticklabels([f"{a}  {b}" for a, b in zip(r.route, r.label)], fontsize=8, color=INK)
    ax.set_xlabel("Aᵉˢˢ (2011 int-$ per person per year)")
    A_w = c[k.format("")] / c.phi
    sec = ax.secondary_xaxis("top", functions=(lambda x: x / A_w, lambda p: p * A_w))
    sec.set_xlabel("φ = floor / current world mean", fontsize=8)
    ax.text(c[k.format("")], len(r) - 0.45, f" central ${c[k.format('')]:,.0f}", fontsize=8, va="bottom")
    ax.grid(axis="y", visible=False)
    ax.set_ylim(-0.6, len(r) - 0.2)
    fig.text(0.02, -0.04, "Dots: route central value. Bars: each route's own range (regional or scenario spread, not a confidence "
             "interval).\nBand: min–max of route centrals. R5's low end is the efficient scenario, which embeds a change in T.",
             fontsize=7, color=INK2, va="top")
    _save(fig, "fig2_floor_triangulation")


def fig3_anatomy(rq3):
    d = rq3[(rq3.A_ess_level == "central") & (rq3.method == "Lorenz fit: singh_maddala")].set_index("group")
    order = ["World", "top1", "next9", "middle40", "next30", "bottom20"]
    names = {"World": "World", "top1": "Top 1 %", "next9": "Next 9 %", "middle40": "Middle 40 %",
             "next30": "Next 30 %", "bottom20": "Bottom 20 %"}
    d = d.loc[order[::-1]]
    y = np.arange(len(d)) + np.array([0, 0, 0, 0, 0, 0.5])
    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    ax.barh(y, d.floor_met_GtCO2, height=0.55, color=BLUE, edgecolor=SURFACE, linewidth=1.5, label="Existing: at or below the floor")
    ax.barh(y, d.above_floor_GtCO2, left=d.floor_met_GtCO2, height=0.55, color=ORANGE, edgecolor=SURFACE, linewidth=1.5,
            label="Materialism: above the floor")
    for yi, (_, r) in zip(y, d.iterrows()):
        tot = r.floor_met_GtCO2 + r.above_floor_GtCO2
        note = f"{r.above_floor_GtCO2:.1f} Gt above"
        if r.unmet_floor_GtCO2 > 0.05:
            note += f";  {r.unmet_floor_GtCO2:.1f} Gt of floor unmet"
        ax.text(tot + 0.4, yi, note, va="center", fontsize=7.5)
    ax.set_yticks(y)
    ax.set_yticklabels([names[g] for g in d.index], color=INK)
    ax.set_xlabel("Gt CO₂ per year (2022 consumption-based CO₂, allocated by the 2019 footprint distribution)")
    ax.set_xlim(0, 50)
    ax.grid(axis="y", visible=False)
    ax.legend(loc="lower right", fontsize=8)
    w = d.loc["World"]
    fig.text(0.02, -0.04, f"Floor = {w.floor_tCO2e_per_cap:.1f} t CO₂e per person (central Aᵉˢˢ at today's T). "
             f"{w.share_of_world_above_floor:.0%} of world CO₂ lies above it.\nGroups of the global population ranked by "
             "individual footprint (Chancel 2022), Singh–Maddala Lorenz fit.", fontsize=7, color=INK2, va="top")
    _save(fig, "fig3_anatomy_of_emissions")


def fig4_frontier(cfg, units, ladder, floor, rq6):
    f = cfg["frontier"]
    s = frontier.sample(ladder, f["ladder_year_main"])
    m = rq6[(rq6.spec == "main") & (rq6.A_ess_level == "central") & rq6.method.str.startswith("quantile")].set_index("reference")
    a, b = m.iloc[0].intercept, m.iloc[0].slope
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ax.scatter(s[A_COL], s.ladder, s=16, color=INK2, alpha=0.45, linewidths=0)
    x = np.logspace(np.log10(s[A_COL].min()), np.log10(s[A_COL].max()), 100)
    ax.plot(x, a + b * np.log(x), color=BLUE, lw=2, label=f"Frontier: quantile regression, τ = {f['tau']}")
    ax.axvline(floor["central"], color=INK2, lw=0.8, ls="-.")
    ax.text(floor["central"] * 1.04, 2.0, "Aᵉˢˢ", fontsize=8, color=INK2)
    lab = s[s.iso_code.isin(f["labelled_countries"]) & (s.iso_code != "USA")]   # the US is labelled as a reference
    ax.scatter(lab[A_COL], lab.ladder, s=24, color=INK, linewidths=0, zorder=4)
    for _, r in lab.iterrows():
        ax.annotate(r.country, (r[A_COL], r.ladder), textcoords="offset points", xytext=(5, -9 if r.iso_code in ("DEU", "JPN", "CHN") else 4), fontsize=7.5)
    for ref, color, dy in (("US", BLUE, 0), ("EU27", ORANGE, 0)):
        r = m.loc[ref]
        ax.annotate("", xy=(r.A_frontier, r.ladder_ref), xytext=(r.A_ref, r.ladder_ref),
                    arrowprops=dict(arrowstyle="->", color=color, lw=1.4))
        ax.plot(r.A_ref, r.ladder_ref, "o", ms=7, color=color, mec=SURFACE, mew=1.2, zorder=5)
        ax.text(r.A_ref * 1.07, r.ladder_ref + (0.12 if ref == "US" else -0.30), f"{ref}: observed Aᵈᵉᵐ = {r.A_dem_observed:.1f}",
                fontsize=7.5, va="center", bbox=dict(facecolor=SURFACE, edgecolor="none", pad=1, alpha=0.85))
    ax.set(xscale="log", xlim=(s[A_COL].min() * 0.8, 2.6e5), xlabel="A, GDP per capita in 2022 (2011 int-$, log scale)",
           ylabel="Life evaluation (Cantril ladder, 2021–2023)")
    ax.legend(loc="upper left", fontsize=8)
    fig.text(0.125, -0.02, f"{len(s)} countries. Arrows: from each reference to the frontier at the same life evaluation. "
             "EU27 = population-weighted mean of members.\nCross-sectional association, not a causal effect; GDP per capita "
             "is a proxy for consumption.", fontsize=7, color=INK2, va="top")
    _save(fig, "fig4_wellbeing_frontier")


def main():
    import floor as fl
    cfg = load.load_config()
    y = cfg["years"]["activity"]
    units = pd.read_csv(load.PROC / f"units_{y}.csv")
    rq2 = pd.read_csv(RES / "rq2_floor.csv")
    floor = fl.floor_values(rq2)
    rq6 = pd.read_csv(RES / "rq6_frontier.csv")
    fig1_feasibility(cfg, units, pd.read_csv(load.PROC / "t_scenarios.csv"), pd.read_csv(load.PROC / "population.csv"),
                     pd.read_csv(load.PROC / "carbon_budget.csv"), floor, rq6)
    fig2_floor(rq2)
    fig3_anatomy(pd.read_csv(RES / "rq3_materialism.csv"))
    fig4_frontier(cfg, units, pd.read_csv(load.PROC / "ladder.csv"), floor, rq6)


if __name__ == "__main__":
    main()
