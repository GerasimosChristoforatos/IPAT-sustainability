"""Five main figures (figures/) and three supplementary figures (figures/si/).

Reads results/ and data/processed/ only (and re-runs the seeded Monte Carlo for S2).
Every figure: PNG (300 dpi), SVG and PDF at journal print size (183 mm double column).
One colour per concept (config: palette). Titles state the message as a sentence and are
partly computed from the results, so they cannot drift from the numbers.
"""
import math
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, LogNorm

import floor as fl
import frontier
import load
import model
import montecarlo
import scenarios

ROOT = Path(__file__).resolve().parents[1]
RES, FIG = ROOT / "results", ROOT / "figures"
SI = FIG / "si"
A_COL, T_COL = "A_int2011usd_per_cap", "T_kgCO2_per_int2011usd"
MM = 1 / 25.4

CFG = load.load_config()
C = CFG["palette"]
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 7.5, "axes.titlesize": 8, "axes.labelsize": 7.5,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.edgecolor": C["ink_secondary"], "axes.labelcolor": C["ink"], "text.color": C["ink"],
    "xtick.color": C["ink_secondary"], "ytick.color": C["ink_secondary"],
    "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.6,
    "axes.facecolor": C["surface"], "figure.facecolor": C["surface"], "savefig.facecolor": C["surface"],
    "axes.grid": False, "legend.frameon": False, "svg.fonttype": "none", "pdf.fonttype": 42,
    "hatch.linewidth": 0.8,
})
W2 = CFG["figure_size_mm"]["double"] * MM


# ----------------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------------
def _title(fig, text, width=92):
    fig.text(0.01, 0.995, "\n".join(textwrap.wrap(text, width)), fontsize=9, fontweight="bold",
             va="top", ha="left")


def _save(fig, name, folder=FIG):
    folder.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "svg", "pdf"):
        fig.savefig(folder / f"{name}.{ext}", dpi=300, bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)


def _pct(r):
    return f"{100 * r:.1f}".rstrip("0").rstrip(".") + " %/yr"


def _load():
    y = CFG["years"]["activity"]
    d = dict(
        units=pd.read_csv(load.PROC / f"units_{y}.csv"),
        pop=pd.read_csv(load.PROC / "population.csv"),
        budget=pd.read_csv(load.PROC / "carbon_budget.csv"),
        ladder=pd.read_csv(load.PROC / "ladder.csv"),
        rq1=pd.read_csv(RES / "rq1_burn_rate.csv"),
        rates=pd.read_csv(RES / "rq1_required_rate.csv"),
        rq2=pd.read_csv(RES / "rq2_floor.csv"),
        rq3=pd.read_csv(RES / "rq3_materialism.csv"),
        rq4=pd.read_csv(RES / "rq4_cost_of_existing.csv"),
        rq5=pd.read_csv(RES / "rq5_required_dematerialism.csv"),
        rq6=pd.read_csv(RES / "rq6_frontier.csv"),
    )
    d["floor"] = fl.floor_values(d["rq2"])
    d["B"] = scenarios.budgets(CFG, d["budget"])
    d["T0"] = scenarios.t0_path(CFG, d["units"])
    d["P"] = float(d["pop"].set_index("label").loc["now", "P_persons"])
    d["obs"] = d["rq6"][(d["rq6"].spec == "main") & (d["rq6"].A_ess_level == "central")
                       & d["rq6"].method.str.startswith("quantile")].set_index("reference")
    return d


TARGET = 1.7
H = CFG["carbon_budget"]["horizon_years"]
RATES = {k: v["rate"] for k, v in CFG["t_scenarios"].items()}
SHORT = {"BAU": "Historical", "Ambitious": "Ambitious", "Extreme": "Extreme"}


def _scenario_axis(ax, y_text=None):
    """Mark the three T scenarios on a secondary top axis (horizontal labels)."""
    sec = ax.secondary_xaxis("top")
    sec.set_xticks([100 * r for r in RATES.values()])
    sec.set_xticklabels([f"{SHORT[k]}\n{100 * r:g} %" for k, r in RATES.items()], fontsize=7)
    sec.tick_params(length=3, color=C["ink_secondary"])
    sec.spines["top"].set_visible(False)
    return sec


# ----------------------------------------------------------------------------
# Fig. 1  The gap
# ----------------------------------------------------------------------------
def fig1_gap(d):
    w = d["rq3"][(d["rq3"].A_ess_level == "central") & (d["rq3"].method == "Lorenz fit: singh_maddala")
                 & (d["rq3"].group == "World")].iloc[0]
    I_pl = model.allowance(d["B"][TARGET], H) / model.KG_PER_GT
    u = d["units"].set_index("unit").loc["World"]
    cost = model.cost_of_existing(d["P"], d["floor"]["central"], u[T_COL]) / model.KG_PER_GT
    today = u.I_kgCO2 / model.KG_PER_GT
    met, above, unmet = w.floor_met_GtCO2, w.above_floor_GtCO2, w.unmet_floor_GtCO2

    fig, ax = plt.subplots(figsize=(W2, 62 * MM))
    y = [1, 0]
    ax.barh(1, met, height=0.5, color=C["floor"], edgecolor=C["surface"], lw=1)
    ax.barh(1, above, left=met, height=0.5, color=C["materialism"], edgecolor=C["surface"], lw=1)
    ax.barh(0, met, height=0.5, color=C["floor"], edgecolor=C["surface"], lw=1)
    ax.barh(0, cost - met, left=met, height=0.5, facecolor=C["surface"], edgecolor=C["floor"], hatch="////", lw=0.8)
    ax.axvline(I_pl, color=C["allowance"], lw=1.6, zorder=5)
    ax.text(I_pl, 1.55, f"Allowance Iᵖˡ = {I_pl:.0f} Gt/yr", ha="center", va="bottom", fontsize=7.5, fontweight="bold")
    # segment labels above / below the bars, never on the allowance line
    ax.text(met / 2 + I_pl * 0.0, 1.31, f"floor  {met:.0f}", ha="center", va="bottom", fontsize=7, color=C["floor"])
    ax.text(met + above / 2, 1.31, f"above the floor (materialism)  {above:.0f}", ha="center", va="bottom",
            fontsize=7, color=C["materialism"])
    ax.text(met + (cost - met) / 2 + 0.6, -0.31, f"floor not yet met  {unmet:.0f}", ha="left", va="top", fontsize=7,
            color=C["floor"])
    ax.text(today + 0.6, 1, f"{today:.0f} Gt/yr  ({today / I_pl:.1f} × Iᵖˡ)", va="center", fontsize=7.5)
    ax.text(cost + 0.6, 0, f"{cost:.0f} Gt/yr  ({cost / I_pl:.1f} × Iᵖˡ)", va="center", fontsize=7.5)
    ax.set_yticks(y)
    ax.set_yticklabels(["Today's emissions\n(2022)", "Cost of existing\n(everyone at Aᵉˢˢ, today's T)"], color=C["ink"])
    ax.set_xlim(0, 52)
    ax.set_ylim(-0.6, 1.8)
    ax.set_xlabel("Gt CO₂ per year")
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    _title(fig, "Even a decent life for all would exceed the allowance at today's technology")
    fig.subplots_adjust(top=0.84)
    _save(fig, "fig1_the_gap")


# ----------------------------------------------------------------------------
# Fig. 2  Where materialism sits
# ----------------------------------------------------------------------------
def fig2_materialism(d):
    r = d["rq3"][(d["rq3"].A_ess_level == "central") & (d["rq3"].method == "Lorenz fit: singh_maddala")].set_index("group")
    order = ["World", "top1", "next9", "middle40", "next30", "bottom20"]
    names = {"World": "World", "top1": "Top 1 %", "next9": "Next 9 %", "middle40": "Middle 40 %",
             "next30": "Next 30 %", "bottom20": "Bottom 20 %"}
    r = r.loc[order[::-1]]
    y = np.arange(len(r)) + np.array([0, 0, 0, 0, 0, 0.6])
    top10 = r.loc[["top1", "next9"], "above_floor_GtCO2"].sum()
    w = r.loc["World"]
    fig, ax = plt.subplots(figsize=(W2, 80 * MM))
    ax.barh(y, r.floor_met_GtCO2, height=0.55, color=C["floor"], edgecolor=C["surface"], lw=1, label="Floor (met)")
    ax.barh(y, r.above_floor_GtCO2, left=r.floor_met_GtCO2, height=0.55, color=C["materialism"],
            edgecolor=C["surface"], lw=1, label="Above the floor (materialism)")
    ax.barh(y, r.unmet_floor_GtCO2, left=r.floor_met_GtCO2 + r.above_floor_GtCO2, height=0.55, facecolor=C["surface"],
            edgecolor=C["floor"], hatch="////", lw=0.8, label="Floor not yet met (deficit)")
    for yi, (_, g) in zip(y, r.iterrows()):
        end = g.floor_met_GtCO2 + g.above_floor_GtCO2 + g.unmet_floor_GtCO2
        txt = f"{g.above_floor_GtCO2:.1f} above" if g.above_floor_GtCO2 >= 0.05 else ""
        if g.unmet_floor_GtCO2 >= 0.05:
            txt += ("; " if txt else "") + f"{g.unmet_floor_GtCO2:.1f} deficit"
        ax.text(end + 0.5, yi, txt, va="center", fontsize=7)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{names[g]}" for g in r.index], color=C["ink"])
    ax.set_xlim(0, 48)
    ax.set_xlabel("Gt CO₂ per year (2022)")
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=3)
    _title(fig, f"{w.share_of_world_above_floor:.0%} of today's CO₂ is consumption above the floor; "
                f"the richest 10 % account for {top10 / w.above_floor_GtCO2:.0%} of it")
    fig.subplots_adjust(top=0.86)
    _save(fig, "fig2_where_materialism_sits")


# ----------------------------------------------------------------------------
# Fig. 3  Technology alone
# ----------------------------------------------------------------------------
def fig3_technology_alone(d):
    B = d["B"][TARGET]
    u = d["units"].set_index("unit")
    r = np.linspace(0, CFG["frames"]["path_rate_grid_max"], CFG["frames"]["path_rate_grid_points"])
    req = d["rates"][(d["rates"].P_label == "now") & (d["rates"].target_degC == TARGET)
                     & (d["rates"].case == "technology alone (A_dem = 1)")].set_index("reference")
    fig, ax = plt.subplots(figsize=(W2, 88 * MM))
    ax.axhspan(H, H + 7, color=C["grid"], lw=0)
    ax.text(0.3, H + 3.5, f"Budget not used up within the {H}-year horizon", va="center", fontsize=7, color=C["ink_secondary"])
    for ref, shade in zip(["World", "EU27", "US"], C["technology_shades"]):
        I0 = d["P"] * u.loc[ref, A_COL] * d["T0"]
        yrs = np.array([model.years_to_exhaustion_path(B, I0, x) if x > 0 else B / I0 for x in r])
        rs = req.loc[ref, "required_rate"]
        m = r <= rs
        ax.plot(100 * r[m], yrs[m], color=shade, lw=2, solid_capstyle="round")
        ax.plot(100 * rs, H, "o", ms=6, color=shade, mec=C["surface"], mew=1, zorder=5, clip_on=False)
        lbl = {"US": "US lifestyle", "EU27": "EU27 lifestyle", "World": "World-average lifestyle"}[ref]
        ax.text(100 * rs, H - 2.2, f"{lbl}\nneeds {_pct(rs)}", ha="center", va="top", fontsize=7)
        ax.text(-0.3, yrs[0], f"{yrs[0]:.0f} yr", ha="right", va="center", fontsize=7, color=C["ink_secondary"])
    for k, x in RATES.items():
        ax.axvline(100 * x, color=C["ink_secondary"], lw=0.6, ls=(0, (2, 2)), zorder=0)
    _scenario_axis(ax)
    ax.set_xlim(0, 100 * CFG["frames"]["path_rate_grid_max"])
    ax.set_ylim(0, H + 7)
    ax.set_yticks(range(0, H + 1, 10))
    ax.set_xlabel("Annual decline of T, carbon intensity of consumption (% per year, sustained)")
    ax.set_ylabel(f"Years until the {TARGET} °C budget is used up")
    us = req.loc["US", "required_rate"]
    _title(fig, f"Technology alone: a US lifestyle for all needs carbon intensity to fall {_pct(us).replace(' %/yr', ' %')} "
                f"a year, about {us / RATES['BAU']:.0f} times the historical rate", width=88)
    fig.subplots_adjust(top=0.8)
    _save(fig, "fig3_technology_alone")


# ----------------------------------------------------------------------------
# Fig. 4  Dematerialism already exists
# ----------------------------------------------------------------------------
def fig4_frontier(d):
    f = CFG["frontier"]
    s = frontier.sample(d["ladder"], f["ladder_year_main"])
    o = d["obs"]
    a, b = o.iloc[0].intercept, o.iloc[0].slope
    fig, ax = plt.subplots(figsize=(W2, 100 * MM))
    ax.scatter(s[A_COL], s.ladder, s=10, color="#9a9994", linewidths=0, zorder=2)
    x = np.logspace(np.log10(s[A_COL].min()), np.log10(s[A_COL].max()), 100)
    ax.plot(x, a + b * np.log(x), color=C["dematerialism"], lw=1.8, zorder=3)
    ax.text(1.25e3, a + b * np.log(1.25e3) + 0.32, f"Frontier (τ = {f['tau']})", color=C["dematerialism"], fontsize=7.5,
            rotation=0)
    ax.axvline(d["floor"]["central"], color=C["floor"], lw=1, ls=(0, (4, 2)), zorder=1)
    ax.text(d["floor"]["central"] * 0.95, 1.75, "Aᵉˢˢ", color=C["floor"], fontsize=7.5, ha="right")
    # labelled countries: (iso, dx pts, dy pts, ha)
    lab = {"CRI": (-6, 6, "right"), "MEX": (-6, -10, "right"), "FIN": (6, 3, "left"), "GTM": (-6, 4, "right"),
           "IND": (6, -2, "left"), "NOR": (6, -9, "left"), "CHN": (6, -8, "left")}
    for iso, (dx, dy, ha) in lab.items():
        r = s[s.iso_code == iso].iloc[0]
        ax.scatter(r[A_COL], r.ladder, s=16, color=C["ink"], linewidths=0, zorder=4)
        ax.annotate(r.country, (r[A_COL], r.ladder), textcoords="offset points", xytext=(dx, dy), ha=ha, fontsize=7)
    for ref, dy in (("US", 8), ("EU27", -11)):
        r = o.loc[ref]
        ax.annotate("", xy=(r.A_frontier, r.ladder_ref), xytext=(r.A_ref, r.ladder_ref),
                    arrowprops=dict(arrowstyle="-|>", color=C["dematerialism"], lw=1.3, mutation_scale=8,
                                    shrinkA=4, shrinkB=0), zorder=5)
        ax.scatter(r.A_ref, r.ladder_ref, s=34, color=C["ink"], edgecolor=C["surface"], linewidths=1, zorder=6)
        ax.annotate(f"{ref}: observed Aᵈᵉᵐ = {r.A_dem_observed:.1f}", (r.A_ref, r.ladder_ref), textcoords="offset points",
                    xytext=(7, dy), fontsize=7.5, va="center")
    rob = d["rq6"][(d["rq6"].A_ess_level == "central") & d["rq6"].spec.str.startswith("robustness")].set_index(["spec", "reference"])
    note = (f"Sensitivity of observed Aᵈᵉᵐ (US / EU27):\n"
            f"τ = {f['tau_robustness']}: {rob.loc[('robustness: tau', 'US')].A_dem_observed:.1f} / "
            f"{rob.loc[('robustness: tau', 'EU27')].A_dem_observed:.1f}\n"
            f"ladder 2023–2025: {rob.loc[('robustness: ladder window', 'US')].A_dem_observed:.1f} / "
            f"{rob.loc[('robustness: ladder window', 'EU27')].A_dem_observed:.1f}\n"
            f"95 % CI, US: {o.loc['US'].A_dem_observed_ci_low:.1f}–{o.loc['US'].A_dem_observed_ci_high:.1f}")
    ax.text(0.99, 0.03, note, transform=ax.transAxes, ha="right", va="bottom", fontsize=7, color=C["ink_secondary"])
    ax.set(xscale="log", xlim=(s[A_COL].min() * 0.8, 1.3e5), ylim=(1.5, 8.1),
           xlabel="A, GDP per person in 2022 (2011 int-$, log scale)", ylabel="Life evaluation (Cantril ladder, 2021–2023)")
    us = o.loc["US"]
    _title(fig, f"Dematerialism already exists: the US life evaluation is reached with "
                f"{(us.A_frontier - us.A_ess) / (us.A_ref - us.A_ess):.0%} of US above-floor consumption")
    fig.subplots_adjust(top=0.88)
    _save(fig, "fig4_dematerialism_exists")


# ----------------------------------------------------------------------------
# Fig. 5  Both levers together
# ----------------------------------------------------------------------------
CAT_NONE, CAT_FAIL, X_MAX_NUM = 0.6, 520.0, 300.0


def _cat_x(v, status):
    if status == model.FLOOR_FAILS:
        return CAT_FAIL
    if status == model.NONE_NEEDED:
        return CAT_NONE
    return min(v, X_MAX_NUM)


def fig5_both_levers(d):
    B = d["B"][TARGET]
    u = d["units"].set_index("unit")
    A_us, A_ess = u.loc["US", A_COL], d["floor"]["central"]
    r = np.linspace(0, CFG["frames"]["path_rate_grid_max"], CFG["frames"]["path_rate_grid_points"])
    dem = np.logspace(0, np.log10(300), 300)
    R, D = np.meshgrid(r, dem)
    F = montecarlo.path_factor_vec(R, H)
    ratio = d["P"] * (A_ess + (A_us - A_ess) / D) * d["T0"] * F / B
    o = d["obs"].loc["US"]
    req = d["rates"][(d["rates"].P_label == "now") & (d["rates"].target_degC == TARGET)
                     & (d["rates"].reference == "US")].set_index("case")

    fig = plt.figure(figsize=(W2, 185 * MM))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1], width_ratios=[1, 0.035], hspace=0.62, wspace=0.03)
    ax = fig.add_subplot(gs[0, 0])
    cax = fig.add_subplot(gs[0, 1])
    grays = LinearSegmentedColormap.from_list("over", ["#f1f0ec", "#8f8e89"])
    over = np.ma.masked_where(ratio <= 1, ratio)
    pm = ax.pcolormesh(100 * R, D, over, cmap=grays, norm=LogNorm(1, 30), shading="auto", rasterized=True)
    ax.contour(100 * R, D, ratio, levels=[1], colors=C["allowance"], linewidths=2)
    cb = fig.colorbar(pm, cax=cax, extend="max")
    cb.set_label("Cumulative CO₂ ÷ budget\n(over budget only)", fontsize=7)
    cb.set_ticks([1, 2, 5, 10, 30]); cb.set_ticklabels(["1", "2", "5", "10", "30"])
    cb.outline.set_visible(False)
    ax.axhspan(o.A_dem_observed_ci_low, o.A_dem_observed_ci_high, color=C["dematerialism"], alpha=0.16, lw=0)
    ax.axhline(o.A_dem_observed, color=C["dematerialism"], lw=1.4)
    ax.text(24.6, o.A_dem_observed_ci_high * 1.06, f"Observed Aᵈᵉᵐ, US: {o.A_dem_observed:.1f} "
            f"(95 % CI {o.A_dem_observed_ci_low:.1f}–{o.A_dem_observed_ci_high:.1f})", ha="right", va="bottom",
            fontsize=7, color=C["dematerialism"])
    ax.axvline(100 * RATES["BAU"], color=C["technology"], lw=1.8)
    for k in ("Ambitious", "Extreme"):
        ax.axvline(100 * RATES[k], color=C["ink_secondary"], lw=0.6, ls=(0, (2, 2)))
    _scenario_axis(ax)
    pts = [(req.loc["technology alone (A_dem = 1)", "required_rate"], 1.0, "technology alone", (0, 9)),
           (req.loc["both levers (observed A_dem, central)", "required_rate"], o.A_dem_observed, "both levers", (8, -12))]
    for rr, dd, lbl, off in pts:
        ax.plot(100 * rr, dd, "o", ms=6, color=C["ink"], mec=C["surface"], mew=1, zorder=6, clip_on=False)
        ax.annotate(f"{lbl}: {_pct(rr)}", (100 * rr, dd), textcoords="offset points", xytext=off,
                    ha="left" if off[0] > 0 else "center", fontsize=7, fontweight="bold")
    ax.text(3.0, 150, "Over budget", fontsize=8, color=C["ink"])
    ax.text(18.5, 150, "Within budget", fontsize=8, color=C["ink"])
    ax.set(yscale="log", xlim=(0, 100 * CFG["frames"]["path_rate_grid_max"]), ylim=(1, 300),
           xlabel="Annual decline of T (% per year, sustained)", ylabel="Aᵈᵉᵐ, dematerialism (log scale)")
    ax.set_yticks([1, 2, 5, 10, 30, 100, 300]); ax.set_yticklabels(["1 (none)", "2", "5", "10", "30", "100", "300"])
    ax.text(-0.09, 1.13, "a", transform=ax.transAxes, fontsize=10, fontweight="bold")
    ax.set_title(f"US lifestyle for all, {TARGET} °C (50 %), cumulative over {H} years", loc="left", fontsize=7.5, pad=26)

    # panel b: required vs observed, path frame
    bx = fig.add_subplot(gs[1, 0])
    rq = d["rq5"][(d["rq5"].frame == "path") & (d["rq5"].P_label == "now")]
    rows = [(t, s) for t in CFG["carbon_budget"]["targets_degC"] for s in ["Today", "BAU", "Ambitious", "Extreme"]]
    names = {"Today": "No further decline", "BAU": f"Historical ({100 * RATES['BAU']:g} %)",
             "Ambitious": f"Ambitious ({100 * RATES['Ambitious']:g} %)", "Extreme": f"Extreme ({100 * RATES['Extreme']:g} %)"}
    ypos = {}
    yv = 0
    for t, s in rows[::-1]:
        ypos[(t, s)] = yv
        yv += 1 if s != "Today" else 1.6
    bx.axvspan(o.A_dem_observed_ci_low, o.A_dem_observed_ci_high, color=C["dematerialism"], alpha=0.16, lw=0)
    bx.axvline(o.A_dem_observed, color=C["dematerialism"], lw=1.4)
    eu_obs = d["obs"].loc["EU27"].A_dem_observed
    bx.axvline(eu_obs, color=C["dematerialism"], lw=1, ls=(0, (3, 2)))
    for (t, s), yy in ypos.items():
        for ref, mk, fc in (("US", "o", C["ink"]), ("EU27", "D", C["surface"])):
            g = rq[(rq.reference == ref) & (rq.target_degC == t) & (rq.T_scenario == s)].iloc[0]
            xv = _cat_x(g.A_dem_required_central, g.status_central)
            bx.plot(xv, yy + (0.14 if ref == "US" else -0.14), mk, ms=5 if mk == "o" else 4.2, mfc=fc, mec=C["ink"],
                    mew=0.9, zorder=5, clip_on=False)
    for xs in (0.8, 380):
        bx.axvline(xs, color=C["ink_secondary"], lw=0.5)
    bx.set_xscale("log")
    bx.set_xlim(0.45, 800)
    bx.set_xticks([CAT_NONE, 1, 2, 5, 10, 30, 100, 300, CAT_FAIL])
    bx.set_xticklabels(["none\nneeded", "1", "2", "5", "10", "30", "100", "≥300", "floor\nfails"])
    bx.minorticks_off()
    bx.set_yticks(list(ypos.values()))
    bx.set_yticklabels([f"{t} °C · {names[s]}" for (t, s) in ypos], fontsize=7)
    bx.tick_params(axis="y", length=0)
    bx.spines["left"].set_visible(False)
    bx.set_xlabel("Required Aᵈᵉᵐ* (cumulative over 50 years, path frame)")
    bx.text(-0.42, 1.06, "b", transform=bx.transAxes, fontsize=10, fontweight="bold")
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], marker="o", ls="", mfc=C["ink"], mec=C["ink"], ms=5, label="Required, US lifestyle"),
               Line2D([], [], marker="D", ls="", mfc=C["surface"], mec=C["ink"], ms=4.2, label="Required, EU27 lifestyle"),
               Line2D([], [], color=C["dematerialism"], lw=1.4, label=f"Observed, US {o.A_dem_observed:.1f} (band: 95 % CI)"),
               Line2D([], [], color=C["dematerialism"], lw=1, ls=(0, (3, 2)), label=f"Observed, EU27 {eu_obs:.1f}")]
    bx.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.42, -0.2), ncol=2)
    both = req.loc["both levers (observed A_dem, central)", "required_rate"]
    _title(fig, f"Neither lever closes the gap alone; together, the dematerialism already observed and a sustained "
                f"{_pct(both).replace(' %/yr', ' %')} annual decline in carbon intensity keep a US lifestyle within "
                f"the {TARGET} °C budget", width=90)
    fig.subplots_adjust(top=0.9)
    _save(fig, "fig5_both_levers")


# ----------------------------------------------------------------------------
# Supplementary
# ----------------------------------------------------------------------------
def s1_floor(d):
    rq2 = d["rq2"]
    r = rq2[rq2.route != "CENTRAL"].iloc[::-1].reset_index(drop=True)
    c = rq2[rq2.route == "CENTRAL"].iloc[0]
    k = "A_ess{}_int2011usd_per_cap"
    fig, ax = plt.subplots(figsize=(W2, 75 * MM))
    ax.axvspan(c[k.format("_low")], c[k.format("_high")], color=C["floor"], alpha=0.10, lw=0)
    ax.axvline(c[k.format("")], color=C["floor"], lw=1.6)
    for i, row in r.iterrows():
        lo = row[k.format("_low")] if not np.isnan(row[k.format("_low")]) else row[k.format("")]
        hi = row[k.format("_high")] if not np.isnan(row[k.format("_high")]) else row[k.format("")]
        ax.plot([lo, hi], [i, i], color=C["ink_secondary"], lw=1.2, solid_capstyle="round")
        ax.plot(row[k.format("")], i, "o", ms=6, color=C["floor"], mec=C["surface"], mew=1, zorder=3)
    ax.set_yticks(range(len(r)))
    ax.set_yticklabels([f"{a}  {b}" for a, b in zip(r.route, r.label)], color=C["ink"])
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("Aᵉˢˢ (2011 int-$ per person per year)")
    A_w = c[k.format("")] / c.phi
    sec = ax.secondary_xaxis("top", functions=(lambda x: x / A_w, lambda p: p * A_w))
    sec.set_xlabel("φ = floor ÷ current world mean", fontsize=7)
    ax.set_ylim(-0.6, len(r) + 0.2)
    ax.text(c[k.format("")] + 80, len(r) - 0.35, f"central ${c[k.format('')]:,.0f}", fontsize=7, color=C["floor"])
    _title(fig, f"S1. Six published routes put the floor Aᵉˢˢ at ${c[k.format('_low')]:,.0f}–{c[k.format('_high')]:,.0f} "
                f"per person, central ${c[k.format('')]:,.0f}")
    fig.subplots_adjust(top=0.78)
    _save(fig, "S1_floor_triangulation", SI)


def s2_montecarlo(d):
    dr = montecarlo.draws(CFG, d["units"], d["budget"], d["pop"], d["floor"])
    temps = CFG["carbon_budget"]["targets_degC"]
    fig, axes = plt.subplots(2, 3, figsize=(W2, 110 * MM), sharex=True, sharey=True)
    styles = [("BAU", C["technology_shades"][0], "-"), ("Ambitious", C["technology_shades"][1], "-"),
              ("Extreme", C["technology_shades"][2], "-"), ("Pooled (rate uniform)", C["ink_secondary"], (0, (3, 2)))]
    for i, ref in enumerate(CFG["frontier"]["references"]):
        for j, t in enumerate(temps):
            ax = axes[i, j]
            for name, col, ls in styles:
                x = np.sort(dr[("path", ref, t, name)])
                ecdf = np.arange(1, len(x) + 1) / len(x)
                fin = np.isfinite(x)
                ax.plot(np.clip(x[fin], 0.1, 1e3), ecdf[fin], color=col, lw=1.4, ls=ls,
                        label={"BAU": "Historical"}.get(name, name))
            ax.axvline(d["obs"].loc[ref].A_dem_observed, color=C["dematerialism"], lw=1.2)
            ax.set_xscale("log"); ax.set_xlim(0.1, 1e3); ax.set_ylim(0, 1.02)
            ax.grid(color=C["grid"], lw=0.5)
            if i == 0:
                ax.set_title(f"{t} °C", fontsize=7.5)
            if j == 0:
                ax.set_ylabel(f"{ref}\nshare of draws ≤ x")
            if i == 1:
                ax.set_xlabel("Required Aᵈᵉᵐ*")
    axes[1, 1].legend(loc="upper center", bbox_to_anchor=(0.5, -0.38), ncol=4)
    _title(fig, "S2. Monte Carlo distribution of required Aᵈᵉᵐ* (path frame, 10,000 draws): curves that end below 1 "
                "mean the floor fails in the remaining draws; purple line = observed Aᵈᵉᵐ", width=96)
    fig.subplots_adjust(top=0.86, hspace=0.25)
    _save(fig, "S2_montecarlo", SI)


def s3_frames(d):
    rq = d["rq5"][(d["rq5"].P_label == "now")]
    fig, axes = plt.subplots(1, 2, figsize=(W2, 85 * MM), sharey=True)
    rows = [(t, s) for t in CFG["carbon_budget"]["targets_degC"] for s in ["Today", "BAU", "Ambitious", "Extreme"]][::-1]
    names = {"Today": "No further decline", "BAU": "Historical", "Ambitious": "Ambitious", "Extreme": "Extreme"}
    for ax, ref in zip(axes, ["US", "EU27"]):
        for yy, (t, s) in enumerate(rows):
            xs = []
            for frame, mk, fc in (("path", "o", C["ink"]), ("static", "o", C["surface"])):
                g = rq[(rq.reference == ref) & (rq.target_degC == t) & (rq.T_scenario == s) & (rq.frame == frame)].iloc[0]
                xv = _cat_x(g.A_dem_required_central, g.status_central)
                xs.append(xv)
                ax.plot(xv, yy, mk, ms=5, mfc=fc, mec=C["ink"], mew=0.9, zorder=5, clip_on=False)
            ax.plot(xs, [yy, yy], color=C["ink_secondary"], lw=0.6, zorder=1)
        ax.axvline(d["obs"].loc[ref].A_dem_observed, color=C["dematerialism"], lw=1.2)
        for xs_ in (0.8, 380):
            ax.axvline(xs_, color=C["ink_secondary"], lw=0.5)
        ax.set_xscale("log"); ax.set_xlim(0.45, 800); ax.minorticks_off()
        ax.set_xticks([CAT_NONE, 1, 3, 10, 30, 100, 300, CAT_FAIL])
        ax.set_xticklabels(["none\nneeded", "1", "3", "10", "30", "100", "≥300", "floor\nfails"])
        ax.set_title(f"{ref} lifestyle", fontsize=7.5)
        ax.set_xlabel("Required Aᵈᵉᵐ*")
        ax.spines["left"].set_visible(False); ax.tick_params(axis="y", length=0)
    axes[0].set_yticks(range(len(rows)))
    axes[0].set_yticklabels([f"{t} °C · {names[s]}" for t, s in rows])
    from matplotlib.lines import Line2D
    fig.legend(handles=[Line2D([], [], marker="o", ls="", mfc=C["ink"], mec=C["ink"], label="Path frame (main)"),
                        Line2D([], [], marker="o", ls="", mfc=C["surface"], mec=C["ink"], label="Static frame (robustness)"),
                        Line2D([], [], color=C["dematerialism"], label="Observed Aᵈᵉᵐ")],
               loc="lower center", bbox_to_anchor=(0.5, -0.06), ncol=3)
    _title(fig, "S3. The static frame understates required Aᵈᵉᵐ* under fast decline, because it ignores emissions on the way "
                "to the 2050 endpoint", width=96)
    fig.subplots_adjust(top=0.84, bottom=0.2)
    _save(fig, "S3_static_vs_path", SI)


def main():
    d = _load()
    fig1_gap(d)
    fig2_materialism(d)
    fig3_technology_alone(d)
    fig4_frontier(d)
    fig5_both_levers(d)
    s1_floor(d)
    s2_montecarlo(d)
    s3_frames(d)


if __name__ == "__main__":
    main()
