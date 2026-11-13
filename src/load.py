"""Raw inputs -> data/processed. Units are carried in column names.

Outputs
  units_<year>.csv      P, GDP, A, I, T for US, EU27, World (validation units)
  countries_<year>.csv  the same per country (for the T anchor and RQ6)
  carbon_budget.csv     remaining carbon budget by target and likelihood (selected variant)
"""
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
RAW, PROC = ROOT / "data" / "raw", ROOT / "data" / "processed"
KG_PER_MT = 1e9      # kg in one million tonnes
KG_PER_GT = 1e12     # kg in one gigatonne


def load_config() -> dict:
    return yaml.safe_load((ROOT / "config.yaml").read_text())


def load_sources() -> dict:
    return yaml.safe_load((ROOT / "data" / "sources.yaml").read_text())


def source_value(key: str) -> float:
    """Return sources.yaml value for 'source_id.value_name'; refuse unverified inputs."""
    sid, name = key.split(".")
    src = load_sources()[sid]
    if src.get("status") != "verified":
        raise ValueError(f"source '{sid}' is not verified; cannot use '{name}'")
    return float(src["values"][name]["value"])


def _owid(cols):
    return pd.read_csv(RAW / "owid-co2-data.csv", usecols=["country", "iso_code", "year", *cols])


def _derive(df: pd.DataFrame) -> pd.DataFrame:
    """Add A [2011 int-$/person/yr], I [kg CO2/yr], T [kg CO2 / 2011 int-$]."""
    out = df.rename(columns={"population": "P_persons", "gdp": "GDP_int2011usd"}).copy()
    out["I_kgCO2"] = out.pop("consumption_co2") * KG_PER_MT
    out["A_int2011usd_per_cap"] = out["GDP_int2011usd"] / out["P_persons"]
    out["T_kgCO2_per_int2011usd"] = out["I_kgCO2"] / out["GDP_int2011usd"]
    return out


def build_countries(cfg: dict) -> pd.DataFrame:
    y = cfg["years"]["activity"]
    d = _owid(["population", "gdp", "consumption_co2"])
    d = d[(d.year == y) & d.iso_code.notna()].dropna(subset=["population", "gdp", "consumption_co2"])
    return _derive(d).sort_values("iso_code").reset_index(drop=True)


def build_units(cfg: dict) -> pd.DataFrame:
    y = cfg["years"]["activity"]
    d = _owid(["population", "gdp", "consumption_co2"])
    d = d[d.year == y]
    rows = []
    for unit, spec in cfg["units"].items():
        r = d[d.country == spec["owid_country"]].iloc[0].copy()
        if spec.get("gdp_from_members"):
            members = d[d.iso_code.isin(cfg["eu27_iso"])]
            if len(members) != len(cfg["eu27_iso"]) or members.gdp.isna().any():
                raise ValueError("EU27 member GDP incomplete")
            r["gdp"] = members.gdp.sum()
        r["unit"] = unit
        rows.append(r)
    out = _derive(pd.DataFrame(rows))
    return out[["unit", "year", "P_persons", "GDP_int2011usd", "I_kgCO2",
                "A_int2011usd_per_cap", "T_kgCO2_per_int2011usd"]].reset_index(drop=True)


def build_budget(cfg: dict) -> pd.DataFrame:
    src = load_sources()["igcc_2025_rcb"]
    d = pd.read_csv(ROOT / src["file"])
    d = d[d.parameters == cfg["carbon_budget"]["igcc_variant"]].copy()
    d["target_degC"] = d.dT_targets.round(1)
    d = d[d.target_degC.isin(cfg["carbon_budget"]["targets_degC"])]
    d = d.rename(columns={"remaining_carbon_budget": "budget_GtCO2"})
    d["base_year"] = cfg["years"]["budget_base"]
    return d[["target_degC", "likelihood", "budget_GtCO2", "base_year"]].reset_index(drop=True)


def main():
    cfg = load_config()
    PROC.mkdir(parents=True, exist_ok=True)
    y = cfg["years"]["activity"]
    build_units(cfg).to_csv(PROC / f"units_{y}.csv", index=False)
    build_countries(cfg).to_csv(PROC / f"countries_{y}.csv", index=False)
    build_budget(cfg).to_csv(PROC / "carbon_budget.csv", index=False)


if __name__ == "__main__":
    main()
