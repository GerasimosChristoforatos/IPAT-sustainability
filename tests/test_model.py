"""Validation identity, closed-form solutions and data-consistency checks."""
import math
import sys
import unittest
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import load   # noqa: E402
import model  # noqa: E402

CFG = load.load_config()


class TestIdentity(unittest.TestCase):
    """P x A x T reproduces observed consumption CO2. T is defined as I/(P*A), so this
    is an identity: it checks units, year alignment and EU27 aggregation, not the model."""

    def test_units_reproduce_observed(self):
        y = CFG["years"]["activity"]
        units = pd.read_csv(ROOT / "data" / "processed" / f"units_{y}.csv")
        raw = pd.read_csv(ROOT / "data" / "raw" / "owid-co2-data.csv",
                          usecols=["country", "year", "consumption_co2"])
        for _, u in units.iterrows():
            name = CFG["units"][u.unit]["owid_country"]
            obs_mt = raw[(raw.country == name) & (raw.year == y)].consumption_co2.iloc[0]
            got_mt = u.P_persons * u.A_int2011usd_per_cap * u.T_kgCO2_per_int2011usd / load.KG_PER_MT
            self.assertAlmostEqual(got_mt, obs_mt, delta=0.5, msg=u.unit)   # within 0.5 Mt

    def test_eu27_members_match_aggregate(self):
        y = CFG["years"]["activity"]
        raw = pd.read_csv(ROOT / "data" / "raw" / "owid-co2-data.csv",
                          usecols=["country", "iso_code", "year", "population", "consumption_co2"])
        raw = raw[raw.year == y]
        mem = raw[raw.iso_code.isin(CFG["eu27_iso"])]
        agg = raw[raw.country == "European Union (27)"].iloc[0]
        self.assertEqual(len(mem), 27)
        self.assertAlmostEqual(mem.population.sum() / agg.population, 1, places=4)
        self.assertAlmostEqual(mem.consumption_co2.sum() / agg.consumption_co2, 1, places=4)

    def test_budget_variant_matches_paper_headline(self):
        b = pd.read_csv(ROOT / "data" / "processed" / "carbon_budget.csv")
        v = b[(b.target_degC == 1.5) & (b.likelihood == 0.5)].budget_GtCO2.iloc[0]
        self.assertEqual(round(v, -1), 130)   # Forster et al. 2026 headline, from start of 2026


class TestClosedForms(unittest.TestCase):
    P, A_ref, A_ess, T = 8.0e9, 50_000.0, 5_000.0, 0.1

    def test_no_dematerialism_is_plain_ipat(self):
        I = model.impact(self.P, self.A_ess, self.A_ref - self.A_ess, 1.0, self.T)
        self.assertAlmostEqual(I, self.P * self.A_ref * self.T)

    def test_limit_is_cost_of_existing(self):
        I = model.impact(self.P, self.A_ess, self.A_ref - self.A_ess, 1e15, self.T)
        self.assertAlmostEqual(I / model.cost_of_existing(self.P, self.A_ess, self.T), 1, places=9)

    def test_required_dem_hits_allowance(self):
        I_pl = 10 * model.KG_PER_GT
        d, status = model.required_dem(self.A_ref, self.A_ess, I_pl, self.P, self.T)
        self.assertEqual(status, model.REQUIRED)
        I = model.impact(self.P, self.A_ess, self.A_ref - self.A_ess, d, self.T)
        self.assertAlmostEqual(I / I_pl, 1, places=12)

    def test_floor_fails(self):
        I_pl = 0.99 * model.cost_of_existing(self.P, self.A_ess, self.T)
        d, status = model.required_dem(self.A_ref, self.A_ess, I_pl, self.P, self.T)
        self.assertEqual((d, status), (math.inf, model.FLOOR_FAILS))

    def test_none_needed(self):
        I_pl = 1.01 * self.P * self.A_ref * self.T
        d, status = model.required_dem(self.A_ref, self.A_ess, I_pl, self.P, self.T)
        self.assertEqual(status, model.NONE_NEEDED)
        self.assertLessEqual(d, 1)

    def test_equal_share(self):
        I_pl = 10 * model.KG_PER_GT
        A = model.equal_share_A(I_pl, self.P, self.T)
        self.assertAlmostEqual(self.P * A * self.T, I_pl)

    def test_observed_dem(self):
        self.assertAlmostEqual(model.observed_dem(50_000, 20_000, 5_000), 3.0)
        self.assertEqual(model.observed_dem(50_000, 4_000, 5_000), math.inf)

    def test_t_endpoint(self):
        self.assertAlmostEqual(model.t_endpoint(1.0, 0.069, 28), 0.931 ** 28)

    def test_path_matches_numerical_sum(self):
        I0, r = 40.0, 0.069
        # continuous integral must lie between the left and right annual Riemann sums
        left = sum(I0 * (1 - r) ** t for t in range(0, 30))
        right = sum(I0 * (1 - r) ** t for t in range(1, 31))
        self.assertTrue(right < model.cumulative_path(I0, r, 30) < left)

    def test_path_exhaustion_inverts_cumulative(self):
        I0, r, B = 40.0, 0.03, 500.0
        t = model.years_to_exhaustion_path(B, I0, r)
        self.assertAlmostEqual(model.cumulative_path(I0, r, t), B, places=8)

    def test_path_never_exhausts(self):
        I0, r = 40.0, 0.115
        bound = I0 / model.decay_constant(r)
        self.assertEqual(model.years_to_exhaustion_path(bound * 1.01, I0, r), math.inf)

    def test_static_exhaustion(self):
        self.assertAlmostEqual(model.years_to_exhaustion_static(500.0, 40.0), 12.5)


class TestMonteCarlo(unittest.TestCase):
    def test_vectorised_matches_scalar(self):
        import numpy as np
        import montecarlo
        A_ess = np.array([4000.0, 6000.0, 9000.0])
        I_pl = np.array([10.0, 3.0, 20.0]) * model.KG_PER_GT
        v = montecarlo.required_dem_vec(50_000.0, A_ess, I_pl, 8.0e9, 0.1)
        for i in range(3):
            self.assertEqual(v[i], model.required_dem(50_000.0, A_ess[i], I_pl[i], 8.0e9, 0.1)[0])


class TestFrontier(unittest.TestCase):
    def test_lp_quantile_regression_matches_statsmodels(self):
        import numpy as np
        import frontier
        lad = pd.read_csv(ROOT / "data" / "processed" / "ladder.csv")
        s = frontier.sample(lad, CFG["frontier"]["ladder_year_main"])
        tau = CFG["frontier"]["tau"]
        lp = frontier.qr_params(np.log(s[frontier.A_COL].values), s.ladder.values, tau)
        irls = frontier.fit_qr(s, tau).params
        # compare the objective (check loss), which is unique even if the minimiser is not
        def loss(b):
            r = s.ladder.values - b[0] - b[1] * np.log(s[frontier.A_COL].values)
            return float(np.sum(np.where(r >= 0, tau * r, (tau - 1) * r)))
        self.assertLessEqual(loss(lp), loss(irls) + 1e-9)
        self.assertAlmostEqual(lp[1], irls[1], delta=0.02)

    def test_frontier_inverse(self):
        import frontier
        self.assertAlmostEqual(frontier.frontier_A(1.0, 0.5, 6.0), math.exp(10.0))
        self.assertEqual(frontier.frontier_A(1.0, -0.1, 6.0), math.inf)


class TestDistribution(unittest.TestCase):
    def test_chancel_groups_reproduce_world_mean(self):
        import floor
        g = floor.chancel_groups(load.load_sources())
        self.assertAlmostEqual(g.pop_share.sum(), 1.0)
        mean = (g.pop_share * g.mean_tCO2e).sum()
        self.assertAlmostEqual(mean, load.load_sources()["chancel_2022"]["values"]["world_mean"]["value"], delta=0.05)

    def test_group_mean_estimate_is_lower_bound(self):
        """Jensen: splitting a group into unequal halves cannot lower above-floor emissions."""
        import numpy as np
        import floor
        f = 2.0
        whole = floor._split(np.array([3.0]), f)[1]
        halves = floor._split(np.array([1.0, 5.0]), f)[1]
        self.assertGreaterEqual(halves, whole)


if __name__ == "__main__":
    unittest.main()
