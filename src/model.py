"""Core identity and closed-form solutions.

    I = P * (A_ess + A_mat / A_dem) * T          constraint  I <= I_pl

Units (SI-consistent throughout this module)
    P      persons
    A_*    2011 int-$ per person per year
    A_dem  dimensionless (1 = no change)
    T      kg CO2 per 2011 int-$
    I      kg CO2 per year
Budgets enter in kg CO2; convert Gt with KG_PER_GT.
"""
import math

KG_PER_GT = 1e12

FLOOR_FAILS = "floor fails"
NONE_NEEDED = "none needed"
REQUIRED = "required"


def impact(P, A_ess, A_mat, A_dem, T):
    """Annual impact I [kg CO2/yr]."""
    return P * (A_ess + A_mat / A_dem) * T


def allowance(budget_kg, horizon_yr):
    """Flat annual allowance I_pl [kg CO2/yr] = remaining budget / horizon."""
    return budget_kg / horizon_yr


def t_endpoint(T0, rate, years):
    """T after `years` of decline at annual `rate` (discrete compounding)."""
    return T0 * (1.0 - rate) ** years


def cost_of_existing(P, A_ess, T):
    """Limit of I as A_dem -> infinity [kg CO2/yr]."""
    return P * A_ess * T


def equal_share_A(I_pl, P, T):
    """Per-person consumption that exactly uses the allowance [int-$/person/yr]."""
    return I_pl / (P * T)


def required_dem(A_ref, A_ess, I_pl, P, T):
    """A_dem* = A_mat / (I_pl/(P*T) - A_ess).

    Returns (value, status). value is math.inf when the floor alone exceeds the
    allowance (status 'floor fails'); status is 'none needed' when value <= 1.
    """
    A_mat = A_ref - A_ess
    headroom = equal_share_A(I_pl, P, T) - A_ess
    if headroom <= 0:
        return math.inf, FLOOR_FAILS
    value = A_mat / headroom
    return value, (NONE_NEEDED if value <= 1 else REQUIRED)


def observed_dem(A_ref, A_frontier, A_ess):
    """Ratio of the reference's A_mat to the A_mat a frontier country needs for the
    same wellbeing. math.inf if the frontier lies at or below the floor."""
    if A_frontier <= A_ess:
        return math.inf
    return (A_ref - A_ess) / (A_frontier - A_ess)


def years_to_exhaustion_static(budget_kg, I):
    """Years until the budget is spent at a constant annual impact I."""
    return budget_kg / I


def decay_constant(rate):
    """Continuous-time decay constant equivalent to a discrete annual rate."""
    return -math.log(1.0 - rate)


def cumulative_path(I_start, rate, years):
    """Cumulative emissions [kg] over `years` when I declines at `rate` from I_start."""
    k = decay_constant(rate)
    return I_start * (1.0 - math.exp(-k * years)) / k


def years_to_exhaustion_path(budget_kg, I_start, rate):
    """Years until the budget is spent when I declines exponentially.
    math.inf if cumulative emissions never reach the budget (I_start/k <= budget)."""
    k = decay_constant(rate)
    x = budget_kg * k / I_start
    return math.inf if x >= 1 else -math.log(1.0 - x) / k
