"""
Momentum *delivery* model -- the coupling between cloud geometry and how much
momentum actually reaches the target.

Why this module exists
----------------------
The v4/v5-plan hit model, P_hit = exp[-d_miss^2/(2 sigma_cloud^2)], scores an
engagement as a Bernoulli hit/miss on whether the target lies inside the cloud.
That is not the physics.  A target of cross-section A_t immersed in a cloud
intercepts only the mass that is geometrically incident on A_t:

    m_delivered = M_launch * F(sigma, d_miss, R_t)

where F is the fraction of the cloud's 2-D areal mass distribution that falls
on the target disc.  For an isotropic Gaussian cloud of per-axis std sigma
whose centroid misses the target centre by d, the exact fraction is a
non-central chi-square CDF (equivalently 1 - Marcum Q_1):

    F = ncx2.cdf( (R_t/sigma)^2 ; df=2, nc=(d/sigma)^2 )

and in the small-target limit R_t << sigma this reduces to the familiar areal
density form  F -> (A_t / (2 pi sigma^2)) * exp(-d^2 / (2 sigma^2)).

Consequences (all DERIVED, and central to this study):
  * A big cloud does NOT buy hit probability for free.  Widening the cloud to
    tolerate a miss distance d costs launched mass as sigma^2.  The mass
    penalty relative to a perfect shot is exactly the cloud-area/target-area
    ratio.
  * For a *deterministic* miss d the optimum is sigma* = d/sqrt(2), giving
        M_launch = m_required * pi * e * d^2 / A_t.
  * "Hit or miss" is the wrong success variable.  A partial delivery is a
    partial orbit shift, so engagements accumulate: the correct metric is the
    distribution of delivered momentum per shot.
"""
import numpy as np
from scipy.stats import ncx2
from scipy.optimize import minimize_scalar


def delivered_fraction(sigma, d_miss, r_target):
    """
    Fraction of cloud mass intercepted by a disc of radius r_target whose
    centre is offset d_miss from the cloud centroid.  Exact (ncx2 / Marcum Q).
    """
    sigma = np.asarray(sigma, float)
    d = np.asarray(d_miss, float)
    rt = np.asarray(r_target, float)
    out = ncx2.cdf((rt / sigma) ** 2, df=2, nc=(d / sigma) ** 2)
    return np.clip(out, 0.0, 1.0)


def delivered_fraction_smalltarget(sigma, d_miss, a_target):
    """Small-target limit, kept for cross-checking the exact expression."""
    return (a_target / (2 * np.pi * sigma ** 2)) * np.exp(-d_miss ** 2 / (2 * sigma ** 2))


def optimal_sigma(d_miss, r_target):
    """sigma minimising launched mass for a deterministic miss distance."""
    if d_miss <= 0:
        return max(r_target, 1e-4)
    lo, hi = 1e-4 * max(d_miss, r_target), 20 * max(d_miss, r_target)
    res = minimize_scalar(
        lambda s: -np.log(max(delivered_fraction(s, d_miss, r_target), 1e-300)),
        bounds=(lo, hi), method="bounded", options={"xatol": 1e-9})
    return float(res.x)


def launch_mass_deterministic(m_required, d_miss, r_target, sigma=None):
    """
    Launched mass needed to deliver m_required to the target given a known
    miss distance.  Returns (M_launch, sigma_used, mass_penalty_factor).
    """
    s = optimal_sigma(d_miss, r_target) if sigma is None else sigma
    f = delivered_fraction(s, d_miss, r_target)
    if f <= 0:
        return np.inf, s, np.inf
    return m_required / f, s, 1.0 / f


def expected_fraction_random_miss(sigma, sigma_miss, r_target):
    """
    Expected delivered fraction when the miss is a 2-D zero-mean Gaussian with
    per-axis std sigma_miss.  Cloud and miss distributions convolve, so the
    expectation is the centred fraction with sigma_eff^2 = sigma^2 + sigma_miss^2.
    """
    s_eff = np.sqrt(np.asarray(sigma, float) ** 2 + np.asarray(sigma_miss, float) ** 2)
    return 1.0 - np.exp(-np.asarray(r_target, float) ** 2 / (2.0 * s_eff ** 2))


def prob_sufficient_delivery(M_launch, m_required, sigma, sigma_miss, r_target,
                             n_samples=200000, rng=None):
    """
    P(delivered mass >= m_required) with a Rayleigh-distributed miss.
    Monte-Carlo over the miss distribution (exact fraction per sample).
    """
    rng = np.random.default_rng(0 if rng is None else rng)
    d = sigma_miss * np.sqrt(rng.chisquare(2, n_samples))
    frac = delivered_fraction(sigma, d, r_target)
    return float(np.mean(M_launch * frac >= m_required)), M_launch * frac


def sigma_for_confidence(m_required, sigma_miss, r_target, M_launch,
                         conf=0.9):
    """
    Given a launched mass budget, the cloud sigma maximising P(sufficient).
    Returns (best_sigma, best_probability).
    """
    grid = np.geomspace(max(r_target, 1e-4), 30 * max(sigma_miss, r_target), 80)
    best, best_p = grid[0], -1.0
    for s in grid:
        # closed form: sufficient iff d <= d_crit where M*F(d_crit)=m_required
        p = _p_sufficient_closed(M_launch, m_required, s, sigma_miss, r_target)
        if p > best_p:
            best, best_p = s, p
    return float(best), float(best_p)


def _p_sufficient_closed(M_launch, m_required, sigma, sigma_miss, r_target):
    """P(delivered >= required) via the critical miss radius (Rayleigh CDF)."""
    f_needed = m_required / M_launch
    if f_needed >= delivered_fraction(sigma, 0.0, r_target):
        return 0.0
    lo, hi = 0.0, 50.0 * max(sigma, sigma_miss, r_target)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if delivered_fraction(sigma, mid, r_target) >= f_needed:
            lo = mid
        else:
            hi = mid
    d_crit = 0.5 * (lo + hi)
    return float(1.0 - np.exp(-d_crit ** 2 / (2.0 * sigma_miss ** 2)))


def min_launch_mass_for_confidence(m_required, sigma_miss, r_target, conf=0.9):
    """
    Smallest launched mass (optimising sigma jointly) that delivers at least
    m_required with probability >= conf.  Returns (M, sigma, penalty).
    """
    lo = m_required
    hi = m_required * 1e12
    for _ in range(200):
        mid = np.sqrt(lo * hi)
        s, p = sigma_for_confidence(m_required, sigma_miss, r_target, mid)
        if p >= conf:
            hi = mid
        else:
            lo = mid
        if hi / lo < 1.0001:
            break
    s, p = sigma_for_confidence(m_required, sigma_miss, r_target, hi)
    return float(hi), float(s), float(hi / m_required)


def penalty_for_confidence(sigma_miss, a_target, conf=0.9):
    """
    Closed form for the launched-mass penalty needed to deliver the required
    mass with probability `conf`, when the miss is a 2-D zero-mean Gaussian of
    per-axis std sigma_miss and the cloud sigma is chosen optimally.

    Small-target limit:  F(sigma,d) = (A_t/2 pi sigma^2) exp(-d^2/2 sigma^2).
    Sufficiency requires d <= d_crit with d_crit^2 = 2 sigma^2 ln(K/sigma^2),
    K = M A_t/(2 pi m_req).  With a Rayleigh miss, P = 1 - exp(-d_crit^2/2 sigma_m^2),
    which is maximised at sigma^2 = K/e, giving

        M/m_req      = -2 pi e sigma_m^2 ln(1-conf) / A_t
        sigma_optimal = sigma_m * sqrt(-ln(1-conf))

    For conf = 0.9 the penalty is 39.3 * sigma_m^2 / A_t.
    """
    return (-2.0 * np.pi * np.e * np.asarray(sigma_miss, float) ** 2
            * np.log(1.0 - conf) / np.asarray(a_target, float))


def sigma_opt_for_confidence(sigma_miss, conf=0.9):
    return np.asarray(sigma_miss, float) * np.sqrt(-np.log(1.0 - conf))


def grain_statistics(M_launch, m_grain, sigma, d_miss, r_target):
    """
    Discrete-grain check: expected number of grains striking the target and the
    Poisson coefficient of variation of the delivered momentum.
    """
    n_total = M_launch / m_grain
    f = delivered_fraction(sigma, d_miss, r_target)
    n_hit = n_total * f
    cv = np.where(n_hit > 0, 1.0 / np.sqrt(np.maximum(n_hit, 1e-12)), np.inf)
    return n_total, n_hit, cv
