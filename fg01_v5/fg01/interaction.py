"""
Resolved cloud-debris interaction physics (INT-1 … INT-5).

Replaces the v5 single-macroscopic-impactor stand-in.  Three pieces:

1. `areal_density` / `mass_efficiency` -- the cloud as a 2-D areal mass-density
   field at the intercept plane, and the fraction of it that lands on the target
   silhouette.

2. `beta_momentum` -- the momentum-enhancement factor derived from crater-ejecta
   scaling instead of assumed.  This is the module's main contribution: v5 had to
   carry a flat U(1.0, 2.5) because no source exists at this scale, and the v6
   plan asks for a physically bounded value for the sub-km/s regime.

3. `ejecta_size_distribution` -- what the many small craters actually throw off,
   for the net-debris check.

--- The beta model ---------------------------------------------------------

beta = 1 + p_ejecta / p_impactor.  A projectile that embeds transfers exactly its
own momentum (beta = 1); ejecta thrown back toward the impactor add to it.

Housen & Holsapple (2011) give the ejecta mass launched faster than speed v, in
the strength regime, as a power law

    M(>v) / m_p = k * (v / v_imp)^(-3 mu)

with mu ~ 0.4-0.55 for competent materials.  Integrating v dM over the ejecta
speed range from the slowest escaping ejecta v_min up to v_imp:

    p_e / p_i = 3 mu k / (1 - 3 mu) * [1 - (v_min/v_imp)^(1-3 mu)]

v_min is set by the target's strength: ejecta slower than ~sqrt(Y/rho) cannot
leave the crater.  For Al-6061 (Y ~ 276 MPa, rho = 2700) that is ~320 m/s --
which is a *large fraction* of a 619 m/s impact speed.  This is exactly why
momentum enhancement is weak here: at sub-km/s the impact barely exceeds the
target's own strength velocity, so almost nothing is ejected fast enough to
matter, and beta sits just above 1.

The same expression at 6 km/s gives the familiar hypervelocity beta ~ 2-3, so
the model reproduces the regime where the literature's numbers come from.  That
is the consistency check, not a claim that those numbers apply here.
"""
import numpy as np
from scipy.stats import ncx2

from .constants import RHO_AL, RHO_RE

# Al-6061 dynamic properties
Y_AL = 276e6            # Pa, yield strength           SOURCED
C_AL = 5100.0           # m/s, bulk sound speed        SOURCED
MU_HH_RANGE = (0.40, 0.55)      # Housen-Holsapple velocity exponent
K_HH_RANGE = (0.2, 0.5)         # ejecta mass coefficient (strength regime)


def v_strength(Y=Y_AL, rho=RHO_AL):
    """Characteristic ejection-speed floor set by target strength [m/s]."""
    return np.sqrt(Y / rho)


def beta_momentum(v_imp, mu=0.45, k=0.3, Y=Y_AL, rho_t=RHO_AL,
                  normal_incidence=True, cos_inc=1.0):
    """
    Momentum-enhancement factor from crater-ejecta scaling.

    Returns beta = 1 + p_ejecta/p_impactor, floored at 1.0 (an impact that
    ejects nothing still delivers its own momentum).

    `cos_inc` applies the oblique-incidence reduction: only the normal component
    of the impact drives cratering, and ejecta momentum scales roughly with it.
    """
    v = np.asarray(v_imp, float)
    v_min = v_strength(Y, rho_t)
    x = np.clip(v_min / np.maximum(v, 1e-9), 0.0, 1.0)
    p3 = 1.0 - 3.0 * mu
    # integral of v dM from v_min to v_imp
    ratio = np.where(
        x >= 1.0, 0.0,
        (3.0 * mu * k / p3) * (1.0 - x ** p3))
    ratio = np.maximum(ratio, 0.0)
    if not normal_incidence:
        ratio = ratio * np.asarray(cos_inc, float)
    return 1.0 + ratio


def beta_envelope(v_imp, Y=Y_AL, rho_t=RHO_AL):
    """(beta_low, beta_central, beta_high) over the mu and k ranges."""
    vals = []
    for mu in MU_HH_RANGE:
        for k in K_HH_RANGE:
            vals.append(beta_momentum(v_imp, mu=mu, k=k, Y=Y, rho_t=rho_t))
    vals = np.array(vals)
    return vals.min(axis=0), beta_momentum(v_imp, 0.45, 0.3, Y, rho_t), vals.max(axis=0)


def penetration_depth(d_g, rho_p, v_imp, Y=Y_AL, rho_t=RHO_AL):
    """
    Crater/penetration depth for a dense sphere into a ductile metal, from the
    standard Poncelet-type resistance law:

        (1/2) rho_p L v^2 = (Y + rho_t v_avg^2) * P     (dimensional balance)

    Implemented in the widely used Poncelet form
        P = (rho_p L / (2 A)) ln(1 + rho_t v^2 / (A Y))  with A ~ 1.5
    which reduces to the hydrodynamic limit P/L -> sqrt(rho_p/rho_t) at high v.
    """
    L = np.asarray(d_g, float)              # sphere: characteristic length ~ d
    v = np.asarray(v_imp, float)
    A = 1.5
    return (rho_p * L / (2 * A * rho_t)) * np.log1p(rho_t * v ** 2 / (A * Y))


def embeds(d_g, rho_p, v_imp, target_thickness, **kw):
    """True if the grain stops inside the target rather than perforating."""
    return penetration_depth(d_g, rho_p, v_imp, **kw) < np.asarray(target_thickness, float)


def crater_volume(m_p, v_imp, Y=Y_AL):
    """Crater volume from the energy/strength balance:  V = KE / Y_dynamic."""
    return 0.5 * np.asarray(m_p, float) * np.asarray(v_imp, float) ** 2 / Y


def ejecta_mass(m_p, v_imp, Y=Y_AL, rho_t=RHO_AL, retention=0.0):
    """
    Mass excavated from the crater [kg].  `retention` is the fraction that
    stays attached (falls back / remains as crater lip); in vacuum with no
    gravity the rest leaves the target.
    """
    return (1.0 - retention) * crater_volume(m_p, v_imp, Y) * rho_t


def ejecta_size_distribution(m_ejecta, l_max, l_min=1e-5, slope=2.7):
    """
    Cumulative number of ejecta fragments larger than l, from a standard
    collisional/cratering power law N(>l) ∝ l^-slope, normalised so the total
    fragment mass equals m_ejecta.

    Returns a callable N(l).  Largest fragment l_max is set by the crater size.
    """
    m_ejecta = float(m_ejecta)
    s = float(slope)

    def mass_integral(lo, hi):
        # dN/dl ∝ l^-(s+1); mass of a fragment ∝ l^3
        if s == 3:
            return np.log(hi / lo)
        return (hi ** (3 - s) - lo ** (3 - s)) / (3 - s)

    norm = m_ejecta / (RHO_AL * (np.pi / 6) * mass_integral(l_min, l_max) * s)

    def N(l):
        l = np.asarray(l, float)
        return np.where(l >= l_max, 0.0, norm * s * (l ** -s - l_max ** -s) / s)

    return N


def areal_density(m_launch, sigma_cloud, x=0.0, y=0.0):
    """Areal MASS density of a Gaussian cloud at offset (x, y) [kg/m^2]."""
    s = np.asarray(sigma_cloud, float)
    r2 = np.asarray(x, float) ** 2 + np.asarray(y, float) ** 2
    return (np.asarray(m_launch, float) / (2 * np.pi * s ** 2)) * np.exp(-r2 / (2 * s ** 2))


def mass_efficiency(sigma_cloud, r_target, d_miss=0.0):
    """
    Fraction of launched cloud mass intercepted by a disc of radius r_target
    offset by d_miss from the cloud centroid.  Exact (non-central chi-square).
    """
    s = np.asarray(sigma_cloud, float)
    return np.clip(ncx2.cdf((np.asarray(r_target, float) / s) ** 2, df=2,
                            nc=(np.asarray(d_miss, float) / s) ** 2), 0.0, 1.0)


def sigma_from_cloud_diameter(d_cloud):
    """
    Convert a quoted cloud 'diameter' to a Gaussian sigma.  Convention: the
    quoted diameter contains 95% of the mass, i.e. d_cloud = 2 * 2.448 * sigma
    (2-D Rayleigh 95th percentile).
    """
    return np.asarray(d_cloud, float) / (2 * 2.4477)


def grains_on_target(m_launch, m_grain, sigma_cloud, r_target, d_miss=0.0):
    """Expected number of grains striking the target, and the Poisson CV."""
    eps = mass_efficiency(sigma_cloud, r_target, d_miss)
    n = np.asarray(m_launch, float) * eps / np.asarray(m_grain, float)
    return n, eps, np.where(n > 0, 1.0 / np.sqrt(np.maximum(n, 1e-12)), np.inf)


def grain_mass(d_g, rho=RHO_RE):
    return (4.0 / 3.0) * np.pi * (np.asarray(d_g, float) / 2.0) ** 3 * rho


def mean_incidence_cos_sphere():
    """
    Mean cos(incidence) for a parallel beam striking a sphere, weighted by the
    projected area:  <cos i> = 2/3.  Grains hit a curved target obliquely, so
    the effective normal-velocity component is reduced by this factor.
    """
    return 2.0 / 3.0
