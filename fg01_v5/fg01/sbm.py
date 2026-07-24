"""
NASA Standard Breakup Model (Johnson et al. 2001) -- fragment production.

AM-13: the effective-mass term is regime dependent and is the single most
error-prone part of the model, so both branches are implemented explicitly and
the branch actually used is returned with every result.

    catastrophic      (E_s >= E_s,c) :  M_eff = m_target + m_projectile
    non-catastrophic  (E_s <  E_s,c) :  M_eff = m_projectile * v_imp[km/s]

    N(>L_c) = 0.1 * M_eff^0.75 * L_c^-1.71        (L_c in metres)

Calibration caveat (reported everywhere it is used): the SBM was fitted to
hypervelocity (>~1-2 km/s) impact and explosion data.  The FG01 v5 design point
sits at 0.5-1.1 km/s, i.e. *below* the calibration range, where the model is an
extrapolation.  Because fragment yield falls steeply with impact speed, using
the SBM here over-predicts fragment counts -- the direction that is
conservative for a debris-generation claim.
"""
import numpy as np


def specific_energy(m_p, v_rel, m_t):
    """E_s = 1/2 m_p v_rel^2 / m_t   [J/kg]."""
    return 0.5 * np.asarray(m_p, float) * np.asarray(v_rel, float) ** 2 / np.asarray(m_t, float)


def regime(m_p, v_rel, m_t, es_c=40.0e3):
    es = specific_energy(m_p, v_rel, m_t)
    return np.where(es >= es_c, "catastrophic", "non-catastrophic"), es


def effective_mass(m_p, v_rel, m_t, es_c=40.0e3):
    """Regime-dependent effective mass [kg] used in the SBM count law."""
    es = specific_energy(m_p, v_rel, m_t)
    cat = es >= es_c
    m_cat = np.asarray(m_t, float) + np.asarray(m_p, float)
    m_non = np.asarray(m_p, float) * (np.asarray(v_rel, float) / 1000.0)
    return np.where(cat, m_cat, m_non), cat, es


def n_fragments(m_p, v_rel, m_t, lc_m, es_c=40.0e3):
    """
    Cumulative number of fragments with characteristic length > lc_m.
    Returns (N, effective_mass, is_catastrophic, E_s).
    """
    m_eff, cat, es = effective_mass(m_p, v_rel, m_t, es_c)
    n = 0.1 * m_eff ** 0.75 * np.asarray(lc_m, float) ** -1.71
    return n, m_eff, cat, es


def crater_ejecta_mass(m_p, v_rel, target_material="Al"):
    """
    Ejecta mass from a sub-catastrophic impact, from the standard cratering
    mass-ratio scaling  M_e/m_p = K * (v/1 km/s)^~1  (Gault-type).

    K is material dependent; for aluminium targets impacted by dense metal at
    0.3-3 km/s, laboratory ejecta-to-projectile mass ratios cluster around
    5-20 per km/s.  Central value K = 10, range 5-20  -- UNVALIDATED at the
    exact combination (Re on Al-6061, sub-km/s), reported as a range.
    """
    k = {"Al": 10.0}.get(target_material, 10.0)
    v_km = np.asarray(v_rel, float) / 1000.0
    return k * np.asarray(m_p, float) * v_km, (5.0 * np.asarray(m_p, float) * v_km,
                                               20.0 * np.asarray(m_p, float) * v_km)
