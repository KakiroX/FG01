"""
Terminal guidance error budget (SIM-4).

Three miss-distance models are carried side by side throughout, because the
choice between them is the single largest driver of the v4 result and must be
visible rather than buried:

M1  "v4 model" (bounding, pessimistic)
        d_miss = sqrt( sigma_pos^2 + (v_rel * tau)^2 )
    The lag term assumes the target sweeps across the aim point at the *impact*
    velocity v_rel and that no extrapolation of the target state is performed.

M2  "co-orbital, no extrapolation"
        the sweep rate is the FG01-to-target relative speed u, not v_rel.
    v_rel is supplied almost entirely by the launcher muzzle velocity; the
    target's motion *relative to the launcher* proceeds at u.  Adds launcher
    pointing error and muzzle-speed dispersion.

M3  "co-orbital + Keplerian extrapolation" (the design model)
        the target state is estimated over a tracking arc and propagated to the
        intercept epoch, so lag enters only through estimation error and
        unmodelled relative acceleration, both computed here.

M2 and M3 are not assumptions of better performance: they follow from the
engagement geometry that the energy gate (SIM-0) already forces, namely a low
relative-velocity, co-planar engagement.  That geometry has its own cost --
FG01 must reach each target's orbit plane -- which is charged in SIM-11.
"""
import numpy as np

from .constants import MU


def sigma_pos_diffraction(R, lam, D, snr=None):
    """
    Angular position uncertainty times range.  The diffraction limit lam/D is
    used as the conservative default; centroiding improves on it roughly as
    (lam/D)/SNR, which is reported when snr is given.
    """
    theta = lam / D
    if snr is not None:
        theta = theta / max(snr, 1.0)
    return theta * np.asarray(R, float)


def sigma_velocity_estimate(sigma_pos, t_obs, n_samples):
    """
    1-sigma velocity error from a least-squares fit of position over a tracking
    arc of length t_obs with n_samples uniformly spaced measurements:
        sigma_v = sigma_pos * sqrt(12 / n) / t_obs      (large-n limit)
    """
    n = max(float(n_samples), 2.0)
    return np.asarray(sigma_pos, float) * np.sqrt(12.0 / n) / float(t_obs)


def relative_accel_scale(a_orbit, sep):
    """
    Magnitude of the Keplerian relative (tidal) acceleration between two nearby
    co-orbital objects separated by `sep`:  ~ n^2 * sep, n = sqrt(mu/a^3).
    This is what an extrapolating guidance law fails to model if it assumes
    straight-line relative motion.
    """
    n = np.sqrt(MU / np.asarray(a_orbit, float) ** 3)
    return 3.0 * n ** 2 * np.asarray(sep, float)   # worst-case (radial) coefficient


def miss_v4(R, tau, v_rel, lam, D):
    """M1: the v4 model."""
    sp = sigma_pos_diffraction(R, lam, D)
    return np.sqrt(sp ** 2 + (np.asarray(v_rel, float) * np.asarray(tau, float)) ** 2)


def miss_coorbital(R, tau, u, v_rel, lam, D, sigma_point, dv_muzzle_frac):
    """
    M2: sweep rate u, plus launcher pointing and muzzle-speed dispersion.

    Muzzle-speed dispersion delta_v/v produces an arrival-time error
    dt = R*delta_v/v^2; during dt the target drifts u*dt relative to the aim
    point, hence the term u * R * (delta_v/v) / v.
    """
    sp = sigma_pos_diffraction(R, lam, D)
    lag = np.asarray(u, float) * np.asarray(tau, float)
    point = np.asarray(sigma_point, float) * np.asarray(R, float)
    timing = (np.asarray(u, float) * np.asarray(R, float)
              * np.asarray(dv_muzzle_frac, float) / np.asarray(v_rel, float))
    return np.sqrt(sp ** 2 + lag ** 2 + point ** 2 + timing ** 2), \
        dict(sigma_pos=sp, lag=lag, pointing=point, timing=timing)


def miss_extrapolating(R, tau, u, v_rel, lam, D, sigma_point, dv_muzzle_frac,
                       t_obs, n_samples, a_orbit, snr=None):
    """
    M3: full error budget with Keplerian extrapolation of the target state.

    t_go = tau + R/v_rel is the interval between the last usable measurement
    and impact.  Contributions:
        position knowledge          sigma_pos
        velocity knowledge          sigma_v * t_go
        unmodelled rel. accel.      0.5 * a_rel * t_go^2
        launcher pointing           sigma_point * R
        muzzle-speed dispersion     u * R * (dv/v) / v
    """
    sp = sigma_pos_diffraction(R, lam, D, snr)
    t_go = np.asarray(tau, float) + np.asarray(R, float) / np.asarray(v_rel, float)
    sv = sigma_velocity_estimate(sp, t_obs, n_samples)
    a_rel = relative_accel_scale(a_orbit, R)
    terms = dict(
        sigma_pos=sp,
        vel_extrap=sv * t_go,
        rel_accel=0.5 * a_rel * t_go ** 2,
        pointing=np.asarray(sigma_point, float) * np.asarray(R, float),
        timing=(np.asarray(u, float) * np.asarray(R, float)
                * np.asarray(dv_muzzle_frac, float) / np.asarray(v_rel, float)),
    )
    total = np.sqrt(sum(v ** 2 for v in terms.values()))
    return total, terms


def p_hit_plan_model(d_miss, sigma_cloud):
    """The plan's Bernoulli hit model, retained for side-by-side reporting."""
    return np.exp(-np.asarray(d_miss, float) ** 2
                  / (2.0 * np.asarray(sigma_cloud, float) ** 2))
