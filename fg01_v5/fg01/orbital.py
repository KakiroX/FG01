"""
Orbital mechanics helpers for FG01 v5.  All DERIVED.

Conventions
-----------
h        : altitude above mean Earth radius, metres
r        : geocentric radius, metres
Delta-h  : *perigee* decrement produced by a single retrograde apogee burn
"""
import numpy as np
from .constants import MU, R_E


def r_of_h(h):
    return R_E + np.asarray(h, dtype=float)


def v_circ(h):
    """Circular orbital speed at altitude h [m/s]."""
    return np.sqrt(MU / r_of_h(h))


def vis_viva(r, a):
    """Speed on an orbit of semi-major axis a at radius r [m/s]."""
    return np.sqrt(MU * (2.0 / r - 1.0 / a))


def period(a):
    return 2.0 * np.pi * np.sqrt(a ** 3 / MU)


def mean_motion(a):
    return np.sqrt(MU / a ** 3)


def dv_shift(h, dh):
    """
    Delta-v of a single retrograde burn at apogee (= the initial circular
    altitude h) that lowers perigee by dh.

    Post-burn ellipse: apogee r = R_E+h, perigee r = R_E+h-dh,
    so a' = R_E + h - dh/2.
    """
    r = r_of_h(h)
    a_new = r - np.asarray(dh, dtype=float) / 2.0
    return v_circ(h) - vis_viva(r, a_new)


def dv_direct_reentry(h, h_reentry=100.0e3):
    """Delta-v to drop perigee straight to h_reentry (the v4 requirement)."""
    r = r_of_h(h)
    r_p = r_of_h(h_reentry)
    a_new = 0.5 * (r + r_p)
    return v_circ(h) - vis_viva(r, a_new)


def dv_hohmann_two_burn(h1, h2):
    """
    Total two-burn Hohmann delta-v between circular altitudes (validation
    cross-check for dv_shift: the first burn of a Hohmann transfer from h to
    h-dh must equal dv_shift(h, dh)).
    """
    r1, r2 = r_of_h(h1), r_of_h(h2)
    a_t = 0.5 * (r1 + r2)
    dv1 = abs(vis_viva(r1, a_t) - np.sqrt(MU / r1))
    dv2 = abs(np.sqrt(MU / r2) - vis_viva(r2, a_t))
    return dv1, dv2, dv1 + dv2


def elements_from_apo_peri(h_a, h_p):
    """(a, e) from apogee/perigee altitudes."""
    ra, rp = r_of_h(h_a), r_of_h(h_p)
    a = 0.5 * (ra + rp)
    e = (ra - rp) / (ra + rp)
    return a, e


def apo_peri_alt(a, e):
    return a * (1 + e) - R_E, a * (1 - e) - R_E


def state_from_elements(a, e, inc, E, raan=0.0, argp=0.0):
    """
    Inertial position/velocity from (a, e, inc) at eccentric anomaly E.
    RAAN and argument of perigee default to zero (they do not affect the
    orbit-averaged drag rates used here).  Returns r_vec, v_vec in m, m/s.
    """
    cE, sE = np.cos(E), np.sin(E)
    r = a * (1 - e * cE)
    # perifocal position
    x_p = a * (cE - e)
    y_p = a * np.sqrt(1 - e ** 2) * sE
    n = np.sqrt(MU / a ** 3)
    Edot = n / (1 - e * cE)
    vx_p = -a * sE * Edot
    vy_p = a * np.sqrt(1 - e ** 2) * cE * Edot

    ci, si = np.cos(inc), np.sin(inc)
    co, so = np.cos(argp), np.sin(argp)
    cO, sO = np.cos(raan), np.sin(raan)
    # rotation perifocal -> inertial
    R11 = cO * co - sO * so * ci
    R12 = -cO * so - sO * co * ci
    R21 = sO * co + cO * so * ci
    R22 = -sO * so + cO * co * ci
    R31 = so * si
    R32 = co * si
    r_vec = np.array([R11 * x_p + R12 * y_p,
                      R21 * x_p + R22 * y_p,
                      R31 * x_p + R32 * y_p])
    v_vec = np.array([R11 * vx_p + R12 * vy_p,
                      R21 * vx_p + R22 * vy_p,
                      R31 * vx_p + R32 * vy_p])
    return r_vec, v_vec, r


def area_to_mass_sphere(d, rho):
    """A/m for a solid sphere: pi r^2 / ((4/3) pi r^3 rho) = 3/(4 rho r)."""
    d = np.asarray(d, dtype=float)
    return 3.0 / (4.0 * rho * (d / 2.0))


def mass_sphere(d, rho):
    d = np.asarray(d, dtype=float)
    return (4.0 / 3.0) * np.pi * (d / 2.0) ** 3 * rho


def xsec_sphere(d):
    d = np.asarray(d, dtype=float)
    return np.pi * (d / 2.0) ** 2
