"""
Orbital-decay propagation (SIM-2 engine).

Two independent propagators:

1. `lifetime()` -- orbit-averaged Gauss variational equations in (a, e).
   The drag acceleration is evaluated at Gauss-Legendre nodes in eccentric
   anomaly around the orbit, with the *vector* atmosphere co-rotation term
   v_rel = v - omega_E x r included, and averaged with the correct time weight
   dt = (1 - e cos E)/n dE.  Density comes from the real NRLMSISE-00 table.
   This is the workhorse: it handles the eccentric post-kick orbit that sweeps
   orders of magnitude in density per revolution.

2. `propagate_cartesian()` -- full 3-D Cartesian RK4 with the same drag model
   plus J2, used to validate (1).  Too slow for multi-decade lifetimes, so it
   is run over bounded windows and on fast-decay cases.

Reentry is declared when perigee altitude < 100 km (the conventional
entry-interface altitude).
"""
import numpy as np

from .constants import MU, R_E, OMEGA_E, J2, CD_SPHERE
from .atmosphere import density

YEAR = 365.25 * 86400.0
H_REENTRY = 100.0e3
T_MAX_DEFAULT = 200.0 * YEAR

# Gauss-Legendre nodes on [0, 2pi) for the orbit average.  Perigee-heavy
# density variation is handled by node count, not by weighting tricks.
_NQ = 96
_x, _w = np.polynomial.legendre.leggauss(_NQ)
_E_NODES = np.pi * (_x + 1.0)          # map [-1,1] -> [0, 2pi]
_E_W = np.pi * _w                      # dE weights


def _rates(a, e, inc, aom, cd, case):
    """
    Orbit-averaged (da/dt, de/dt) [m/s, 1/s] for drag with co-rotation.

    a   : semi-major axis [m]
    e   : eccentricity
    inc : inclination [rad]
    aom : A/m [m^2/kg]
    cd  : drag coefficient
    """
    e = max(float(e), 0.0)
    E = _E_NODES
    cE, sE = np.cos(E), np.sin(E)
    one_m_ecE = 1.0 - e * cE
    r = a * one_m_ecE
    alt = r - R_E

    # position/velocity in the orbital plane, then rotated by inclination
    x_p = a * (cE - e)
    y_p = a * np.sqrt(max(1.0 - e * e, 1e-16)) * sE
    n = np.sqrt(MU / a ** 3)
    Edot = n / one_m_ecE
    vx_p = -a * sE * Edot
    vy_p = a * np.sqrt(max(1.0 - e * e, 1e-16)) * cE * Edot

    ci, si = np.cos(inc), np.sin(inc)
    rx, ry, rz = x_p, y_p * ci, y_p * si
    vx, vy, vz = vx_p, vy_p * ci, vy_p * si

    # atmospheric co-rotation: v_atm = omega_E zhat x r
    vax, vay, vaz = -OMEGA_E * ry, OMEGA_E * rx, np.zeros_like(rz)
    wx, wy, wz = vx - vax, vy - vay, vz - vaz
    w = np.sqrt(wx * wx + wy * wy + wz * wz)

    rho = density(alt, case)
    k = -0.5 * rho * cd * aom * w
    ax, ay, az = k * wx, k * wy, k * wz

    # RSW decomposition
    rn = np.sqrt(rx * rx + ry * ry + rz * rz)
    ur = np.stack([rx / rn, ry / rn, rz / rn])
    hx = ry * vz - rz * vy
    hy = rz * vx - rx * vz
    hz = rx * vy - ry * vx
    hn = np.sqrt(hx * hx + hy * hy + hz * hz)
    uw = np.stack([hx / hn, hy / hn, hz / hn])
    us = np.stack([uw[1] * ur[2] - uw[2] * ur[1],
                   uw[2] * ur[0] - uw[0] * ur[2],
                   uw[0] * ur[1] - uw[1] * ur[0]])
    aR = ax * ur[0] + ay * ur[1] + az * ur[2]
    aS = ax * us[0] + ay * us[1] + az * us[2]

    # true anomaly from eccentric anomaly
    sqrt1me2 = np.sqrt(max(1.0 - e * e, 1e-16))
    sin_nu = sqrt1me2 * sE / one_m_ecE
    cos_nu = (cE - e) / one_m_ecE
    p = a * (1.0 - e * e)

    dadt = 2.0 / (n * sqrt1me2) * (e * sin_nu * aR + (p / r) * aS)
    dedt = sqrt1me2 / (n * a) * (sin_nu * aR + (cos_nu + cE) * aS)

    # time-average: (1/2pi) * integral f * (1 - e cosE) dE
    wt = _E_W * one_m_ecE / (2.0 * np.pi)
    return float(np.sum(dadt * wt)), float(np.sum(dedt * wt))


def lifetime(a0, e0, aom, inc_deg=85.0, cd=CD_SPHERE, case="nominal",
             t_max=T_MAX_DEFAULT, record=False):
    """
    Time to reentry (perigee altitude < 100 km) in seconds.

    Returns (t_reentry, history) with history = None unless record=True.
    If the object has not reentered by t_max, returns (np.inf, history).
    """
    inc = np.deg2rad(inc_deg)
    a, e, t = float(a0), float(e0), 0.0
    hist = [] if record else None

    while t < t_max:
        rp = a * (1.0 - e)
        if rp - R_E < H_REENTRY:
            break
        da, de = _rates(a, e, inc, aom, cd, case)
        if record:
            hist.append((t, a, e, rp - R_E, a * (1 + e) - R_E))
        if da == 0.0 or not np.isfinite(da):
            return np.inf, hist
        # adaptive step: limit change in a to min(2 km, 0.2% of a) and in
        # perigee radius to 2 km, and never step more than 1 year at once.
        dt_a = min(2000.0, 0.002 * a) / abs(da)
        drp = abs(da * (1 - e) - a * de)
        dt_p = 2000.0 / drp if drp > 0 else np.inf
        step = min(dt_a, dt_p, 1.0 * YEAR, t_max - t)
        if step <= 0:
            break
        # RK2 (midpoint) on (a, e)
        a_m = a + 0.5 * step * da
        e_m = max(e + 0.5 * step * de, 0.0)
        if a_m * (1 - e_m) - R_E < H_REENTRY:
            step *= 0.5
            a_m = a + 0.5 * step * da
            e_m = max(e + 0.5 * step * de, 0.0)
        da2, de2 = _rates(a_m, e_m, inc, aom, cd, case)
        a += step * da2
        e = max(e + step * de2, 0.0)
        t += step
        if a - R_E < H_REENTRY:
            break
    else:
        return np.inf, hist
    return t, hist


def lifetime_years(*args, **kwargs):
    t, _ = lifetime(*args, **kwargs)
    return t / YEAR


def _accel(rv, aom, cd, case, use_j2=True):
    r = rv[:3]
    v = rv[3:]
    rn = np.linalg.norm(r)
    alt = rn - R_E
    # two-body + J2
    a_grav = -MU * r / rn ** 3
    if use_j2:
        z2r2 = (r[2] / rn) ** 2
        f = 1.5 * J2 * MU * R_E ** 2 / rn ** 5
        a_grav = a_grav + np.array([f * r[0] * (5 * z2r2 - 1),
                                    f * r[1] * (5 * z2r2 - 1),
                                    f * r[2] * (5 * z2r2 - 3)])
    # drag with co-rotation
    v_atm = np.array([-OMEGA_E * r[1], OMEGA_E * r[0], 0.0])
    w = v - v_atm
    wn = np.linalg.norm(w)
    rho = float(density(alt, case))
    a_drag = -0.5 * rho * cd * aom * wn * w
    return np.concatenate([v, a_grav + a_drag])


def propagate_cartesian(r0, v0, aom, cd=CD_SPHERE, case="nominal",
                        t_end=30 * 86400.0, dt=20.0, stop_alt=H_REENTRY,
                        use_j2=True):
    """
    Full 3-D Cartesian RK4 propagation (validation propagator).
    Returns (t, a, e, alt_perigee_estimate, reentered_bool).

    `use_j2=False` is used when comparing against the orbit-averaged
    propagator: J2 produces short-period oscillations of the *osculating*
    semi-major axis of order J2*a*(R_E/p)^2 (several km at LEO) which are not
    secular and would otherwise swamp the drag signal being compared.
    """
    rv = np.concatenate([np.asarray(r0, float), np.asarray(v0, float)])
    t = 0.0
    while t < t_end:
        step = min(dt, t_end - t)
        k1 = _accel(rv, aom, cd, case, use_j2)
        k2 = _accel(rv + 0.5 * step * k1, aom, cd, case, use_j2)
        k3 = _accel(rv + 0.5 * step * k2, aom, cd, case, use_j2)
        k4 = _accel(rv + step * k3, aom, cd, case, use_j2)
        rv = rv + (step / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        t += step
        rn = np.linalg.norm(rv[:3])
        if rn - R_E < stop_alt:
            return t, *elements_from_state(rv), True
    return t, *elements_from_state(rv), False


def elements_from_state(rv):
    r, v = rv[:3], rv[3:]
    rn = np.linalg.norm(r)
    vn = np.linalg.norm(v)
    energy = 0.5 * vn ** 2 - MU / rn
    a = -MU / (2 * energy)
    h = np.cross(r, v)
    e_vec = np.cross(v, h) / MU - r / rn
    e = np.linalg.norm(e_vec)
    return a, e, a * (1 - e) - R_E
