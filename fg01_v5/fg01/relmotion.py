"""
Engagement geometry and relative motion (INT-0).

This module exists to test the single most attackable claim in the project:
that the target's transverse sweep across the aim point is the *platform-target
relative speed* u, not the *impact* speed v_rel.

The decisive quantity is not the magnitude of the relative velocity but its
DIRECTION relative to the target's own velocity vector, because only the
retrograde component lowers the orbit.

    w   = v_projectile - v_target          (impact relative velocity)
    Ds  = beta * m_p * w / m_t             (velocity change of the target)
    useful (retrograde) part = -Ds . u_t

Two ways to obtain a given |w|:

  (a) CO-ORBITAL: platform matches the target's orbit and the launcher supplies
      all of w, fired anti-parallel to the target velocity.
          cos(psi) = 1  -- every unit of w is retrograde.

  (b) CROSSING: platform and target are on planes separated by angle theta and
      w comes from the geometry, |w| = 2 v sin(theta/2).
          cos(psi) = sin(theta/2)  -- see `retrograde_efficiency`.

The useful-delta-v per unit specific energy is then

    co-orbital : 2*beta/|w|
    crossing   : beta/v_orbital           (independent of theta)

so the co-orbital advantage is exactly 2*v_orbital/|w|.  At v_orbital = 7.45 km/s
and |w| = 619 m/s that is a factor of 24.  This is a purely kinematic result and
it decides the architecture before any guidance argument is made.
"""
import numpy as np

from .constants import MU, R_E


def plane_angle(inc1, inc2, draan):
    """
    Angle between two orbital planes [rad] from inclinations and the RAAN
    difference (spherical law of cosines for the angle between orbit normals).
    """
    i1, i2, dO = map(np.asarray, (inc1, inc2, draan))
    c = np.cos(i1) * np.cos(i2) + np.sin(i1) * np.sin(i2) * np.cos(dO)
    return np.arccos(np.clip(c, -1.0, 1.0))


def v_rel_crossing(v_orbital, theta):
    """|w| for two circular orbits of equal speed whose planes differ by theta."""
    return 2.0 * np.asarray(v_orbital, float) * np.sin(np.asarray(theta, float) / 2.0)


def retrograde_efficiency(theta):
    """
    cos(psi): the fraction of the relative-velocity vector that acts retrograde
    on the target, for a crossing engagement of plane angle theta.

    Derivation. Equal speeds v, velocity vectors separated by theta:
        w . (-u_t) = v(1 - cos theta),   |w| = 2 v sin(theta/2)
        cos psi = v(1-cos theta) / (2 v sin(theta/2))
                = 2 sin^2(theta/2) / (2 sin(theta/2)) = sin(theta/2)
    """
    return np.sin(np.asarray(theta, float) / 2.0)


def useful_dv_per_specific_energy(beta, w, cos_psi):
    """
    (retrograde delta-v delivered) / (specific energy deposited)  [ (m/s) / (J/kg) ].

    Ds_useful = beta m_p |w| cos_psi / m_t ;  E_s = 0.5 m_p |w|^2 / m_t
      -> ratio = 2 beta cos_psi / |w|
    """
    return 2.0 * np.asarray(beta, float) * np.asarray(cos_psi, float) / np.asarray(w, float)


def transverse_rate_crossing(v_orbital, theta):
    """
    Speed at which a crossing target sweeps across the line of sight.  For a
    crossing geometry essentially the whole relative velocity is transverse at
    the moment of closest approach.
    """
    return v_rel_crossing(v_orbital, theta)


def engagement_window(v_transverse, r_track):
    """Time the target spends inside tracking range r_track [s]."""
    v = np.maximum(np.asarray(v_transverse, float), 1e-9)
    return 2.0 * np.asarray(r_track, float) / v


def miss_no_lead(v_transverse, tau, tof):
    """Miss if the shot is aimed at the measured position with no lead."""
    return np.asarray(v_transverse, float) * (np.asarray(tau, float)
                                              + np.asarray(tof, float))


def miss_with_lead(sigma_v, tau, tof, a_rel=0.0):
    """
    Miss when the target state is extrapolated to the intercept epoch.
    Latency enters only through the *estimation* error, not through the
    transverse rate itself.
    """
    t_go = np.asarray(tau, float) + np.asarray(tof, float)
    return np.sqrt((np.asarray(sigma_v, float) * t_go) ** 2
                   + (0.5 * np.asarray(a_rel, float) * t_go ** 2) ** 2)


def sigma_v_from_track(sigma_pos, t_obs, rate_hz):
    """1-sigma velocity estimate from a least-squares fit over a tracking arc."""
    n = np.maximum(np.asarray(t_obs, float) * float(rate_hz), 2.0)
    return np.asarray(sigma_pos, float) * np.sqrt(12.0 / n) / np.asarray(t_obs, float)


def coorbital_drift_rate(dh_m, h_m):
    """
    Along-track drift rate between two co-planar circular orbits separated by
    dh in altitude:  du = |dv_circ/dh| * dh = (v/2a) * dh.
    This is the u that a co-orbital engagement actually sees.
    """
    a = R_E + np.asarray(h_m, float)
    v = np.sqrt(MU / a)
    return v * np.asarray(dh_m, float) / (2.0 * a)


def dv_plane_change(v_orbital, theta):
    """Impulsive delta-v to rotate the velocity vector through theta."""
    return 2.0 * np.asarray(v_orbital, float) * np.sin(np.asarray(theta, float) / 2.0)


def nodal_rate(a, inc, j2=1.08262668e-3):
    """J2 secular nodal regression [rad/s]."""
    a = np.asarray(a, float)
    n = np.sqrt(MU / a ** 3)
    return -1.5 * j2 * (R_E / a) ** 2 * n * np.cos(np.asarray(inc, float))


def differential_nodal_rate(a1, a2, inc1, inc2):
    """Relative RAAN drift rate between two orbits [rad/s]."""
    return nodal_rate(a2, inc2) - nodal_rate(a1, inc1)
