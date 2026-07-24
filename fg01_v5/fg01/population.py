"""
Debris-population and breakup-cloud model (ECO-2).

Real cluster parameters, with the provenance of each number recorded.  The
central difficulty of this module is stated up front: the >=1 cm population is
*not* catalogued.  Roughly 1.2 million objects >=1 cm are estimated to be in
orbit against ~40,000 tracked, so every member count below is a modelled
estimate, not an observation, and is tagged UNVALIDATED accordingly.

The generator produces a synthetic fragment cloud from a parent breakup using
the NASA SBM delta-v distribution, then ages it under differential J2 nodal
regression, which is what turns a compact cloud into a shell and drives the
whole economic result.
"""
import numpy as np

from .constants import MU, R_E, J2, RHO_AL
from .relmotion import nodal_rate, plane_angle

# --- Major LEO breakup clouds ------------------------------------------------
# n_ge_1cm entries are MODELLED estimates (UNVALIDATED); tracked counts and
# orbital parameters are SOURCED.  See docs/CITATION_LOG.md.
CLUSTERS = {
    "Fengyun-1C": dict(
        year=2007, h_km=865.0, inc_deg=98.8, n_tracked=3500, n_ge_1cm=35000,
        note="2007 ASAT test, sun-synchronous. Largest single debris event.",
        source="SOURCED: 3500 tracked, 865 km, i=98.8 deg; "
               "n_ge_1cm 35000-40000 is a MODELLED estimate"),
    "Cosmos-2251": dict(
        year=2009, h_km=789.0, inc_deg=74.0, n_tracked=1668, n_ge_1cm=17000,
        note="2009 accidental collision with Iridium-33.",
        source="SOURCED: 1668 catalogued; n_ge_1cm scaled from tracked count"),
    "Iridium-33": dict(
        year=2009, h_km=789.0, inc_deg=86.4, n_tracked=628, n_ge_1cm=6500,
        note="2009 accidental collision with Cosmos-2251.",
        source="SOURCED: 628 catalogued; n_ge_1cm scaled"),
    "Cosmos-1408": dict(
        year=2021, h_km=480.0, inc_deg=82.5, n_tracked=1500, n_ge_1cm=15000,
        note="2021 ASAT test. Low altitude -- decays fast unaided.",
        source="SOURCED: >=1500 tracked, 465-490 km, i=82.5 deg"),
}

# ESA MASTER, August 2024 (SOURCED)
MASTER_TOTAL_GE_1CM = 1.2e6
MASTER_TOTAL_GE_10CM = 54000
MASTER_TRACKED = 40000


def sbm_delta_v(n, lc_m, rng, catastrophic=True):
    """
    NASA SBM fragment ejection-speed distribution.  log10(dv) is normal with
    mean mu(chi) and sigma 0.4, where chi = log10(A/m).  For the size range of
    interest this yields ejection speeds of tens to hundreds of m/s.
    """
    chi = np.log10(3.0 / (4.0 * RHO_AL * (np.asarray(lc_m, float) / 2.0)))
    mu = (0.9 * chi + 2.9) if catastrophic else (0.9 * chi + 2.9)
    return 10 ** rng.normal(mu, 0.4, n)


def synth_cloud(cluster, n_frag, rng, lc_m=0.01):
    """
    Generate a synthetic fragment cloud: apply an isotropic SBM delta-v to the
    parent state and convert to orbital elements (a, e, i, RAAN).
    """
    c = CLUSTERS[cluster] if isinstance(cluster, str) else cluster
    a0 = R_E + c["h_km"] * 1e3
    v0 = np.sqrt(MU / a0)
    inc0 = np.radians(c["inc_deg"])

    dv = sbm_delta_v(n_frag, lc_m, rng)
    # isotropic directions
    u = rng.normal(size=(n_frag, 3))
    u /= np.linalg.norm(u, axis=1, keepdims=True)
    dvv = dv[:, None] * u

    # parent state at the ascending node, velocity in the orbital plane
    r_vec = np.array([a0, 0.0, 0.0])
    v_vec = np.array([0.0, v0 * np.cos(inc0), v0 * np.sin(inc0)])
    v_new = v_vec[None, :] + dvv

    r_n = np.linalg.norm(r_vec)
    v_n = np.linalg.norm(v_new, axis=1)
    energy = 0.5 * v_n ** 2 - MU / r_n
    a = -MU / (2 * energy)
    h_vec = np.cross(np.broadcast_to(r_vec, v_new.shape), v_new)
    h_n = np.linalg.norm(h_vec, axis=1)
    e_vec = np.cross(v_new, h_vec) / MU - r_vec / r_n
    e = np.linalg.norm(e_vec, axis=1)
    inc = np.arccos(np.clip(h_vec[:, 2] / h_n, -1, 1))
    raan = np.arctan2(h_vec[:, 0], -h_vec[:, 1])

    ok = (a > 0) & (a * (1 - e) > R_E + 150e3) & np.isfinite(a)
    return dict(a=a[ok], e=e[ok], inc=inc[ok], raan=raan[ok] % (2 * np.pi),
                n_kept=int(ok.sum()), n_requested=n_frag)


def age_cloud(cloud, years):
    """Advance RAAN by differential J2 nodal regression over `years`."""
    rate = nodal_rate(cloud["a"], cloud["inc"])
    out = dict(cloud)
    out["raan"] = (cloud["raan"] + rate * years * 365.25 * 86400.0) % (2 * np.pi)
    return out


def engageable_mask(cloud, a_plat, inc_plat, raan_plat, theta_max_rad,
                    dh_band_m=50e3):
    """
    Fragments a platform can engage: plane angle within theta_max AND altitude
    within a band the platform can reach cheaply.
    """
    th = plane_angle(inc_plat, cloud["inc"], cloud["raan"] - raan_plat)
    in_plane = th <= theta_max_rad
    in_alt = np.abs(cloud["a"] - a_plat) <= dh_band_m
    return in_plane & in_alt, th


def scale_to_population(cloud, n_total_ge_1cm):
    """Weight each synthetic fragment so the cloud represents the real count."""
    return n_total_ge_1cm / max(cloud["n_kept"], 1)
