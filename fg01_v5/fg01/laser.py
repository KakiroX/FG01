"""
Laser-ablation debris removal: the competing method, modelled on the same
footing as FG01 so the comparison is like-for-like.

The v5/v6 cost comparison used published figures of $100-500/kg for ground-based
laser ablation, tagged "theoretical".  Those figures were never re-derived, and
FG01 was declared to lose to them by five orders of magnitude.  This module
builds the laser engagement from physics so that both systems are judged by the
same standard.

Three things distinguish laser ablation from a momentum kick, and all three are
implemented explicitly:

1. THE FLUENCE THRESHOLD.  Ablation is a threshold process.  Below F_th the
   surface is heated and nothing is ejected, so the delivered impulse is
   exactly zero -- not small, zero.  Spreading a laser spot does not buy
   tolerance the way spreading a grain cloud does; it buys nothing until the
   fluence recovers.  This is the formal statement of "miss entirely, no
   effect" versus "partial effect".

2. THE BEAM IS AIMED, THE IMPULSE IS NOT.  Ablation recoil is directed along
   the surface normal / beam axis, i.e. away from the illuminator.  For a
   ground station that is mostly radial.  Only the component opposite the
   target's velocity lowers the orbit, and that component is largest exactly
   where the range and the airmass are worst.

3. THE ATMOSPHERE.  Turbulence limits the achievable spot to the Fried
   parameter unless adaptive optics restores it, and both the transmission and
   the seeing degrade as the airmass rises -- which is precisely the geometry
   the useful-impulse term prefers.

Momentum coupling model (Phipps et al., laser propulsion literature): for
nanosecond pulses on metals the coupling coefficient C_m rises from zero at the
ablation threshold to a maximum near an optimum fluence and falls off as
Phi^-1/2 beyond it.  C_m,max ~ 1e-4 N.s/J for aluminium.
"""
import numpy as np

C_LIGHT = 2.99792458e8
SOLAR_CONST = 1361.0            # W/m^2 at 1 AU
M_SUN_V = -26.74                # apparent V magnitude of the Sun

# --- ablation parameters for Al-6061, ns pulses -----------------------------
F_TH = 1.0e4            # J/m^2  = 1 J/cm^2   ablation threshold   SOURCED range 0.5-5
F_OPT = 5.0e4           # J/m^2  = 5 J/cm^2   optimum coupling     SOURCED
CM_MAX = 1.0e-4         # N.s/J  peak momentum coupling            SOURCED


def spot_diameter(R, D, lam, beam_quality=2.0, strehl=1.0):
    """
    Far-field spot diameter [m] at range R for aperture D.

    2.44*lam*R/D is the first-null diameter of an ideal circular aperture;
    beam_quality (M^2-like) and the adaptive-optics Strehl ratio degrade it.
    """
    d = 2.44 * np.asarray(lam, float) * np.asarray(R, float) \
        * np.asarray(beam_quality, float) / np.asarray(D, float)
    return d / np.sqrt(np.clip(strehl, 1e-3, 1.0))


def airmass(elev_deg):
    """Kasten-Young airmass."""
    e = np.asarray(elev_deg, float)
    return 1.0 / (np.sin(np.radians(e))
                  + 0.50572 * (e + 6.07995) ** -1.6364)


def atmospheric_transmission(elev_deg, tau_zenith=0.15):
    return np.exp(-tau_zenith * airmass(elev_deg))


def fried_parameter(elev_deg, r0_zenith=0.15, lam=532e-9, lam_ref=500e-9):
    """Fried parameter r0 [m], scaled for wavelength and airmass."""
    return (r0_zenith * (np.asarray(lam, float) / lam_ref) ** 1.2
            * airmass(elev_deg) ** -0.6)


def ao_strehl(D, r0, n_actuators=1000, loop_hz=1000.0, wind_m_s=20.0,
              theta_off_arcsec=0.0, theta_iso_arcsec=2.5):
    """
    Strehl ratio for a real adaptive-optics system, Marechal approximation over
    the standard error budget.  Deliberately generous to the laser: the
    dominant term is the *fitting* error set by actuator spacing, not the full
    uncorrected D/r0 variance, so a well-equipped 5-10 m telescope reaches
    Strehl 0.5-0.8 rather than zero.

        sigma^2_fit    = 0.34 (d_act / r0)^(5/3),  d_act = D / sqrt(N_act)
        sigma^2_servo  = (f_G / f_loop)^(5/3),     f_G = 0.43 v_wind / r0
        sigma^2_aniso  = (theta / theta_0)^(5/3)

    n_actuators = 1000 across a 5 m aperture is comparable to deployed
    facility-class AO.
    """
    D = np.asarray(D, float)
    r0 = np.asarray(r0, float)
    d_act = D / np.sqrt(float(n_actuators))
    s2_fit = 0.34 * (d_act / r0) ** (5.0 / 3.0)
    f_g = 0.43 * wind_m_s / r0
    s2_servo = (f_g / float(loop_hz)) ** (5.0 / 3.0)
    s2_aniso = (float(theta_off_arcsec) / float(theta_iso_arcsec)) ** (5.0 / 3.0)
    return np.exp(-np.clip(s2_fit + s2_servo + s2_aniso, 0, 50))


def fluence(E_pulse, spot_d, transmission=1.0):
    """Peak fluence [J/m^2] for pulse energy E over a spot of diameter spot_d."""
    area = np.pi * (np.asarray(spot_d, float) / 2.0) ** 2
    return np.asarray(E_pulse, float) * np.asarray(transmission, float) / area


def coupling_coefficient(F, f_th=F_TH, f_opt=F_OPT, cm_max=CM_MAX):
    """
    Momentum coupling C_m [N.s/J] as a function of fluence.

    Zero below threshold -- this is the cliff.  Linear ramp from threshold to
    the optimum, then the standard Phi^-1/2 roll-off.
    """
    F = np.asarray(F, float)
    ramp = cm_max * np.clip((F - f_th) / (f_opt - f_th), 0.0, 1.0)
    roll = cm_max * np.sqrt(f_opt / np.maximum(F, 1e-30))
    return np.where(F < f_th, 0.0, np.where(F <= f_opt, ramp, roll))


def impulse_per_pulse(E_pulse, spot_d, a_target, transmission=1.0, **kw):
    """
    Impulse [N.s] delivered to a target of cross-section a_target.

    Only the fraction of the beam intercepted by the target ablates, and the
    coupling is evaluated at the local fluence -- so a spot that is too large
    fails twice over: less energy intercepted AND a lower (possibly zero)
    coupling coefficient.
    """
    F = fluence(E_pulse, spot_d, transmission)
    cm = coupling_coefficient(F, **kw)
    spot_area = np.pi * (np.asarray(spot_d, float) / 2.0) ** 2
    e_on_target = (np.asarray(E_pulse, float) * np.asarray(transmission, float)
                   * np.minimum(np.asarray(a_target, float) / spot_area, 1.0))
    return cm * e_on_target, F, cm


# --------------------------------------------------------------------------
# Ground-station engagement geometry
# --------------------------------------------------------------------------
def pass_geometry(h_m, phi_rad, r_e=6.371e6):
    """
    Geometry of a circular-orbit target seen from a ground station.

    phi is the geocentric central angle between station and target; negative
    means the target is approaching (it has not yet reached the station's
    meridian).

    Returns (slant_range, elevation_deg, retrograde_efficiency).

    Retrograde efficiency is the fraction of the ablation impulse -- directed
    along the line of sight, away from the station -- that opposes the target's
    velocity.  With the station at (r_e, 0) and the target at
    (r cos phi, r sin phi) moving in +phi:

        LOS_hat . v_hat = r_e sin(phi) / R

    so the useful (retrograde) fraction is -r_e sin(phi)/R, positive while the
    target approaches and negative once it recedes.  It vanishes at zenith.
    """
    r = r_e + np.asarray(h_m, float)
    phi = np.asarray(phi_rad, float)
    R = np.sqrt(r_e ** 2 + r ** 2 - 2 * r_e * r * np.cos(phi))
    sin_elev = np.clip((r * np.cos(phi) - r_e) / R, -1.0, 1.0)
    elev = np.degrees(np.arcsin(sin_elev))
    retro = -r_e * np.sin(phi) / R
    return R, elev, retro


def angular_rate(h_m, R, mu=3.986004418e14, r_e=6.371e6):
    """Apparent transverse angular rate of the target as seen from the ground."""
    r = r_e + np.asarray(h_m, float)
    v = np.sqrt(mu / r)
    return v / np.asarray(R, float)


def point_ahead_angle(v_transverse, R):
    """Point-ahead angle for a round-trip light time 2R/c."""
    return 2.0 * np.asarray(v_transverse, float) / C_LIGHT


# --------------------------------------------------------------------------
# Detection photometry -- can the target even be seen?
# --------------------------------------------------------------------------
def apparent_magnitude(d_target, R, albedo=0.1, phase=1.0):
    """Apparent V magnitude of a diffuse sphere at range R."""
    a_t = np.pi * (np.asarray(d_target, float) / 2.0) ** 2
    intensity = albedo * SOLAR_CONST * a_t * phase / np.pi     # W/sr
    irradiance = intensity / np.asarray(R, float) ** 2
    return M_SUN_V - 2.5 * np.log10(irradiance / SOLAR_CONST)


def photon_rate(d_target, R, D_receive, albedo=0.1, lam=550e-9,
                throughput=0.4, qe=0.8):
    """Detected photo-electrons per second from a sunlit diffuse sphere."""
    a_t = np.pi * (np.asarray(d_target, float) / 2.0) ** 2
    intensity = albedo * SOLAR_CONST * a_t / np.pi
    irradiance = intensity / np.asarray(R, float) ** 2
    a_rec = np.pi * (np.asarray(D_receive, float) / 2.0) ** 2
    power = irradiance * a_rec
    e_photon = 6.62607015e-34 * C_LIGHT / lam
    return power / e_photon * throughput * qe


def resolved_ratio(d_target, R, D, lam):
    """
    Target angular size divided by the diffraction limit of the observing
    aperture.  >1 means the target is resolved; <1 means it is an unresolved
    point source and its position must be inferred by centroiding.
    """
    theta_t = np.asarray(d_target, float) / np.asarray(R, float)
    theta_d = 1.22 * np.asarray(lam, float) / np.asarray(D, float)
    return theta_t / theta_d
