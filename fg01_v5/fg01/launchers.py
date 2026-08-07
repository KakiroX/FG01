"""
Launcher physics: chemical (gunpowder / smokeless) vs. electromagnetic
(railgun, coilgun), for the FG01 momentum-kick projectile.

Every model is a lumped-parameter energy/force balance validated against a
demonstrated benchmark. Sources carry the same tags as the rest of the project:
DERIVED / SOURCED / UNVALIDATED. Citation keys map to docs/CITATION_LOG.md.

The three launcher families are compared on the metrics that actually decide the
choice for an orbital platform:
  - launched system mass (barrel/rails/coils + energy store), because it is
    lifted to orbit;
  - muzzle-velocity repeatability sigma_v/v, because it feeds the SIM-4/INT-4
    aim-error budget through the term u*R*(sigma_v/v)/v_rel and therefore sets
    launched projectile mass quadratically;
  - consumable dependence, because a chemical gun depletes a propellant
    magazine while an EM launcher recharges from the platform's solar array;
  - vacuum / optical-contamination compatibility, because the launcher sits
    metres from a precision terminal sensor;
  - reusable shot life, because the platform fires tens to thousands of times.

Units: SI (m, kg, s, J, Pa, A) unless a name says otherwise.
"""
import numpy as np

from .constants import G0

# ---------------------------------------------------------------------------
# Chemical-propellant properties.  Impetus (force constant) lambda = R*T_v/M
# [J/kg]; gas ratio of specific heats gamma; adiabatic flame temperature T_v.
# SOURCED [carlucci_ballistics] (Carlucci & Jacobson, Ballistics, 2nd ed.);
# black-powder value SOURCED [kirshenbaum].
# ---------------------------------------------------------------------------
PROPELLANTS = {
    "black_powder":   dict(impetus=0.30e6, gamma=1.22, Tv=2000.0,
                           rho=1700.0, note="classical gunpowder; low, dirty, hygroscopic"),
    "single_base_NC": dict(impetus=0.98e6, gamma=1.24, Tv=2600.0,
                           rho=1600.0, note="nitrocellulose smokeless"),
    "double_base":    dict(impetus=1.10e6, gamma=1.25, Tv=3000.0,
                           rho=1650.0, note="NC+NG smokeless, hotter, more erosive"),
}

# Gun ballistic (thermodynamic) efficiency: fraction of propellant potential
# reaching projectile KE.  SOURCED [carlucci_ballistics] typical 0.20-0.35.
GUN_BALLISTIC_ETA = (0.20, 0.35)
GUN_PEAK_TO_MEAN = 2.2          # peak/space-mean pressure ratio, SOURCED (Lagrange)

# Charge-metering and thermal contributions to muzzle-velocity spread.
CHARGE_METER_TOL_KG = 1.0e-3 * 1e-3     # 1 mg absolute metering tolerance (fine)
PROP_TEMP_SENSITIVITY = 0.0012          # fractional dv/v per K, SOURCED (~0.12%/K)
MAGAZINE_THERMAL_K = 5.0                # controlled magazine +-5 K

# ---------------------------------------------------------------------------
# Electromagnetic-launcher properties.
# ---------------------------------------------------------------------------
RAIL_LPRIME = 0.45e-6           # H/m inductance gradient, SOURCED [mcnab2003]
RAIL_ETA = (0.20, 0.35)         # railgun efficiency, SOURCED [mcnab2003, fair]
COIL_ETA = (0.30, 0.50)         # coilgun efficiency, SOURCED [turman1996, kaye]
CAP_ENERGY_DENSITY = (1.5e3, 5.0e3)     # J/kg pulse capacitors, SOURCED range
RAIL_RHO_CU = 8960.0
# muzzle-velocity repeatability, EM: set by capacitor charge-voltage control.
EM_VOLTAGE_TOL = 0.001          # 0.1% voltage control -> ~0.1% velocity (v ~ V)
RAIL_LIFE_SHOTS = (200, 2000)   # rail erosion life at moderate energy, SOURCED
COIL_LIFE_SHOTS = (1e5, 1e7)    # no sliding contact -> bearing/switch limited


def muzzle_energy(m_p, v):
    return 0.5 * np.asarray(m_p, float) * np.asarray(v, float) ** 2


# ---------------------------------------------------------------------------
# Chemical gun
# ---------------------------------------------------------------------------
def chemical_gun(m_p, v, propellant="single_base_NC", eta_bal=0.28,
                 p_peak_MPa=200.0, bore_d=None):
    """
    Interior-ballistics lumped model.

    Propellant mass from the energy balance
        0.5 m_p v^2 = eta_bal * m_prop * lambda/(gamma-1)
    barrel length from the space-mean pressure driving the projectile, and
    barrel mass from a thick-wall (Lame) chamber sized to the peak pressure.
    """
    p = PROPELLANTS[propellant]
    lam, gam = p["impetus"], p["gamma"]
    E = float(muzzle_energy(m_p, v))
    avail_per_kg = lam / (gam - 1.0)
    m_prop = E / (eta_bal * avail_per_kg)

    if bore_d is None:                       # size bore to the projectile
        rho_proj = 19250.0                   # tungsten slug
        bore_d = 2.0 * (3.0 * m_p / (4.0 * np.pi * rho_proj)) ** (1 / 3)
        bore_d = max(bore_d, 4e-3)
    A = np.pi * (bore_d / 2) ** 2
    p_peak = p_peak_MPa * 1e6
    p_mean = p_peak / GUN_PEAK_TO_MEAN
    a_mean = p_mean * A / m_p
    L_barrel = v ** 2 / (2.0 * a_mean)

    # thick-wall chamber: hoop stress, steel yield 1 GPa, SF 2 -> allowable 500 MPa
    sigma_allow = 500e6
    t_wall = (bore_d / 2) * (p_peak / sigma_allow)      # thin-wall estimate
    t_wall = max(t_wall, 2e-3)
    v_steel = np.pi * ((bore_d / 2 + t_wall) ** 2 - (bore_d / 2) ** 2) * L_barrel
    m_barrel = v_steel * 7850.0
    # breech + chamber mass ~ 40% of barrel
    m_gun = m_barrel * 1.4

    # muzzle-velocity repeatability
    charge_frac_tol = CHARGE_METER_TOL_KG / max(m_prop, 1e-9)
    sig_v_charge = 0.5 * charge_frac_tol                # v ~ sqrt(charge)
    sig_v_thermal = PROP_TEMP_SENSITIVITY * MAGAZINE_THERMAL_K
    sig_v_frac = float(np.hypot(sig_v_charge, sig_v_thermal))

    return dict(
        m_prop_kg=m_prop, m_prop_g=m_prop * 1e3,
        E_muzzle_J=E, E_chem_J=m_prop * avail_per_kg,
        eta_overall=E / (m_prop * avail_per_kg),
        bore_mm=bore_d * 1e3, L_barrel_m=L_barrel,
        p_peak_MPa=p_peak_MPa, m_launcher_kg=m_gun,
        m_energy_store_kg=m_prop,          # the charge IS the store
        sigma_v_frac=sig_v_frac,
        recoil_impulse_Ns=m_p * v * 1.3,   # +30% for ejected gas momentum
        consumable=True, clean_vacuum=False,
        shot_life="magazine-limited",
        note=p["note"])


# ---------------------------------------------------------------------------
# Railgun
# ---------------------------------------------------------------------------
def railgun(m_p, v, eta=0.28, L_rail=None, lprime=RAIL_LPRIME,
            cap_density=3.0e3):
    """
    Rail launcher.  Force F = 0.5 L' I^2; choose rail length (or an
    acceleration limit) and solve for current and stored energy.
    """
    E_kin = float(muzzle_energy(m_p, v))
    if L_rail is None:
        a = 5.0e4                            # 5e4 g? no: 5e4 m/s^2 modest
        L_rail = v ** 2 / (2 * a)
    a = v ** 2 / (2 * L_rail)
    F = m_p * a
    I = np.sqrt(2 * F / lprime)
    E_store = E_kin / eta
    m_cap = E_store / cap_density
    # rails: two copper bars, cross-section sized to carry I at ~ 1e9 A/m^2 pulse
    j = 1.0e9
    a_rail = I / j
    m_rails = 2 * a_rail * L_rail * RAIL_RHO_CU
    # containment structure ~ same as rails (magnetic pressure)
    m_struct = m_rails
    m_launcher = m_rails + m_struct
    return dict(
        E_muzzle_J=E_kin, E_store_J=E_store, eta=eta,
        current_kA=I / 1e3, L_rail_m=L_rail, accel_m_s2=a,
        m_launcher_kg=m_launcher, m_energy_store_kg=m_cap,
        m_rails_kg=m_rails,
        sigma_v_frac=EM_VOLTAGE_TOL,
        recoil_impulse_Ns=m_p * v,
        consumable=False, clean_vacuum=True,
        shot_life="rail-erosion limited (200-2000 full-power shots)",
        note="sliding plasma/solid armature; erosion at rail-armature interface")


# ---------------------------------------------------------------------------
# Coilgun (multi-stage induction/reluctance)
# ---------------------------------------------------------------------------
def coilgun(m_p, v, eta=0.40, l_stage=0.30, accel_g=1.0e4,
            cap_density=1.5e3, bore_d=None):
    """
    Multi-stage coilgun.  No sliding contact -> clean, long-lived, precisely
    timed.  Reuses the v5 SIM-6 sizing with the launcher-trade metrics added.
    """
    E_kin = float(muzzle_energy(m_p, v))
    E_store = E_kin / eta
    a = accel_g * G0
    L_barrel = v ** 2 / (2 * a)
    n_stage = max(1, int(np.ceil(L_barrel / l_stage)))
    m_cap = E_store / cap_density
    if bore_d is None:
        rho_proj = 19250.0
        bore_d = max(2.0 * (3.0 * m_p / (4.0 * np.pi * rho_proj)) ** (1 / 3), 4e-3)
    turns_per_m = 200.0
    a_wire = 3.31e-6
    l_wire = np.pi * bore_d * turns_per_m * L_barrel
    m_cu = l_wire * a_wire * RAIL_RHO_CU
    m_launcher = m_cu * 1.6                  # +core/structure
    return dict(
        E_muzzle_J=E_kin, E_store_J=E_store, eta=eta,
        L_barrel_m=L_barrel, n_stages=n_stage, bore_mm=bore_d * 1e3,
        m_launcher_kg=m_launcher, m_energy_store_kg=m_cap, m_cu_kg=m_cu,
        sigma_v_frac=EM_VOLTAGE_TOL * 0.7,   # staged timing averages out ripple
        recoil_impulse_Ns=m_p * v,
        consumable=False, clean_vacuum=True,
        shot_life="1e5-1e7 shots (no sliding contact)",
        note="reluctance/induction stages, capacitor-timed")


def charge_time(E_store_J, p_solar_W=5000.0):
    """Time to recharge the energy store from the platform's solar array."""
    return E_store_J / p_solar_W


def magazine_shots(m_magazine_kg, m_prop_per_shot_kg):
    return m_magazine_kg / m_prop_per_shot_kg
