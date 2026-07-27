"""
CART-01 -- Cartridge charge law: how much rhenium powder to launch.

STANDALONE AND ADDITIVE. This module imports from fg01/ read-only and writes
only cart01_* outputs. It modifies no existing simulation and is deliberately
NOT registered in run_all.py, so a results section being written against the
existing outputs is unaffected.

--------------------------------------------------------------------------
PROBLEM
--------------------------------------------------------------------------
The cartridge charge is adjustable in flight between 50 g and 100 g of rhenium
powder. Given (a) an approximate debris mass and (b) the debris orbit altitude,
how many grams should be loaded?

--------------------------------------------------------------------------
DERIVATION
--------------------------------------------------------------------------
1. Required orbit change. To lower perigee by Delta_h from a circular orbit at
   altitude h (vis-viva, single retrograde burn at apogee):

       dv(h) = sqrt(mu/r) - sqrt(mu*(2/r - 1/a')),   a' = r - Delta_h/2

2. Momentum balance. The mass that actually STRIKES the target must supply the
   target's momentum change:

       m_incident = gamma * m_t * dv(h) / (beta * v_rel)                    (1)

   gamma is the engineering margin, beta the momentum-enhancement factor
   (INT-2: beta = 1.201 at 619 m/s, and velocity dependent).

3. Areal capture. Only the fraction of the cloud incident on the target's
   silhouette lands (INT-1/INT-4). At the cloud width that minimises launched
   mass for a given aim error and reliability:

       P(m_t) = max( K / A_t(m_t), 1 ),   K = -2*pi*e*sigma_miss^2*ln(1-conf)

   so the LAUNCHED powder mass is

       m_L = m_incident * P(m_t)                                            (2)

4. Target geometry. Debris as a solid Al-6061 sphere (AM-2):

       d_t = (6 m_t / (pi rho_Al))^(1/3),   A_t = pi/4 * d_t^2

--------------------------------------------------------------------------
THE RESULTING SCALING -- two regimes
--------------------------------------------------------------------------
Because A_t ~ m_t^(2/3), the launched mass scales as

       m_L ~ m_t^(1/3)      for SMALL targets   (areal penalty dominates)
       m_L ~ m_t            for LARGE targets   (penalty floors at 1)

The crossover is where the target silhouette equals K: at sigma_miss = 5.775 mm
and 95% reliability that is A_t = 1.706e-3 m^2, i.e. d_t = 4.7 cm, m_t = 143 g.

--------------------------------------------------------------------------
THE CONSTRAINT THAT MATTERS MOST
--------------------------------------------------------------------------
Delivering the momentum of eq. (1) puts the impact at

       E_s / E_s,c = gamma * dv * v_rel / (2 * beta * E_s,c)

which is INDEPENDENT of target mass -- 0.67 at the 619 m/s design point. So a
correctly-sized shot is always sub-catastrophic by the same margin.

The danger is the 50 g FLOOR. For a small target the required launched mass is
only a few grams, so firing the minimum cartridge delivers far more momentum
than needed and drives E_s over the catastrophic threshold. The fix is not to
refuse the shot but to DEFOCUS: widen the cloud so that only the required mass
lands and the surplus misses. This module computes the required dispersion.
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd

from fg01.constants import RHO_AL, RHO_RE, GAMMA, ES_C, ES_C_RANGE, OUT
from fg01.orbital import dv_shift, dv_direct_reentry, mass_sphere, xsec_sphere
from fg01.interaction import beta_envelope, mean_incidence_cos_sphere
from fg01.io_utils import save_table, save_json, banner

# ---- cartridge and launcher limits ---------------------------------------
M_MIN = 0.050               # kg, minimum rhenium charge
M_MAX = 0.100               # kg, maximum rhenium charge
M_CANISTER = 0.019          # kg, aluminium canister/armature (COIL-01)
M_ACCEL_MAX = 0.100         # kg, total accelerated mass the coilgun is sized for

# ---- engagement design constants ----------------------------------------
V_REL = 619.0               # m/s, design intercept velocity (CHK-01)
SIGMA_MISS = 0.005775       # m, 1-sigma aim error (SIM-4 / INT-0 budget)
CONF = 0.95                 # per-shot delivery reliability (INT-4)
DH = 200e3                  # m, orbit-lowering target


def beta_of(v_rel=V_REL):
    """INT-2 momentum enhancement, oblique-corrected for a curved target."""
    _, ctr, _ = beta_envelope(float(v_rel))
    return 1.0 + (float(ctr) - 1.0) * mean_incidence_cos_sphere()


def K_areal(sigma_miss=SIGMA_MISS, conf=CONF):
    """Areal constant: penalty = K/A_t (small-target limit), INT-4."""
    return -2.0 * np.pi * np.e * sigma_miss ** 2 * np.log(1.0 - conf)


def target_geometry(m_t):
    """Diameter and silhouette of a solid Al-6061 sphere of mass m_t."""
    d_t = (6.0 * np.asarray(m_t, float) / (np.pi * RHO_AL)) ** (1.0 / 3.0)
    return d_t, np.pi / 4.0 * d_t ** 2


def charge_law(m_t, h_km, v_rel=V_REL, dh=DH, gamma=GAMMA,
               sigma_miss=SIGMA_MISS, conf=CONF, es_c=ES_C):
    """
    THE CHARGE LAW.

    Returns a dict with the powder mass to load and everything needed to
    justify it. m_t in kg, h_km in km.
    """
    m_t = float(m_t)
    dv = float(dv_shift(float(h_km) * 1e3, dh))
    dv_floor = float(dv_direct_reentry(float(h_km) * 1e3))
    beta = beta_of(v_rel)
    d_t, a_t = target_geometry(m_t)

    # (1) mass that must strike the target
    m_inc_req = gamma * m_t * dv / (beta * v_rel)
    # sub-catastrophic ceiling on incident mass
    m_inc_max = 2.0 * es_c * m_t / v_rel ** 2

    # (2) areal penalty at the mass-optimal cloud width
    penalty = max(K_areal(sigma_miss, conf) / a_t, 1.0)
    m_L_ideal = m_inc_req * penalty

    # clamp to the cartridge range
    if m_L_ideal < M_MIN:
        mode = "DEFOCUS (charge floored at 50 g)"
        m_L = M_MIN
        n_shots = 1
        # widen the cloud so only m_inc_req lands.
        # small-target limit: eps = A_t / (2 pi sigma^2)  ->  sigma = sqrt(A_t/(2 pi eps))
        eps_needed = m_inc_req / m_L
        sigma_cloud = np.sqrt(a_t / (2.0 * np.pi * eps_needed))
    elif m_L_ideal <= M_MAX:
        mode = "NOMINAL"
        m_L = m_L_ideal
        n_shots = 1
        eps_needed = m_inc_req / m_L
        sigma_cloud = sigma_miss * np.sqrt(-np.log(1.0 - conf))
    else:
        mode = "MULTI-SHOT (charge capped at 100 g)"
        m_L = M_MAX
        n_shots = int(np.ceil(m_L_ideal / M_MAX))
        eps_needed = min(m_inc_req / (m_L * n_shots), 1.0)
        sigma_cloud = sigma_miss * np.sqrt(-np.log(1.0 - conf))

    m_inc_actual = m_L * n_shots * min(eps_needed, 1.0)
    es = 0.5 * m_inc_actual * v_rel ** 2 / m_t

    return dict(
        m_debris_kg=m_t, h_km=float(h_km), d_debris_cm=d_t * 100,
        dv_required_m_s=dv, beta=beta,
        m_incident_required_g=m_inc_req * 1e3,
        m_incident_max_subcat_g=m_inc_max * 1e3,
        areal_penalty=penalty,
        m_powder_ideal_g=m_L_ideal * 1e3,
        m_powder_g=m_L * 1e3, n_shots=n_shots, mode=mode,
        cloud_sigma_mm=sigma_cloud * 1e3,
        cloud_diameter_cm=sigma_cloud * 2 * 2.4477 * 100,
        Es_kJ_kg=es / 1e3, Es_over_threshold=es / es_c,
        subcatastrophic=bool(es < es_c),
        subcat_worst_case=bool(es < ES_C_RANGE[0]),
        total_accel_mass_g=(m_L + M_CANISTER) * 1e3,
        within_coilgun_limit=bool(m_L + M_CANISTER <= M_ACCEL_MAX),
        self_disposal_ok=bool(v_rel >= dv_floor),
        v_rel_m_s=float(v_rel))


def envelope_bounds(h_km, v_rel=V_REL, gamma=GAMMA, dh=DH,
                    sigma_miss=SIGMA_MISS, conf=CONF):
    """
    Debris-mass range a SINGLE cartridge can service at altitude h.

    Solves m_L(m_t) = M_MIN and = M_MAX. In the large-target regime the penalty
    is 1 and the relation is linear, so it inverts directly; in the small-target
    regime it goes as m_t^(1/3). Solved numerically to cover both.
    """
    def f(m_t):
        return charge_law(m_t, h_km, v_rel, dh, gamma, sigma_miss,
                          conf)["m_powder_ideal_g"] / 1e3

    def solve(target):
        lo, hi = 1e-6, 50.0
        for _ in range(200):
            mid = np.sqrt(lo * hi)
            if f(mid) < target:
                lo = mid
            else:
                hi = mid
        return np.sqrt(lo * hi)

    return solve(M_MIN), solve(M_MAX)


def run():
    banner("CART-01  Cartridge charge law (50-100 g adjustable)")

    beta = beta_of()
    K = K_areal()
    d_x, _ = target_geometry(0.143)
    print(f"\n[0] CONSTANTS")
    print(f"    v_rel                    {V_REL:.0f} m/s")
    print(f"    beta (INT-2)             {beta:.4f}")
    print(f"    aim error sigma_miss     {SIGMA_MISS*1e3:.3f} mm")
    print(f"    reliability              {CONF:.0%}")
    print(f"    areal constant K         {K:.4e} m^2")
    print(f"    penalty=1 crossover      A_t = {K:.3e} m^2 -> "
          f"d = {2*np.sqrt(K/np.pi)*100:.1f} cm, m = "
          f"{float(mass_sphere(2*np.sqrt(K/np.pi), RHO_AL))*1e3:.0f} g")

    # ---- the closed form, with numbers -----------------------------------
    print(f"\n[1] THE CHARGE LAW\n")
    print("    m_powder = clamp[ gamma * m_t * dv(h) / (beta * v_rel)"
          " * P(m_t) , 50 g, 100 g ]")
    print("    P(m_t)   = max( K / A_t , 1 ),   K = -2*pi*e*sigma^2*ln(1-conf)")
    print("    A_t      = pi/4 * (6 m_t / (pi rho_Al))^(2/3)")
    print(f"\n    with gamma={GAMMA}, beta={beta:.4f}, v_rel={V_REL:.0f} m/s, "
          f"K={K:.4e} m^2:")
    for h in (600, 900, 1200):
        dv = float(dv_shift(h * 1e3, DH))
        cp = GAMMA * dv / (beta * V_REL)
        print(f"      h={h:4d} km:  dv={dv:5.2f} m/s   "
              f"m_powder = {cp:.5f} * m_t * P(m_t)")

    # ---- single-cartridge envelope ---------------------------------------
    print(f"\n[2] DEBRIS-MASS ENVELOPE FOR ONE CARTRIDGE\n")
    rows = []
    for h in (400, 500, 600, 700, 800, 900, 1000, 1100, 1200):
        lo, hi = envelope_bounds(h)
        d_lo, _ = target_geometry(lo)
        d_hi, _ = target_geometry(hi)
        rows.append(dict(h_km=h,
                         dv_required_m_s=float(dv_shift(h * 1e3, DH)),
                         m_debris_min_g=lo * 1e3, m_debris_max_g=hi * 1e3,
                         d_debris_min_cm=d_lo * 100, d_debris_max_cm=d_hi * 100))
    d_env = pd.DataFrame(rows)
    save_table(d_env, "cart01_envelope",
               "Debris-mass range one 50-100 g cartridge can service, by "
               "altitude. Below the minimum the charge floors at 50 g and the "
               "cloud must be defocused; above the maximum multiple shots are "
               "needed.")
    print(d_env.to_string(index=False, float_format=lambda x: f"{x:.2f}"))

    # ---- lookup table -----------------------------------------------------
    print(f"\n[3] LOOKUP TABLE: powder grams by debris mass and altitude\n")
    masses = [0.0014, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0]
    alts = [400, 600, 800, 900, 1000, 1200]
    rows = []
    for m_t in masses:
        for h in alts:
            r = charge_law(m_t, h)
            rows.append(r)
    d_lut = pd.DataFrame(rows)
    save_table(d_lut, "cart01_lookup",
               "Charge law evaluated over debris mass and altitude. "
               "m_powder_g is what to load; mode says whether the cartridge "
               "floor or ceiling is binding.")
    piv = d_lut.pivot_table(index="m_debris_kg", columns="h_km",
                            values="m_powder_g")
    print("    powder to load (g):")
    print(piv.to_string(float_format=lambda x: f"{x:.1f}"))
    piv2 = d_lut.pivot_table(index="m_debris_kg", columns="h_km",
                             values="n_shots")
    print("\n    shots required:")
    print(piv2.to_string(float_format=lambda x: f"{x:.0f}"))

    # ---- the small-target hazard -----------------------------------------
    print(f"\n[4] THE 50 g FLOOR IS A HAZARD FOR SMALL DEBRIS\n")
    rows = []
    for m_t in (0.0014, 0.005, 0.02, 0.05, 0.1, 0.15, 0.24, 0.36):
        # what a 50 g cartridge does if fired at the mass-optimal cloud width
        d_t, a_t = target_geometry(m_t)
        pen = max(K_areal() / a_t, 1.0)
        m_inc = M_MIN / pen
        es = 0.5 * m_inc * V_REL ** 2 / m_t
        r = charge_law(m_t, 900)
        rows.append(dict(m_debris_g=m_t * 1e3, d_debris_cm=d_t * 100,
                         incident_if_focused_g=m_inc * 1e3,
                         Es_if_focused_kJ_kg=es / 1e3,
                         times_threshold=es / ES_C,
                         catastrophic_if_focused=bool(es >= ES_C),
                         required_incident_g=r["m_incident_required_g"],
                         defocus_cloud_cm=r["cloud_diameter_cm"],
                         Es_after_defocus_kJ_kg=r["Es_kJ_kg"]))
    d_haz = pd.DataFrame(rows)
    save_table(d_haz, "cart01_small_target_hazard",
               "Firing the 50 g minimum at a small target with a normally-"
               "focused cloud is catastrophic. The required defocus is given.")
    print(d_haz.to_string(index=False, float_format=lambda x: f"{x:.2f}"))

    # ---- validation against the existing modules -------------------------
    print(f"\n[5] VALIDATION against INT-3/INT-4 (existing outputs, unchanged)\n")
    m_t_1cm = float(mass_sphere(0.01, RHO_AL))
    r = charge_law(m_t_1cm, 900)
    a_t = float(xsec_sphere(0.01))
    pen_ref = K_areal() / a_t
    print(f"    1 cm target at 900 km:")
    print(f"      areal penalty  this module {r['areal_penalty']:.3f}  "
          f"vs INT-4 21.727  ({abs(r['areal_penalty']-21.727)/21.727:.2%})")
    m_on_int3 = GAMMA * float(dv_shift(900e3, DH)) * m_t_1cm / (beta * V_REL)
    print(f"      mass on target this module {r['m_incident_required_g']:.4f} g "
          f"vs INT-3 0.19709 g  "
          f"({abs(r['m_incident_required_g']-0.19709)/0.19709:.2%})")
    print(f"      ideal powder   {r['m_powder_ideal_g']:.3f} g "
          f"vs INT-4 launched 4.282 g "
          f"({abs(r['m_powder_ideal_g']-4.282)/4.282:.2%})")

    # ---- IS DOSING WORTH IT AT ALL? --------------------------------------
    print(f"\n[7] CONSTANT CHARGE vs DOSED -- is metering worth the mechanism?\n")
    M_MAG = 30.0                # kg magazine (SIM-11 reference platform)
    N_ENG_DECADE = 24           # ECO-1 realised engagements per platform-decade
    SHOTS_PER_ENG = 2           # allow for re-engagement / partial delivery
    need = N_ENG_DECADE * SHOTS_PER_ENG

    rows = []
    for label, charge in (("constant 100 g", 0.100),
                          ("constant 81 g (fits 100 g total)", 0.081),
                          ("constant 50 g (floor)", 0.050),
                          ("dosed 50-100 g, envelope avg", 0.075)):
        shots = M_MAG / charge
        m_tot = charge + M_CANISTER
        # muzzle energy scales linearly with mass at fixed velocity
        e_muzzle = 0.5 * m_tot * 400.0 ** 2
        e_cap = e_muzzle / 0.196          # COIL-01 efficiency
        # largest single-shot target (penalty = 1 regime), 900 km
        cp = GAMMA * float(dv_shift(900e3, DH)) / (beta_of() * V_REL)
        m_t_max = charge / cp
        # smallest target that is safe WITHOUT defocusing
        m_t_safe = 0.5 * charge * V_REL ** 2 / ES_C
        rows.append(dict(
            scheme=label, charge_g=charge * 1e3,
            total_accel_g=m_tot * 1e3,
            within_100g_total=bool(m_tot <= M_ACCEL_MAX),
            cap_energy_kJ=e_cap / 1e3,
            cap_mass_kg=e_cap / 1500.0,
            shots_per_30kg=shots,
            shots_needed_decade=need,
            margin_x=shots / need,
            max_single_shot_target_g=m_t_max * 1e3,
            min_safe_target_no_defocus_g=m_t_safe * 1e3))
    d_cd = pd.DataFrame(rows)
    save_table(d_cd, "cart01_constant_vs_dosed",
               "Constant charge vs dosed. Magazine 30 kg (SIM-11), engagement "
               "demand 24 per decade x2 shots (ECO-1). Shows that shot count is "
               "not the binding constraint, so metering buys little.")
    print(d_cd.to_string(index=False, float_format=lambda x: f"{x:.1f}"))
    print(f"\n    magazine gives 300-600 shots; ECO-1 demand is ~{need} per "
          f"decade -> {300/need:.0f}-{600/need:.0f}x margin either way.")
    print(f"    So dosing does NOT buy needed shots. Its only real value would be")
    print(f"    rhenium economy: 100 g vs 75 g avg = 25 g/shot = "
          f"${25*7.283/1000:.2f}/shot, i.e. ${need*25*7.283/1000:.0f} per decade.")

    concl = {
        "formula": ("m_powder = clamp[ gamma*m_t*dv(h)/(beta*v_rel) * P(m_t),"
                    " 50 g, 100 g ];  P = max(K/A_t, 1);"
                    " K = -2*pi*e*sigma_miss^2*ln(1-conf)"),
        "dosing_verdict": (
            "Charge metering is NOT required. Because delivered momentum is set "
            "by cloud width times charge, dispersion control alone can size any "
            "shot, and the magazine already provides 300-600 shots against an "
            "ECO-1 demand of ~48 per decade. A CONSTANT 81 g charge is "
            "recommended: it keeps total accelerated mass at exactly the 100 g "
            "the coilgun is sized for (42 kJ, 28 kg of capacitors), removes the "
            "metering mechanism, and still covers single-shot targets to 581 g. "
            "A constant 100 g charge would push total mass to 119 g and require "
            "50 kJ / 34 kg of capacitors -- paying 6 kg permanently to avoid a "
            "mechanism it does not need. Dispersion control remains mandatory "
            "either way."),
        "constants": dict(gamma=GAMMA, beta=beta, v_rel=V_REL,
                          sigma_miss_mm=SIGMA_MISS * 1e3, conf=CONF,
                          K_m2=K),
        "scaling": ("m_powder ~ m_t^(1/3) below 143 g debris (areal penalty "
                    "dominates), ~m_t above it (penalty floors at 1)"),
        "envelope_900km_g": list(np.array(envelope_bounds(900)) * 1e3),
        "Es_over_threshold_when_correctly_sized": float(
            GAMMA * float(dv_shift(900e3, DH)) * V_REL / (2 * beta * ES_C)),
        "small_target_warning": (
            "Below ~240 g debris the 50 g cartridge floor delivers more "
            "momentum than the target can absorb sub-catastrophically if the "
            "cloud is normally focused. Defocus per cart01_small_target_hazard."),
        "canister_note": (
            f"50-100 g of powder plus the {M_CANISTER*1e3:.0f} g aluminium "
            f"canister gives {(M_MIN+M_CANISTER)*1e3:.0f}-"
            f"{(M_MAX+M_CANISTER)*1e3:.0f} g accelerated. The COIL-01 design "
            "was sized for 100 g TOTAL, so a full 100 g powder charge is 19% "
            "over that and needs ~50 kJ rather than 42 kJ."),
    }
    save_json(concl, "cart01_conclusions")
    print(f"\n[6] NOTE ON TOTAL ACCELERATED MASS")
    print(f"    powder 50-100 g + canister {M_CANISTER*1e3:.0f} g = "
          f"{(M_MIN+M_CANISTER)*1e3:.0f}-{(M_MAX+M_CANISTER)*1e3:.0f} g accelerated")
    print(f"    COIL-01 was sized for 100 g TOTAL -> a full 100 g powder charge")
    print(f"    is 19% over and needs ~50 kJ rather than 42 kJ.")
    return d_lut, d_env, concl


if __name__ == "__main__":
    run()
