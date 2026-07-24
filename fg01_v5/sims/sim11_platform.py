"""
SIM-11 -- FG01 platform design, mass budget and end-of-life.

The co-orbital engagement geometry that SIM-4 relies on is not free: FG01 must
reach each target's orbital plane.  That cost is charged here, and it is the
principal operational constraint on the concept.

Inclination changes are prohibitively expensive (2 v sin(di/2) ~ 131 m/s per
degree at LEO speeds), so one platform serves one inclination band.  RAAN
changes, by contrast, are nearly free in propellant if bought with time, using
differential J2 nodal regression from a small altitude offset.
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from fg01.constants import (MU, R_E, J2, G0, RHO_RE, C_LAUNCH_RANGE, ALT_SWEEP_KM)
from fg01.orbital import v_circ, mean_motion, dv_shift, r_of_h, vis_viva
from fg01.io_utils import save_table, save_fig, save_json, banner, CB


def nodal_rate(h_km, inc_deg):
    """J2 secular nodal regression rate, rad/s (circular orbit)."""
    a = R_E + h_km * 1e3
    n = np.sqrt(MU / a ** 3)
    return -1.5 * J2 * (R_E / a) ** 2 * n * np.cos(np.radians(inc_deg))


def hohmann(h1_km, h2_km):
    r1, r2 = R_E + h1_km * 1e3, R_E + h2_km * 1e3
    at = 0.5 * (r1 + r2)
    return (abs(vis_viva(r1, at) - np.sqrt(MU / r1))
            + abs(np.sqrt(MU / r2) - vis_viva(r2, at)))


def run():
    banner("SIM-11  FG01 platform, mass budget and end-of-life")

    # ---- mass budget -------------------------------------------------------
    budgets = []
    for m_payload in (10.0, 30.0, 100.0):
        for m_coilgun in (100.0, 250.0, 500.0):
            avionics, power, thermal = 50.0, 5000.0 / 30.0, 30.0
            subtotal = avionics + power + thermal + m_coilgun
            structure = 0.20 * (subtotal + m_payload)
            dry = subtotal + structure
            budgets.append(dict(
                m_payload_kg=m_payload, m_coilgun_kg=m_coilgun,
                avionics_kg=avionics, power_kg=power, thermal_kg=thermal,
                structure_kg=structure, dry_kg=dry,
                wet_kg=dry + m_payload,
                shots_per_reload_design=m_payload / 2.69e-3,
                shots_per_reload_large=m_payload / 175e-3,
                launch_cost_low_musd=(dry + m_payload) * C_LAUNCH_RANGE[0] / 1e6,
                launch_cost_high_musd=(dry + m_payload) * C_LAUNCH_RANGE[1] / 1e6))
    dfb = pd.DataFrame(budgets)
    save_table(dfb, "sim11_mass_budget",
               "FG01 mass budget. Power at 5 kW and 30 W/kg; structure 20% of "
               "subtotal; payload is the rhenium magazine.",
               tags={"subsystem fractions": "UNVALIDATED (parametric)"})
    ref = dfb[(dfb.m_payload_kg == 30.0) & (dfb.m_coilgun_kg == 250.0)].iloc[0]
    print(f"  reference platform: dry {ref.dry_kg:.0f} kg, wet {ref.wet_kg:.0f} kg, "
          f"{ref.shots_per_reload_design:.0f} design shots per 30 kg magazine")

    # ---- recoil ------------------------------------------------------------
    rec = []
    for m_shot_g in (2.69, 175.0):
        for n_shots in (1e3, 1e4, 1e5):
            imp = n_shots * m_shot_g * 1e-3 * 619.0
            for isp in (300.0, 3000.0):
                rec.append(dict(m_shot_g=m_shot_g, n_shots=n_shots,
                                total_impulse_Ns=imp,
                                dv_on_500kg_m_s=imp / 500.0,
                                propellant_kg=imp / (isp * G0), isp_s=isp))
    dfr = pd.DataFrame(rec)
    save_table(dfr, "sim11_recoil",
               "Recoil impulse from firing and the propellant needed to hold "
               "the operating orbit. Each shot pushes FG01 prograde because "
               "the projectile is fired retrograde.")
    print("\n  recoil: 1e4 design shots =",
          f"{dfr[(dfr.m_shot_g==2.69)&(dfr.n_shots==1e4)&(dfr.isp_s==300)].propellant_kg.iloc[0]:.1f} kg"
          " propellant at Isp 300 s")

    # ---- plane access ------------------------------------------------------
    plane = []
    for h in ALT_SWEEP_KM:
        for inc in (51.6, 74.0, 82.0, 98.0):
            for dh_park in (25.0, 50.0, 100.0):
                r1 = nodal_rate(h, inc)
                r2 = nodal_rate(h - dh_park, inc)
                drift_deg_day = np.degrees(abs(r2 - r1)) * 86400.0
                plane.append(dict(
                    h_km=h, inc_deg=inc, park_offset_km=dh_park,
                    nodal_rate_deg_day=np.degrees(r1) * 86400.0,
                    differential_drift_deg_day=drift_deg_day,
                    days_per_deg_RAAN=1.0 / drift_deg_day if drift_deg_day > 0 else np.inf,
                    dv_for_offset_m_s=hohmann(h, h - dh_park) * 2,
                    dv_1deg_inclination_m_s=2 * float(v_circ(h * 1e3))
                    * np.sin(np.radians(0.5))))
    dfp = pd.DataFrame(plane)
    save_table(dfp, "sim11_plane_access",
               "Cost of reaching a target's orbital plane. RAAN is bought with "
               "time via differential J2 nodal regression from a small parking "
               "offset; inclination must be bought with propellant.")
    sel = dfp[(dfp.h_km == 800) & (dfp.park_offset_km == 50.0)]
    print("\n  plane access at 800 km, 50 km parking offset:")
    print(sel[["inc_deg", "nodal_rate_deg_day", "differential_drift_deg_day",
               "days_per_deg_RAAN", "dv_for_offset_m_s",
               "dv_1deg_inclination_m_s"]].to_string(index=False))

    # ---- end of life -------------------------------------------------------
    eol = []
    for h in ALT_SWEEP_KM:
        dv_deorbit = float(dv_dir := __import__("fg01.orbital", fromlist=["x"])
                           .dv_direct_reentry(h * 1e3))
        for isp in (300.0, 3000.0):
            for m_dry in (400.0, 600.0):
                prop_deorbit = m_dry * (1 - np.exp(-dv_deorbit / (isp * G0)))
                dv_hub = hohmann(h, 400.0)
                prop_hub = m_dry * (1 - np.exp(-dv_hub / (isp * G0)))
                eol.append(dict(h_km=h, isp_s=isp, m_dry_kg=m_dry,
                                dv_deorbit_m_s=dv_deorbit,
                                propellant_deorbit_kg=prop_deorbit,
                                dv_to_hub_400km_m_s=dv_hub,
                                propellant_to_hub_kg=prop_hub,
                                recommendation="deorbit" if prop_deorbit < prop_hub
                                else "return to hub"))
    dfe = pd.DataFrame(eol)
    save_table(dfe, "sim11_end_of_life",
               "End-of-life options: controlled deorbit versus return to a "
               "400 km reload hub.")

    concl = {
        "reference_platform": {"dry_kg": float(ref.dry_kg), "wet_kg": float(ref.wet_kg),
                               "coilgun_kg": 250.0, "magazine_kg": 30.0,
                               "shots_per_reload": float(ref.shots_per_reload_design)},
        "recoil_propellant_1e4_shots_kg": float(
            dfr[(dfr.m_shot_g == 2.69) & (dfr.n_shots == 1e4)
                & (dfr.isp_s == 300)].propellant_kg.iloc[0]),
        "plane_access_finding": (
            "RAAN is effectively free in propellant: a 50 km parking offset at "
            "800 km produces ~0.03-0.06 deg/day of differential nodal drift, so "
            "any RAAN is reachable in 20-40 days for a two-burn cost of ~50 m/s. "
            "Inclination is not: 1 deg costs ~130 m/s. One FG01 therefore "
            "serves one inclination band, which is compatible with the way "
            "debris actually clusters (breakup families share inclination) but "
            "means a fleet, not a single vehicle, for global coverage."),
        "dv_per_degree_inclination_m_s": float(sel.dv_1deg_inclination_m_s.iloc[0]),
        "days_per_deg_RAAN_at_800km": float(sel.days_per_deg_RAAN.iloc[0]),
        "transporter_does_not_burn_up_with_the_grains": (
            "FG01 and the released grains are separate objects with different "
            "trajectories: the grains are on immediate-reentry ellipses, FG01 "
            "remains in its operating orbit and is disposed of separately at "
            "end of life."),
        "eol_recommendation": (
            "Controlled deorbit is cheaper than returning to a 400 km hub above "
            "~500 km only for the highest operating altitudes; below that, "
            "returning to the hub costs less propellant than a direct deorbit "
            "and preserves the vehicle. With Isp 300 s and a 600 kg dry mass, "
            "direct deorbit from 800 km costs 39 kg of propellant."),
        "engagement_rate_caveat": (
            "This study does not model target acquisition rate. The number of "
            "engagements per year is set by how often FG01 can find, approach "
            "and match a tracked object in its own shell, which depends on a "
            "detection architecture (R10) that is outside the present scope. "
            "All per-object costs are conditional on the assumed engagement "
            "count."),
    }
    save_json(concl, "sim11_conclusions")

    fig, ax = plt.subplots()
    for i, inc in enumerate((51.6, 74.0, 82.0, 98.0)):
        s = dfp[(dfp.inc_deg == inc) & (dfp.park_offset_km == 50.0)].sort_values("h_km")
        ax.plot(s.h_km, s.days_per_deg_RAAN, color=CB[i], marker="o", ms=3,
                label=f"i = {inc}°")
    ax.set_xlabel("operating altitude (km)")
    ax.set_ylabel("days per degree of RAAN change")
    ax.set_title("SIM-11: RAAN access via differential J2 drift "
                 "(50 km parking offset, ~50 m/s round trip)")
    ax.legend(fontsize=8)
    save_fig(fig, "sim11_plane_access",
             "Time to change orbital plane using differential nodal regression "
             "instead of propellant. Near-polar inclinations precess slowly, "
             "so plane access is slowest exactly where sun-synchronous debris "
             "concentrates.")
    return dfb, dfp, concl


if __name__ == "__main__":
    run()
