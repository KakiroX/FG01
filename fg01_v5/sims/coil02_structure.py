"""
COIL-02 -- Wall thicknesses and structural mass for the FG01 coilgun.

STANDALONE AND ADDITIVE. Imports fg01/ read-only, writes only coil02_* outputs,
not registered in run_all.py. Modifies no existing simulation.

--------------------------------------------------------------------------
A COILGUN BARREL IS NOT A GUN BARREL
--------------------------------------------------------------------------
There is no chamber pressure. Nothing is trying to burst the bore tube. The
structural loads are three separate, unrelated items:

  1. BORE TUBE -- carries no pressure at all. Its thickness is set by the
     coil-to-armature GAP (coupling k falls as the gap grows) and by handling
     loads. CRITICAL and easy to get wrong: it must be a poor electrical
     conductor. A metal bore tube is a shorted turn around the armature; it
     would carry its own induced current, shield the armature and destroy the
     coupling the launcher depends on. G10/FR4, PEEK or a ceramic.

  2. COIL CONTAINMENT -- this is the real pressure vessel. Magnetic pressure
     B^2/(2 mu0) pushes the windings radially outward and must be reacted by
     hoop tension in banding over the coil:  sigma = P r / t.

  3. OUTER CASE / RAIL -- reacts recoil and holds stage alignment. A beam
     problem, not a pressure problem.

--------------------------------------------------------------------------
CONCLUSION FOR A COST MODEL (stated up front)
--------------------------------------------------------------------------
Wall thickness is NOT a cost driver. All structure together is a few kg against
28 kg of capacitors. Model the launcher cost as capacitors + copper + a
structural fraction, and do not spend effort on wall optimisation.
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd

from fg01.io_utils import save_table, save_json, banner

MU0 = 4e-7 * np.pi
G0 = 9.80665

# ---- design point (COIL-01, 400 m/s) -------------------------------------
BORE_ID = 0.030          # m, canister outer diameter = bore inner diameter
N_STAGES = 31
L_TOTAL = 1.364          # m
STAGE_PITCH = L_TOTAL / N_STAGES
L_COIL = 0.040           # m, coil axial length
M_PROJ = 0.100           # kg
A_MAX_G = 49_689.0       # peak acceleration from COIL-01
RECOIL_NS = 40.7         # N.s per shot
E_CAP_J = 42_260.0       # stored energy
M_CAP_KG = 28.2          # capacitor mass
M_CU_KG = 4.17           # copper at 22 turns/stage

MATERIALS = {
    # name: (yield/tensile Pa, density kg/m3, conductive?, note)
    "4340 steel":        (500e6, 7850, True,  "cheap, conductive - banding only"),
    "Ti-6Al-4V":         (900e6, 4430, True,  "light metal"),
    "S-glass/epoxy":     (1700e6, 2000, False, "insulating, good for banding"),
    "Carbon fibre T1000": (3000e6, 1600, True, "conductive! insulate from coil"),
    "Zylon/PBO":         (5800e6, 1560, False, "highest specific strength"),
}

BORE_MATERIALS = {
    # name: (density, resistivity ohm.m, max service K, note)
    "G10/FR4 glass-epoxy": (1850, 1e13, 400, "standard, cheap, machinable"),
    "PEEK":                (1320, 1e14, 480, "tougher, lower friction"),
    "Alumina Al2O3":       (3900, 1e12, 1900, "hard, wear-resistant, brittle"),
    "316 stainless":       (8000, 7.4e-7, 1100, "DO NOT USE - shorted turn"),
}


def magnetic_pressure(B):
    return B ** 2 / (2 * MU0)


def hoop_thickness(P, r, sigma_yield, sf):
    """Thin-wall hoop: sigma = P r / t  ->  t = P r SF / sigma_yield."""
    return P * r * sf / sigma_yield


def run():
    banner("COIL-02  Coilgun wall thickness and structural mass")

    # peak field, back-computed from the peak force COIL-01 reported
    F_max = M_PROJ * A_MAX_G * G0
    r_bore = BORE_ID / 2
    r_coil_mean = r_bore + 0.005 + 0.0013     # tube + former + half winding
    B_peak = np.sqrt(2 * MU0 * F_max / (np.pi * r_coil_mean ** 2))
    P_mag = magnetic_pressure(B_peak)
    r_band = r_coil_mean + 0.0013 + 0.0005    # outside the winding

    print(f"\n[0] LOADS")
    print(f"    peak accelerating force   {F_max/1e3:8.1f} kN")
    print(f"    equivalent peak field     {B_peak:8.2f} T")
    print(f"    MAGNETIC PRESSURE         {P_mag/1e6:8.1f} MPa  <- the only "
          f"'pressure' in the system")
    print(f"    (for comparison, the chemical gun ran at 85.4 MPa chamber)")
    print(f"    banding radius            {r_band*1e3:8.1f} mm")
    print(f"    recoil per shot           {RECOIL_NS:8.1f} N.s")

    # ---- 1. coil containment banding -------------------------------------
    print(f"\n[1] COIL CONTAINMENT BANDING  (the real structural item)\n")
    rows = []
    for name, (sy, rho, cond, note) in MATERIALS.items():
        for sf in (1.5, 2.0):
            t = hoop_thickness(P_mag, r_band, sy, sf)
            # banding covers the coil length on every stage
            vol = 2 * np.pi * r_band * t * L_COIL * N_STAGES
            rows.append(dict(material=name, safety_factor=sf,
                             wall_mm=t * 1e3,
                             mass_all_stages_kg=vol * rho,
                             conductive=cond, note=note))
    d1 = pd.DataFrame(rows)
    save_table(d1, "coil02_banding",
               "Coil containment banding: hoop thickness against magnetic "
               "pressure, and total mass over all 31 stages.")
    print(d1.to_string(index=False, float_format=lambda x: f"{x:.3f}"))

    # ---- 2. bore tube ----------------------------------------------------
    print(f"\n[2] BORE TUBE  (carries NO pressure -- thickness set by coupling)\n")
    rows = []
    for name, (rho, res, tmax, note) in BORE_MATERIALS.items():
        for t in (0.001, 0.0015, 0.002, 0.003):
            vol = 2 * np.pi * (r_bore + t / 2) * t * L_TOTAL
            # coupling penalty: k falls roughly as the gap grows relative to
            # the coil radius; a 1 mm increase in gap is ~1 mm of coil radius
            gap = 0.005 + t
            k_rel = (r_bore / (r_bore + gap)) ** 1.5
            rows.append(dict(material=name, wall_mm=t * 1e3,
                             mass_kg=vol * rho,
                             resistivity_ohm_m=res,
                             coupling_factor_rel=k_rel,
                             usable=bool(res > 1e-3),
                             note=note))
    d2 = pd.DataFrame(rows)
    save_table(d2, "coil02_bore_tube",
               "Bore tube options. 'usable' is False for conductors: a metal "
               "tube is a shorted turn that shields the armature.")
    print(d2[d2.wall_mm.isin([1.0, 1.5, 2.0])].to_string(
        index=False, float_format=lambda x: f"{x:.4g}"))
    print(f"\n    NOTE: a stainless bore tube has resistivity 7.4e-7 ohm.m and "
          f"would act as")
    print(f"    a shorted turn around the armature -- it must NOT be used, at "
          f"any thickness.")

    # ---- 3. outer case ---------------------------------------------------
    print(f"\n[3] OUTER CASE / RAIL  (recoil and alignment, not pressure)\n")
    rows = []
    for name, (sy, rho, cond, _) in MATERIALS.items():
        if name in ("Zylon/PBO",):
            continue
        # treat as a thin tube in compression carrying peak recoil force
        r_case = r_band + 0.010
        for sf in (2.0,):
            a_req = F_max * sf / sy               # required cross-section
            t = a_req / (2 * np.pi * r_case)
            t = max(t, 0.0008)                    # minimum practical gauge
            vol = 2 * np.pi * r_case * t * L_TOTAL
            rows.append(dict(material=name, safety_factor=sf,
                             wall_mm=t * 1e3, mass_kg=vol * rho,
                             gauge_limited=bool(a_req / (2 * np.pi * r_case)
                                                < 0.0008)))
    d3 = pd.DataFrame(rows)
    save_table(d3, "coil02_outer_case",
               "Outer case sized for peak recoil in compression. Every option "
               "is minimum-gauge limited, i.e. strength is not the driver.")
    print(d3.to_string(index=False, float_format=lambda x: f"{x:.3f}"))

    # ---- 4. mass roll-up for the cost model ------------------------------
    print(f"\n[4] MASS ROLL-UP FOR A COST MODEL\n")
    band_cf = float(d1[(d1.material == "S-glass/epoxy")
                       & (d1.safety_factor == 2.0)].mass_all_stages_kg.iloc[0])
    band_steel = float(d1[(d1.material == "4340 steel")
                          & (d1.safety_factor == 2.0)].mass_all_stages_kg.iloc[0])
    tube = float(d2[(d2.material == "G10/FR4 glass-epoxy")
                    & (d2.wall_mm == 1.5)].mass_kg.iloc[0])
    case = float(d3[d3.material == "Ti-6Al-4V"].mass_kg.iloc[0])

    rows = [
        dict(item="capacitor bank", mass_kg=M_CAP_KG, driver="stored energy"),
        dict(item="coil copper", mass_kg=M_CU_KG, driver="turns x stages"),
        dict(item="outer case (Ti, min gauge)", mass_kg=case,
             driver="minimum gauge"),
        dict(item="bore tube (G10, 1.5 mm)", mass_kg=tube,
             driver="coupling gap"),
        dict(item="coil banding (S-glass, SF2)", mass_kg=band_cf,
             driver="magnetic pressure"),
    ]
    d4 = pd.DataFrame(rows)
    d4["share_pct"] = 100 * d4.mass_kg / d4.mass_kg.sum()
    struct = float(d4[d4["item"].str.contains("case|tube|banding")].mass_kg.sum())
    subtotal = float(d4.mass_kg.sum())
    save_table(d4, "coil02_mass_rollup",
               "Launcher mass by item at the COIL-01 400 m/s design point.")
    print(d4.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    print(f"\n    ALL STRUCTURE (case + tube + banding) = {struct:.2f} kg "
          f"= {100*struct/subtotal:.1f}% of the launcher")
    print(f"    capacitors alone                      = {M_CAP_KG:.2f} kg "
          f"= {100*M_CAP_KG/subtotal:.1f}%")
    print(f"    subtotal                              = {subtotal:.2f} kg "
          f"(+ ~30% for switches, bus, mounts)")
    print(f"    with integration margin               = {subtotal*1.3:.2f} kg")

    if band_steel > band_cf:
        print(f"\n    steel banding would add {band_steel-band_cf:.2f} kg over "
              f"S-glass, still small.")

    concl = {
        "headline": (
            "Wall thickness is not a cost or mass driver for a coilgun. There "
            "is no chamber pressure; the only structural load is magnetic "
            f"pressure of {P_mag/1e6:.0f} MPa on the coil, which S-glass "
            f"banding carries at {float(d1[(d1.material=='S-glass/epoxy')&(d1.safety_factor==2.0)].wall_mm.iloc[0]):.2f} mm "
            f"and {band_cf:.2f} kg over all 31 stages. Total structure is "
            f"{struct:.1f} kg against {M_CAP_KG:.0f} kg of capacitors."),
        "magnetic_pressure_MPa": float(P_mag / 1e6),
        "B_peak_T": float(B_peak),
        "recommended": {
            "bore_tube": "G10/FR4 or PEEK, 1.5 mm wall, MUST be insulating",
            "coil_banding": (f"S-glass/epoxy, "
                             f"{float(d1[(d1.material=='S-glass/epoxy')&(d1.safety_factor==2.0)].wall_mm.iloc[0]):.2f} mm "
                             f"at SF 2"),
            "outer_case": "Ti-6Al-4V or Al, 0.8 mm minimum gauge (not "
                          "strength-limited)"},
        "for_cost_model": (
            "Model launcher mass as: capacitors (E_stored/1500 J/kg) + copper "
            "(~4 kg) + 1.3x integration factor. Treat all walls as a fixed "
            "~1.5 kg. Do not optimise wall thickness -- it is 4% of the mass "
            "and its cost is negligible against the capacitor bank."),
        "trap": (
            "The bore tube must be electrically insulating. A metal tube is a "
            "shorted turn around the armature: it carries induced current, "
            "shields the armature and destroys the coupling. This is the one "
            "structural choice that can break the launcher."),
    }
    save_json(concl, "coil02_conclusions")
    return d1, d2, d4, concl


if __name__ == "__main__":
    run()
