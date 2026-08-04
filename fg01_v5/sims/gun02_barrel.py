"""
GUN-02 -- Barrel sizing for the FG01 chemical launcher at 619 m/s.

STANDALONE AND ADDITIVE. Imports gun01_ballistics read-only, writes only
gun02_* outputs, not registered in run_all.py.

Answers three questions for the 30 mm / 100 g / 619 m/s design point:
  1. barrel LENGTH  -- where do diminishing returns set in?
  2. barrel WALL    -- hoop stress, and the taper (pressure falls with travel,
                       so a constant wall is mostly dead mass)
  3. CHARGE         -- grams of propellant, and the chamber that holds it

The taper matters: peak pressure occurs within the first ~2 cm of travel, so
the chamber needs the full wall while the muzzle end needs almost nothing.
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import numpy as np
import pandas as pd

from gun01_ballistics import (charge_for_velocity, simulate, PROPELLANTS,
                              P_FRICTION_SMOOTH)
from fg01.io_utils import save_table, save_json, banner

V_TARGET = 619.0
BORE = 0.030
M_PROJ = 0.100
PROP = "Retumbo"

MATERIALS = {
    # yield Pa, density kg/m3, note
    "4340 steel":   (500e6, 7850, "baseline, cheap, tough"),
    "Ti-6Al-4V":    (900e6, 4430, "best metal strength/weight"),
    "Inconel 718":  (1030e6, 8190, "hot section, erosion resistant"),
    "AerMet 100":   (1720e6, 7890, "ultra-high-strength steel"),
}
SF = 2.0                    # burst safety factor on yield
T_MIN = 0.0015              # m, minimum practical wall (handling, threads)
T_LINER = 0.0010            # m, Inconel erosion liner at the chamber end


def run():
    banner("GUN-02  Barrel length, wall thickness and charge at 619 m/s")

    # ---- 1. LENGTH -------------------------------------------------------
    print("\n[1] BARREL LENGTH SWEEP (30 mm bore, 619 m/s, Retumbo)\n")
    rows = []
    for L in (0.5, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0):
        mc, r = charge_for_velocity(V_TARGET, d_bore=BORE, L_barrel=L,
                                    propellant=PROP)
        if not mc:
            continue
        rows.append(dict(L_m=L, charge_g=mc * 1e3, P_max_MPa=r["P_max_MPa"],
                         x_at_Pmax_cm=r["x_at_Pmax"] * 100,
                         a_max_g=r["a_max_g"], t_ms=r["t_muzzle_ms"],
                         eff_pct=r["efficiency_pct"],
                         recoil_Ns=r["recoil_impulse_Ns"],
                         under_cap=bool(r["P_max_MPa"] < 100.0)))
    d1 = pd.DataFrame(rows)
    d1["charge_saved_g"] = -d1.charge_g.diff()
    save_table(d1, "gun02_length",
               "Barrel length sweep at 619 m/s. Longer barrel buys a smaller "
               "charge and lower pressure, with diminishing returns.")
    print(d1.to_string(index=False, float_format=lambda x: f"{x:.2f}"))

    feas = d1[d1.under_cap]
    L_min = float(feas.L_m.min())
    knee = feas[feas.charge_saved_g < 1.0]
    L_opt = float(knee.L_m.iloc[0]) if len(knee) else 1.0
    print(f"\n    shortest barrel under 100 MPa: {L_min:.1f} m")
    print(f"    knee (charge saving < 1 g per step): {L_opt:.1f} m")

    # ---- 2. the design point and its pressure profile --------------------
    mc, r = charge_for_velocity(V_TARGET, d_bore=BORE, L_barrel=L_opt,
                                propellant=PROP)
    print(f"\n[2] DESIGN POINT: {BORE*1e3:.0f} mm bore, {L_opt:.1f} m barrel, "
          f"{mc*1e3:.1f} g {PROP}")
    for k in ("v_muzzle", "P_max_MPa", "x_at_Pmax", "t_muzzle_ms", "KE_J",
              "efficiency_pct", "a_max_g", "recoil_impulse_Ns", "V0_cm3"):
        print(f"    {k:24s} {r[k]:.4g}")

    # ---- 3. WALL THICKNESS, tapered -------------------------------------
    # envelope of pressure the wall at station x must contain: the max pressure
    # ever seen at that station, i.e. P(t) once the projectile has passed x.
    x_hist, P_hist = r["x"], r["P"]
    stations = np.linspace(0.0, L_opt, 60)
    P_env = np.array([P_hist[x_hist >= s].max() if (x_hist >= s).any()
                      else P_hist.max() for s in stations])
    # chamber (behind the projectile start) always sees full peak
    P_env = np.maximum.accumulate(P_env[::-1])[::-1]

    print(f"\n[4] WALL THICKNESS (hoop sigma = P*r/t, SF {SF:.0f} on yield)\n")
    rows = []
    for name, (sy, rho, note) in MATERIALS.items():
        t_ch = max(r["P_max_MPa"] * 1e6 * (BORE / 2) * SF / sy, T_MIN)
        t_prof = np.maximum(P_env * (BORE / 2) * SF / sy, T_MIN)
        # mass of the tapered barrel
        r_i = BORE / 2
        m_taper = np.trapezoid(rho * np.pi * ((r_i + t_prof) ** 2 - r_i ** 2),
                               stations)
        m_const = rho * np.pi * ((r_i + t_ch) ** 2 - r_i ** 2) * L_opt
        rows.append(dict(material=name, chamber_wall_mm=t_ch * 1e3,
                         muzzle_wall_mm=float(t_prof[-1]) * 1e3,
                         mass_tapered_kg=m_taper,
                         mass_constant_kg=m_const,
                         saving_pct=100 * (1 - m_taper / m_const),
                         note=note))
    d2 = pd.DataFrame(rows)
    save_table(d2, "gun02_wall",
               f"Barrel wall at SF {SF:.0f} on yield, chamber and muzzle ends, "
               "with tapered vs constant-wall mass.")
    print(d2.to_string(index=False, float_format=lambda x: f"{x:.3f}"))

    # taper profile for the chosen material
    pick = "Ti-6Al-4V"
    sy = MATERIALS[pick][0]
    t_prof = np.maximum(P_env * (BORE / 2) * SF / sy, T_MIN)
    prof = pd.DataFrame(dict(
        station_cm=stations * 100, P_envelope_MPa=P_env / 1e6,
        wall_mm=t_prof * 1e3, OD_mm=(BORE + 2 * t_prof) * 1e3))
    save_table(prof, "gun02_taper",
               f"Wall taper along the barrel for {pick} at SF {SF:.0f}.")
    print(f"\n    taper profile ({pick}):")
    print(prof.iloc[::8].to_string(index=False, float_format=lambda x: f"{x:.2f}"))

    # ---- 4. CHAMBER AND CARTRIDGE ---------------------------------------
    print(f"\n[5] CHAMBER AND CARTRIDGE\n")
    V0 = r["V0_cm3"] / 1e6
    l_ch = V0 / (np.pi * (BORE / 2) ** 2)
    rows = []
    for wall_mm, mat in ((1.5, "brass"), (1.2, "steel"), (0.0, "combustible")):
        if wall_mm > 0:
            rho_c = 8500 if mat == "brass" else 7850
            r_i, r_o = BORE / 2, BORE / 2 + wall_mm / 1e3
            v = np.pi * (r_o ** 2 - r_i ** 2) * l_ch + np.pi * r_o ** 2 * 0.008
            m = v * rho_c
        else:
            m = 0.25 * mc     # nitrocellulose-bonded case, ~25% of charge
        rows.append(dict(case_type=mat, wall_mm=wall_mm,
                         case_mass_g=m * 1e3,
                         round_mass_g=(mc + m) * 1e3 + M_PROJ * 1e3,
                         rounds_per_10kg=10.0 / (mc + m + M_PROJ)))
    d3 = pd.DataFrame(rows)
    save_table(d3, "gun02_cartridge",
               "Cartridge case options. Round mass includes the 100 g "
               "projectile.")
    print(f"    chamber volume {V0*1e6:.1f} cm3, length {l_ch*100:.1f} cm "
          f"(loading density 600 kg/m3)")
    print(d3.to_string(index=False, float_format=lambda x: f"{x:.2f}"))

    m_ti = float(d2[d2.material == "Ti-6Al-4V"].mass_tapered_kg.iloc[0])
    concl = {
        "design": dict(bore_mm=BORE * 1e3, length_m=L_opt,
                       charge_g=mc * 1e3, propellant=PROP,
                       P_max_MPa=r["P_max_MPa"],
                       v_muzzle=r["v_muzzle"],
                       chamber_volume_cm3=r["V0_cm3"],
                       chamber_length_cm=l_ch * 100),
        "wall": dict(
            material="Ti-6Al-4V",
            chamber_wall_mm=float(d2[d2.material == "Ti-6Al-4V"].chamber_wall_mm.iloc[0]),
            muzzle_wall_mm=float(d2[d2.material == "Ti-6Al-4V"].muzzle_wall_mm.iloc[0]),
            barrel_mass_kg=m_ti,
            taper_saving_pct=float(d2[d2.material == "Ti-6Al-4V"].saving_pct.iloc[0]),
            liner="1 mm Inconel 718 over the first ~15 cm for erosion"),
        "note_on_taper": (
            "Peak pressure occurs within the first few cm of travel, so the "
            "chamber end needs the full wall and the muzzle end reaches the "
            "1.5 mm practical minimum. Tapering saves a large fraction of the "
            "barrel mass at no performance cost."),
    }
    save_json(concl, "gun02_conclusions")
    return d1, d2, prof, concl


if __name__ == "__main__":
    run()
