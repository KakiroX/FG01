"""
LNCH-02 -- SIM-11 platform budget re-run with the CHEMICAL launcher.

STANDALONE AND ADDITIVE. sim11_platform.py is NOT modified; this file re-runs
the same platform model with the launcher architecture swapped, and writes only
lnch02_* outputs. Not registered in run_all.py.

SIM-11 was built around an electromagnetic launcher and therefore carries three
masses that exist ONLY to serve it:

    launcher hardware   100-500 kg   (parametric placeholder)
    power at 5 kW       166.7 kg     (5 kW is the capacitor CHARGING power --
                                      COIL-01 needs 95.3 kJ in 21 s = 4.5 kW)
    thermal             30 kg        (the coilgun must dissipate ~80% of its
                                      stored energy INSIDE the vehicle)

A chemical launcher removes all three: it needs no charging power, and it vents
most of its waste heat overboard with the propellant gas. What it adds instead
is a consumable -- propellant and case -- that grows linearly with shot count.

The question this file answers is therefore NOT the one LNCH-01 answered.
LNCH-01 compared the launchers alone and found a mass crossover at 545 shots.
Here the whole platform is compared, including the power and thermal plant each
architecture drags behind it.

Round mass accounting (identical convention to LNCH-01):
    the 100 g projectile is PAYLOAD in both architectures, so it cancels;
    the chemical consumable is propellant + case only.
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from fg01.constants import (R_E, G0, C_LAUNCH_RANGE, ALT_SWEEP_KM)
from fg01.orbital import dv_direct_reentry
from fg01.io_utils import save_table, save_fig, save_json, banner, CB

# ---- launcher design points (GUN-02, COIL-01/02/03, both at 619 m/s) --------
V_MUZZLE = 619.0
M_PROJ = 0.100                    # kg, on target -- payload in BOTH cases

# chemical, GUN-02
M_BARREL = 0.724                  # kg, Ti-6Al-4V tapered
M_BREECH = 6.70                   # kg, breech, mount, autoloader   UNVALIDATED
M_POWDER = 27.21e-3               # kg per round
CASES = {"brass": 139.1e-3, "steel": 111.0e-3, "combustible": 6.8e-3}
RECOIL_CHEM = 67.06               # N.s per shot, includes propellant gas
Q_BARREL = 23.7e3                 # J per shot into the barrel

# coilgun, COIL-01 at 619 m/s
M_CAPS = 63.55
M_CU = 11.90
K_INTEG = 1.30                    # switches, bus, mounts
M_COIL_LAUNCHER = K_INTEG * (M_CAPS + M_CU)
RECOIL_COIL = 62.34               # N.s per shot
E_STORED = 95.32e3                # J per shot
ETA_COIL = 0.2038
Q_COIL = E_STORED * (1 - ETA_COIL)   # J per shot to dissipate inside the vehicle

# ---- platform model, unchanged in form from SIM-11 --------------------------
M_AVIONICS = 50.0
W_PER_KG = 30.0                   # power plant specific mass
F_STRUCTURE = 0.20
P_COIL_W = 5000.0                 # SIM-11's value = the charging power
P_CHEM_W = 500.0                  # housekeeping, comms, ADCS only
M_THERMAL_COIL = 30.0             # SIM-11's value
# thermal plant scaled by waste heat actually rejected by the vehicle
M_THERMAL_CHEM = M_THERMAL_COIL * (Q_BARREL / Q_COIL)


def platform(launcher_kg, power_W, thermal_kg, consumable_kg, magazine_kg):
    power_kg = power_W / W_PER_KG
    subtotal = M_AVIONICS + power_kg + thermal_kg + launcher_kg + consumable_kg
    structure = F_STRUCTURE * (subtotal + magazine_kg)
    dry = subtotal + structure
    return dict(avionics_kg=M_AVIONICS, power_kg=power_kg, thermal_kg=thermal_kg,
                launcher_kg=launcher_kg, consumable_kg=consumable_kg,
                structure_kg=structure, dry_kg=dry, wet_kg=dry + magazine_kg)


def run():
    banner("LNCH-02  SIM-11 platform budget with the CHEMICAL launcher")

    n_ref = 48          # ECO-1 realised demand, shots per platform-decade
    case_ref = "combustible"

    # ---- 1. reference platform, both architectures ----------------------
    print("\n[1] REFERENCE PLATFORM  (48 shots, 30 kg rhenium magazine)\n")
    rows = []
    for case, m_case in CASES.items():
        p = platform(M_BARREL + M_BREECH, P_CHEM_W, M_THERMAL_CHEM,
                     n_ref * (M_POWDER + m_case), 30.0)
        p.update(architecture=f"chemical ({case} case)", round_g=(M_POWDER + m_case) * 1e3)
        rows.append(p)
    p = platform(M_COIL_LAUNCHER, P_COIL_W, M_THERMAL_COIL, 0.0, 30.0)
    p.update(architecture="coilgun", round_g=0.0)
    rows.append(p)
    df1 = pd.DataFrame(rows)[["architecture", "round_g", "launcher_kg", "power_kg",
                              "thermal_kg", "consumable_kg", "avionics_kg",
                              "structure_kg", "dry_kg", "wet_kg"]]
    df1["launch_cost_low_MUSD"] = df1.wet_kg * C_LAUNCH_RANGE[0] / 1e6
    df1["launch_cost_high_MUSD"] = df1.wet_kg * C_LAUNCH_RANGE[1] / 1e6
    save_table(df1, "lnch02_mass_budget",
               "FG01 platform mass budget with the chemical launcher, against "
               "the coilgun baseline SIM-11 assumed. 48 shots, 30 kg magazine. "
               "The 100 g projectile is payload in both cases and cancels.",
               tags={"launcher masses": "DERIVED (GUN-02, COIL-01)",
                     "breech/autoloader mass": "UNVALIDATED (parametric)",
                     "subsystem fractions": "UNVALIDATED (SIM-11 convention)"})
    print(df1.to_string(index=False, float_format=lambda x: f"{x:.2f}"))

    chem = df1[df1.architecture.str.contains(case_ref)].iloc[0]
    coil = df1[df1.architecture == "coilgun"].iloc[0]
    print(f"\n    chemical ({case_ref}) dry {chem.dry_kg:.1f} kg  vs  "
          f"coilgun dry {coil.dry_kg:.1f} kg   -> "
          f"{coil.dry_kg / chem.dry_kg:.2f}x")
    print(f"    saving {coil.dry_kg - chem.dry_kg:.1f} kg dry, "
          f"${(coil.wet_kg - chem.wet_kg) * C_LAUNCH_RANGE[1] / 1e6:.2f}M "
          f"at the high launch price")

    # ---- 2. where the saving comes from ---------------------------------
    print("\n[2] WHERE THE SAVING COMES FROM\n")
    src = pd.DataFrame([
        dict(item="launcher hardware", coilgun_kg=M_COIL_LAUNCHER,
             chemical_kg=M_BARREL + M_BREECH,
             why="63.6 kg of capacitors vs a 0.72 kg barrel"),
        dict(item="power plant", coilgun_kg=P_COIL_W / W_PER_KG,
             chemical_kg=P_CHEM_W / W_PER_KG,
             why="5 kW exists only to charge the bank in 21 s"),
        dict(item="thermal", coilgun_kg=M_THERMAL_COIL, chemical_kg=M_THERMAL_CHEM,
             why=f"{Q_COIL/1e3:.0f} kJ/shot rejected internally vs "
                 f"{Q_BARREL/1e3:.1f} kJ into the barrel"),
        dict(item="consumable (48 shots)", coilgun_kg=0.0,
             chemical_kg=n_ref * (M_POWDER + CASES[case_ref]),
             why="propellant + combustible case"),
    ])
    src["delta_kg"] = src.coilgun_kg - src.chemical_kg
    save_table(src, "lnch02_saving_breakdown",
               "Line-by-line source of the platform mass difference at 48 shots.")
    print(src.to_string(index=False, float_format=lambda x: f"{x:.2f}"))
    print(f"\n    subsystem delta {src.delta_kg.sum():.1f} kg, before the "
          f"{F_STRUCTURE:.0%} structure multiplier")

    # ---- 3. platform-level crossover ------------------------------------
    print("\n[3] PLATFORM-LEVEL CROSSOVER\n")
    n = np.arange(0, 12001)
    rows = []
    for case, m_case in CASES.items():
        m_ch = np.array([platform(M_BARREL + M_BREECH, P_CHEM_W, M_THERMAL_CHEM,
                                  k * (M_POWDER + m_case), 30.0)["dry_kg"]
                         for k in n])
        m_co = platform(M_COIL_LAUNCHER, P_COIL_W, M_THERMAL_COIL, 0.0,
                        30.0)["dry_kg"]
        over = np.where(m_ch > m_co)[0]
        n_x = int(over[0]) if len(over) else np.inf
        rows.append(dict(case_type=case, round_g=(M_POWDER + m_case) * 1e3,
                         crossover_shots=n_x,
                         margin_over_ECO1=n_x / n_ref if np.isfinite(n_x) else np.inf,
                         launcher_only_crossover=545))
    df3 = pd.DataFrame(rows)
    save_table(df3, "lnch02_crossover",
               "Shot count at which the chemical platform becomes heavier than "
               "the coilgun platform, compared with the launcher-only crossover "
               "of 545 shots found in LNCH-01.")
    print(df3.to_string(index=False, float_format=lambda x: f"{x:.2f}"))
    print(f"\n    ECO-1 realised demand is {n_ref} shots per platform-decade")

    # ---- 4. recoil ------------------------------------------------------
    print("\n[4] RECOIL AND STATION-KEEPING\n")
    rows = []
    for arch, imp in (("chemical", RECOIL_CHEM), ("coilgun", RECOIL_COIL)):
        for k in (48, 300, 1000):
            tot = k * imp
            for isp in (300.0, 3000.0):
                rows.append(dict(architecture=arch, impulse_per_shot_Ns=imp,
                                 n_shots=k, total_impulse_Ns=tot,
                                 dv_on_500kg_m_s=tot / 500.0,
                                 isp_s=isp, propellant_kg=tot / (isp * G0)))
    df4 = pd.DataFrame(rows)
    save_table(df4, "lnch02_recoil",
               "Recoil impulse and the propellant needed to hold the operating "
               "orbit. The chemical launcher's impulse includes the propellant "
               "gas, which the projectile-only figure (61.9 N.s) omits.")
    print(df4[df4.isp_s == 300.0].to_string(index=False,
                                            float_format=lambda x: f"{x:.3f}"))
    print(f"\n    chemical recoil is {100*(RECOIL_CHEM/RECOIL_COIL - 1):.1f}% "
          f"higher per shot; at {n_ref} shots that is "
          f"{n_ref*(RECOIL_CHEM-RECOIL_COIL)/(300*G0)*1e3:.1f} g extra propellant")

    # ---- 5. end of life, at the new dry mass ----------------------------
    print("\n[5] END OF LIFE AT THE NEW DRY MASS\n")
    rows = []
    for h in ALT_SWEEP_KM:
        dv = float(dv_direct_reentry(h * 1e3))
        for arch, m_dry in (("chemical", float(chem.dry_kg)),
                            ("coilgun", float(coil.dry_kg))):
            for isp in (300.0, 3000.0):
                rows.append(dict(architecture=arch, h_km=h, m_dry_kg=m_dry,
                                 isp_s=isp, dv_deorbit_m_s=dv,
                                 propellant_deorbit_kg=m_dry *
                                 (1 - np.exp(-dv / (isp * G0)))))
    df5 = pd.DataFrame(rows)
    save_table(df5, "lnch02_end_of_life",
               "Controlled-deorbit propellant at each architecture's dry mass. "
               "The lighter platform is cheaper to dispose of as well as to "
               "launch.")
    sel = df5[(df5.h_km == 800) & (df5.isp_s == 300.0)]
    print(sel.to_string(index=False, float_format=lambda x: f"{x:.2f}"))

    # ---- 6. figure ------------------------------------------------------
    fig, ax = plt.subplots()
    for i, (case, m_case) in enumerate(CASES.items()):
        m_ch = np.array([platform(M_BARREL + M_BREECH, P_CHEM_W, M_THERMAL_CHEM,
                                  k * (M_POWDER + m_case), 30.0)["dry_kg"]
                         for k in n])
        ax.plot(n, m_ch, color=CB[i], lw=1.6, label=f"chemical, {case} case")
    m_co = platform(M_COIL_LAUNCHER, P_COIL_W, M_THERMAL_COIL, 0.0, 30.0)["dry_kg"]
    ax.axhline(m_co, color=CB[3], lw=1.6, ls="--", label="coilgun (flat)")
    ax.axvline(n_ref, color="0.4", lw=1.0, ls=":")
    ax.annotate(f"ECO-1 demand\n{n_ref} shots", (n_ref, m_co * 0.45),
                xytext=(6, 0), textcoords="offset points", fontsize=8, color="0.3")
    ax.set_xlabel("shots per platform")
    ax.set_ylabel("platform dry mass (kg)")
    ax.set_xlim(0, 12000)
    ax.set_title("LNCH-02: platform dry mass vs shot count, chemical vs coilgun")
    ax.legend(fontsize=8, loc="upper left")
    save_fig(fig, "lnch02_crossover",
             "Platform dry mass against shot count. The coilgun is flat because "
             "it has no consumable; the chemical platform starts far lower and "
             "climbs. Including the power and thermal plant each architecture "
             "requires moves the crossover from 545 shots (launcher only) to "
             "several thousand.")

    concl = {
        "question": ("SIM-11 sized the FG01 platform around an electromagnetic "
                     "launcher. This re-run swaps in the chemical launcher "
                     "recommended by LNCH-01 and re-derives the platform."),
        "reference_platform_chemical": {
            "case": case_ref, "n_shots": n_ref,
            "dry_kg": float(chem.dry_kg), "wet_kg": float(chem.wet_kg),
            "launcher_kg": float(chem.launcher_kg),
            "power_kg": float(chem.power_kg), "thermal_kg": float(chem.thermal_kg)},
        "reference_platform_coilgun": {
            "dry_kg": float(coil.dry_kg), "wet_kg": float(coil.wet_kg),
            "launcher_kg": float(coil.launcher_kg),
            "power_kg": float(coil.power_kg), "thermal_kg": float(coil.thermal_kg)},
        "dry_mass_ratio": float(coil.dry_kg / chem.dry_kg),
        "headline": (
            "Swapping the coilgun for the chemical launcher removes three masses "
            "at once, not one: the launcher itself, the 5 kW power plant that "
            "existed only to charge the capacitor bank, and most of the thermal "
            "plant. The launcher hardware saving is the smallest of the three "
            "once the power plant is counted."),
        "crossover_finding": (
            "LNCH-01 put the LAUNCHER-only mass crossover at 545 shots. At "
            "platform level the crossover moves to " +
            ", ".join(f"{r.case_type} {int(r.crossover_shots)}"
                      for r in df3.itertuples()) +
            " shots, because the coilgun platform also carries ~150 kg of power "
            "and thermal plant the chemical platform does not. Against ECO-1's "
            "48 shots per platform-decade the margin is one to two orders of "
            "magnitude."),
        "recoil_finding": (
            f"Chemical recoil is {RECOIL_CHEM:.1f} N.s per shot against the "
            f"coilgun's {RECOIL_COIL:.1f}, {100*(RECOIL_CHEM/RECOIL_COIL-1):.0f}% "
            "higher because the propellant gas carries momentum the projectile "
            "does not. At 48 shots the difference is a few grams of "
            "station-keeping propellant -- immaterial."),
        "unchanged_by_this_swap": (
            "Plane access and the RAAN-versus-inclination finding in SIM-11 are "
            "properties of the orbit, not the launcher, and are unaffected. "
            "SIM-11 remains the authority on them."),
        "caveats": [
            "Breech, mount and autoloader mass (6.7 kg) is a parametric estimate, "
            "UNVALIDATED. It is the largest single uncertainty in the chemical "
            "column, but it would have to be wrong by a factor of ~20 to change "
            "the conclusion.",
            "Thermal plant is scaled linearly by rejected heat from SIM-11's "
            "30 kg. Radiator mass is not in general linear in heat load; this is "
            "a first-order treatment, UNVALIDATED.",
            "The 500 W housekeeping figure for the chemical platform is assumed, "
            "not derived. A detection architecture (R10) could raise it, but it "
            "would raise it for both architectures equally.",
            "Live propellant on an unattended long-duration platform is a "
            "programmatic and safety question this model does not evaluate. It "
            "remains the most likely single reason to choose the coilgun anyway.",
        ],
    }
    save_json(concl, "lnch02_conclusions")

    print("\n" + "=" * 70)
    print(f"  chemical platform  {chem.dry_kg:6.1f} kg dry")
    print(f"  coilgun platform   {coil.dry_kg:6.1f} kg dry   "
          f"({coil.dry_kg/chem.dry_kg:.2f}x)")
    print(f"  platform crossover {int(df3[df3.case_type==case_ref].crossover_shots.iloc[0])} "
          f"shots (launcher-only was 545)")
    print("=" * 70)
    return df1, df3, concl


run = run

if __name__ == "__main__":
    run()
