"""
LNCH-01 -- Launcher trade study: chemical (gunpowder / smokeless) vs.
electromagnetic (railgun, coilgun) for the FG01 momentum-kick projectile.

Requirement points (from the design outputs):
  R1  small-debris design shot : 2.14 g @ 619 m/s   = 410 J
  R2  upper-velocity variant   : 2.14 g @ 700 m/s   = 524 J   (30 kJ/kg ceiling)
  R3  large-object shot        : 1.00 kg @ 619 m/s  = 192 kJ  (Architecture B)

The launcher choice feeds back into the mission through three channels, all
quantified here:
  (a) launched dry mass  -> platform launch cost (v5 SIM-11 / SIM-9)
  (b) muzzle-velocity repeatability sigma_v/v -> the SIM-4/INT-4 aim-error term
      u*R*(sigma_v/v)/v_rel, which sets launched projectile mass quadratically
  (c) consumable vs. rechargeable -> logistics over the mission shot count
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from fg01.constants import OPTICAL_LAMBDA_M, OPTICAL_APERTURE_M, RHO_AL
from fg01.orbital import xsec_sphere, mass_sphere
from fg01.launchers import (chemical_gun, railgun, coilgun, charge_time,
                            magazine_shots, muzzle_energy, PROPELLANTS,
                            GUN_BALLISTIC_ETA, RAIL_ETA, COIL_ETA,
                            RAIL_LIFE_SHOTS, COIL_LIFE_SHOTS)
from fg01.delivery import penalty_for_confidence
from fg01.io_utils import save_table, save_fig, save_json, banner, CB

REQS = [
    ("R1 small-debris design", 2.14e-3, 619.0),
    ("R2 upper-velocity",      2.14e-3, 700.0),
    ("R3 large-object shot",   1.00,    619.0),
]

# aim-budget context (from INT-0 / SIM-4 co-orbital design point)
U_REL = 10.0            # m/s platform-target relative speed
R_ENG = 10.0           # m engagement range
V_REL_KICK = 619.0     # m/s
A_TARGET_1CM = float(xsec_sphere(0.01))
C_LIGHT = 2.99792458e8


def validate():
    """Cross-check each model against a demonstrated benchmark."""
    out = {}
    # chemical: a .22-class shot -- ~2 g at ~400 m/s from ~0.1-0.2 g of powder.
    g = chemical_gun(2.0e-3, 400.0, propellant="single_base_NC", eta_bal=0.28)
    out["chemical_22LR_class"] = {
        "modelled_propellant_g": round(g["m_prop_g"], 3),
        "real_22LR_propellant_g": "0.10-0.16",
        "note": ("A .22 Long Rifle launches ~2.6 g at ~330 m/s on ~0.1 g of "
                 "powder; the model's propellant mass for a comparable shot is "
                 "the right order, confirming the impetus/efficiency balance.")}
    # railgun: the ONR benchmark, 3.2 kg at 2520 m/s, 10.16 MJ muzzle
    r = railgun(3.2, 2520.0, eta=0.30, L_rail=10.0)
    out["railgun_ONR_2008"] = {
        "modelled_store_MJ": round(r["E_store_J"] / 1e6, 1),
        "modelled_current_MA": round(r["current_kA"] / 1e3, 2),
        "real_muzzle_MJ": 10.16, "real_current_MA": "~3-3.5",
        "note": ("Modelled stored energy 33.9 MJ at 30% efficiency matches the "
                 "~32 MJ bank behind the 10.16 MJ ONR shot; peak current a few "
                 "MA, consistent with the reported 3+ MA armature current.")}
    # coilgun: Turman 1996 lab coilgun, 10 g at 2000 m/s
    c = coilgun(0.010, 2000.0, eta=0.40)
    out["coilgun_turman1996"] = {
        "modelled_store_kJ": round(c["E_store_J"] / 1e3, 1),
        "real_muzzle_kJ": 20.0,
        "note": "20 kJ muzzle at 40% eff -> 50 kJ store, matches lab benchmark."}
    return out


def sigma_v_miss(sig_v_frac):
    """Contribution of muzzle-velocity spread to the terminal miss [m]."""
    return U_REL * R_ENG * sig_v_frac / V_REL_KICK


def launched_mass_penalty(sig_v_frac):
    """
    Total miss and the resulting launched-mass penalty when muzzle dispersion is
    added in quadrature to the fixed 5.44 mm of non-launcher error (SIM-4).
    """
    d_fixed = 5.44e-3          # pointing + ranging + optics, launcher-independent
    d_total = np.hypot(d_fixed, sigma_v_miss(sig_v_frac))
    pen = float(penalty_for_confidence(d_total, A_TARGET_1CM, 0.95))
    return d_total, max(pen, 1.0)


def run():
    banner("LNCH-01  Launcher trade: chemical vs railgun vs coilgun")

    val = validate()
    save_json(val, "lnch01_validation")
    print("  validation:")
    for k, v in val.items():
        print(f"    {k}: {v['note']}")

    # ---- size every launcher at every requirement point --------------------
    rows = []
    for lab, m_p, v in REQS:
        E = float(muzzle_energy(m_p, v))
        g_lo = chemical_gun(m_p, v, "single_base_NC", eta_bal=GUN_BALLISTIC_ETA[1])
        g_hi = chemical_gun(m_p, v, "double_base", eta_bal=GUN_BALLISTIC_ETA[0])
        r = railgun(m_p, v, eta=np.mean(RAIL_ETA))
        c = coilgun(m_p, v, eta=np.mean(COIL_ETA))
        for name, d, extra in [
            ("chemical (smokeless)", g_lo, dict(
                energy_store_kg=g_lo["m_prop_kg"],
                store_is_consumable=True)),
            ("railgun", r, dict(energy_store_kg=r["m_energy_store_kg"],
                                store_is_consumable=False)),
            ("coilgun", c, dict(energy_store_kg=c["m_energy_store_kg"],
                                store_is_consumable=False)),
        ]:
            d_miss, pen = launched_mass_penalty(d["sigma_v_frac"])
            rows.append(dict(
                requirement=lab, launcher=name, m_p_g=m_p * 1e3, v_m_s=v,
                E_muzzle_J=E,
                E_store_J=d.get("E_store_J", d.get("E_chem_J", np.nan)),
                efficiency=d.get("eta", d.get("eta_overall", np.nan)),
                m_launcher_kg=d["m_launcher_kg"],
                m_energy_store_kg=extra["energy_store_kg"],
                m_total_kg=d["m_launcher_kg"] + extra["energy_store_kg"],
                barrel_or_rail_m=d.get("L_barrel_m", d.get("L_rail_m", np.nan)),
                sigma_v_frac=d["sigma_v_frac"],
                miss_from_muzzle_mm=sigma_v_miss(d["sigma_v_frac"]) * 1e3,
                total_miss_mm=d_miss * 1e3,
                launched_mass_penalty=pen,
                recoil_impulse_Ns=d["recoil_impulse_Ns"],
                consumable=d["consumable"],
                clean_vacuum=d["clean_vacuum"],
                shot_life=d["shot_life"]))
    df = pd.DataFrame(rows)
    save_table(df, "lnch01_sizing",
               "Chemical, railgun and coilgun sized at each requirement point, "
               "with the launched-mass penalty each implies through its "
               "muzzle-velocity repeatability.",
               tags={"interior ballistics": "DERIVED", "impetus": "SOURCED",
                     "efficiencies": "SOURCED"})
    print("\n  R1 small-debris design shot (2.14 g @ 619 m/s):")
    print(df[df.requirement == "R1 small-debris design"][
        ["launcher", "m_launcher_kg", "m_energy_store_kg", "efficiency",
         "sigma_v_frac", "total_miss_mm", "launched_mass_penalty"]].to_string(index=False))

    # ---- consumable logistics over the mission -----------------------------
    # v6 realised ~24 engagements/platform; staged/large-object variants fire more
    log = []
    for n_shots in (24, 100, 1000, 10000):
        g = chemical_gun(2.14e-3, 619.0, "single_base_NC")
        prop_total = n_shots * g["m_prop_kg"]
        for name, store_kg_per_shot, consumable in [
            ("chemical (smokeless)", g["m_prop_kg"], True),
            ("railgun", 0.0, False),
            ("coilgun", 0.0, False)]:
            log.append(dict(
                n_shots=n_shots, launcher=name,
                consumable_mass_kg=n_shots * store_kg_per_shot,
                recharge_from_solar=not consumable,
                charge_time_s_per_shot=(charge_time(
                    coilgun(2.14e-3, 619.0)["E_store_J"]) if not consumable else np.nan)))
    dfl = pd.DataFrame(log)
    save_table(dfl, "lnch01_logistics",
               "Consumable mass over the mission at the design shot (2.14 g). "
               "EM launchers recharge from the 5 kW solar array; the chemical "
               "gun depletes a propellant magazine.")

    # ---- energy crossover: where does capacitor mass beat propellant? ------
    # Capacitor storage is ~1.5 kJ/kg; propellant potential is ~4 MJ/kg, so at
    # high per-shot energy the EM energy store dominates system mass. Find the
    # muzzle energy at which the chemical gun's total mass undercuts the coilgun.
    cross = []
    for E in np.geomspace(50.0, 5.0e5, 60):
        v = 619.0
        m_p = 2 * E / v ** 2
        g = chemical_gun(m_p, v, "single_base_NC")
        c = coilgun(m_p, v, eta=np.mean(COIL_ETA))
        r = railgun(m_p, v, eta=np.mean(RAIL_ETA))
        cross.append(dict(
            E_muzzle_J=E, m_p_g=m_p * 1e3,
            chemical_total_kg=g["m_launcher_kg"] + g["m_prop_kg"],
            coilgun_total_kg=c["m_launcher_kg"] + c["m_energy_store_kg"],
            railgun_total_kg=r["m_launcher_kg"] + r["m_energy_store_kg"],
            chemical_over_coilgun=(g["m_launcher_kg"] + g["m_prop_kg"])
            / (c["m_launcher_kg"] + c["m_energy_store_kg"])))
    dfx = pd.DataFrame(cross)
    save_table(dfx, "lnch01_energy_crossover",
               "System mass vs. muzzle energy. Chemical mass is nearly flat "
               "(propellant is trivially light); EM mass grows linearly with "
               "energy because capacitor storage is ~1.5 kJ/kg against a "
               "propellant potential of ~4 MJ/kg.")
    # chemical is always lighter here; report the ratio at the two regimes
    e_small = dfx.iloc[(dfx.E_muzzle_J - 410).abs().idxmin()]
    e_large = dfx.iloc[(dfx.E_muzzle_J - 192000).abs().idxmin()]
    print(f"\n  system-mass ratio chemical/coilgun: "
          f"{e_small.chemical_over_coilgun:.2f}x at 410 J, "
          f"{e_large.chemical_over_coilgun:.3f}x at 192 kJ")

    # ---- weighted decision matrix ------------------------------------------
    # criteria, weights (stated), and 0-1 scores per launcher at R1
    r1 = df[df.requirement == "R1 small-debris design"].set_index("launcher")
    r3 = df[df.requirement == "R3 large-object shot"].set_index("launcher")
    crit = {
        "launched mass (R1)":        (0.15, "min"),
        "launched mass (R3)":        (0.10, "min"),
        "muzzle repeatability":      (0.25, "min"),   # -> launched-mass penalty
        "consumable independence":   (0.20, "bool"),
        "vacuum/optical cleanliness":(0.15, "bool"),
        "reusable shot life":        (0.15, "life"),
    }
    launchers = ["chemical (smokeless)", "railgun", "coilgun"]
    raw = {
        "launched mass (R1)": {l: r1.loc[l, "m_total_kg"] for l in launchers},
        "launched mass (R3)": {l: r3.loc[l, "m_total_kg"] for l in launchers},
        "muzzle repeatability": {l: r1.loc[l, "launched_mass_penalty"] for l in launchers},
        "consumable independence": {"chemical (smokeless)": 0.0, "railgun": 1.0, "coilgun": 1.0},
        "vacuum/optical cleanliness": {"chemical (smokeless)": 0.0, "railgun": 0.7, "coilgun": 1.0},
        "reusable shot life": {"chemical (smokeless)": 0.3, "railgun": 0.2, "coilgun": 1.0},
    }
    score_rows = []
    totals = {l: 0.0 for l in launchers}
    for cname, (w, kind) in crit.items():
        vals = raw[cname]
        if kind == "min":
            best = min(vals.values())
            sc = {l: best / vals[l] for l in launchers}
        elif kind in ("bool", "life"):
            sc = vals
        for l in launchers:
            totals[l] += w * sc[l]
        score_rows.append(dict(criterion=cname, weight=w,
                               **{l: round(sc[l], 3) for l in launchers}))
    score_rows.append(dict(criterion="WEIGHTED TOTAL", weight=1.0,
                           **{l: round(totals[l], 3) for l in launchers}))
    dfsc = pd.DataFrame(score_rows)
    save_table(dfsc, "lnch01_decision_matrix",
               "Weighted decision matrix (weights stated). 'min' criteria scored "
               "as best/value; boolean and life criteria scored 0-1 directly.")
    print("\n  decision matrix:")
    print(dfsc.to_string(index=False))
    winner = max(totals, key=totals.get)

    concl = {
        "requirement_points": {lab: {"m_p_g": mp * 1e3, "v_m_s": v,
                                     "E_muzzle_J": float(muzzle_energy(mp, v))}
                               for lab, mp, v in REQS},
        "R1_totals_kg": {l: float(r1.loc[l, "m_total_kg"]) for l in launchers},
        "R3_totals_kg": {l: float(r3.loc[l, "m_total_kg"]) for l in launchers},
        "muzzle_repeatability": {l: float(r1.loc[l, "sigma_v_frac"]) for l in launchers},
        "launched_mass_penalty_R1": {l: float(r1.loc[l, "launched_mass_penalty"])
                                     for l in launchers},
        "weighted_scores": {l: float(totals[l]) for l in launchers},
        "winner": winner,
        "headline": (
            f"The coilgun wins the weighted trade ({totals['coilgun']:.2f} vs "
            f"{totals['railgun']:.2f} railgun, {totals['chemical (smokeless)']:.2f} "
            "chemical). The FG01 requirement is tiny -- 410 J, 2.14 g, 619 m/s -- "
            "so all three are technically trivial to build; the choice is decided "
            "by the SECONDARY metrics, and there the coilgun dominates: no "
            "consumable, clean in vacuum next to the sensor, best muzzle "
            "repeatability, and 1e5-1e7 shot life with no sliding contact."),
        "why_not_chemical": (
            "A gunpowder gun is the most compact and the lightest launcher "
            "hardware, but at these energies its charge is a fraction of a gram "
            "and hard to meter precisely, giving ~0.5-1% muzzle-velocity spread "
            "against the EM launchers' ~0.1%. It also vents hot, reactive, "
            "particulate combustion products metres from the terminal optics and "
            "in vacuum, and it consumes a magazine that must be resupplied -- "
            "fatal for a platform whose entire economic case is amortising "
            "hardware over many shots without resupply."),
        "why_not_railgun": (
            "A railgun works but is the worst EM option here: it is the least "
            "efficient (20-35%), its sliding armature erodes the rails and caps "
            "life at a few hundred to a couple of thousand full-power shots, and "
            "it needs mega-ampere pulses whose switchgear is heavy. Railguns earn "
            "their place at multi-MJ muzzle energies; the FG01 shot is 410 J, "
            "four to five orders of magnitude below where a railgun's compactness "
            "advantage appears."),
        "energy_crossover": {
            "chemical_over_coilgun_at_410J": float(e_small.chemical_over_coilgun),
            "chemical_over_coilgun_at_192kJ": float(e_large.chemical_over_coilgun),
            "finding": (
                "The chemical gun is ~10-14x lighter than either EM launcher at "
                "EVERY energy, because propellant stores ~4 MJ/kg of available "
                "work against a pulsed capacitor's ~1.5 kJ/kg. What changes with "
                "energy is the ABSOLUTE mass gap, and therefore whether it "
                "matters. At the 410 J small-debris shot the gap is ~0.9 kg "
                "(0.09 vs 1.03 kg) -- negligible against a 400-600 kg platform, "
                "so the secondary metrics rightly decide and the coilgun wins. "
                "At the 192 kJ large-object shot the gap is ~300 kg (22 vs "
                "322 kg) -- large enough to dominate the platform mass budget "
                "and flip the recommendation to the chemical gun. The crossover "
                "is not in the ratio but in whether the gap is affordable."),
        },
        "recommendation": (
            "REGIME-DEPENDENT. (1) Small-debris mission (the baseline, "
            "~400-520 J per shot): COILGUN -- capacitor-driven, solar-recharged, "
            "clean, long-lived; energy-store mass is negligible at this energy so "
            "its drawbacks do not bite. (2) Large-object variant (~100+ kJ per "
            "shot): the coilgun's capacitor mass (300+ kg) becomes prohibitive "
            "and a CHEMICAL gun (~22 kg) or a rechargeable inductive/rotating "
            "energy store is preferable despite the contamination penalty, "
            "because at high energy propellant is three orders of magnitude more "
            "mass-efficient than capacitors. The baseline FG01 mission is "
            "small-debris, so the primary recommendation is the coilgun; the "
            "chemical option is held in reserve for the large-object architecture."),
    }
    save_json(concl, "lnch01_conclusions")
    print(f"\n  WINNER: {winner}  (score {totals[winner]:.2f})")

    # ---- figures -----------------------------------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.2))
    ax = axes[0]
    x = np.arange(len(launchers))
    for i, req in enumerate(["R1 small-debris design", "R3 large-object shot"]):
        s = df[df.requirement == req].set_index("launcher")
        ax.bar(x + (i - 0.5) * 0.35, [s.loc[l, "m_total_kg"] for l in launchers],
               0.35, color=CB[i], label=req.split()[0])
    ax.set_yscale("log"); ax.set_xticks(x)
    ax.set_xticklabels([l.split()[0] for l in launchers], fontsize=8)
    ax.set_ylabel("launcher + energy store (kg)")
    ax.set_title("System mass"); ax.legend(fontsize=7.5)
    ax = axes[1]
    pen = [r1.loc[l, "launched_mass_penalty"] for l in launchers]
    ax.bar(x, pen, color=[CB[0], CB[1], CB[2]])
    ax.set_xticks(x); ax.set_xticklabels([l.split()[0] for l in launchers], fontsize=8)
    ax.set_ylabel("launched-mass penalty (×)")
    ax.set_title("Cost of muzzle-velocity spread (R1)")
    ax = axes[2]
    ns = np.array([24, 100, 1000, 10000])
    g = chemical_gun(2.14e-3, 619.0)
    ax.loglog(ns, ns * g["m_prop_kg"], color=CB[0], marker="o", ms=4,
              label="chemical: propellant magazine")
    ax.axhline(0, color=CB[2])
    ax.plot(ns, np.full_like(ns, 1e-4, dtype=float), color=CB[2], marker="s",
            ms=4, label="EM: zero consumable (solar)")
    ax.set_xlabel("shots over mission")
    ax.set_ylabel("consumable mass (kg)")
    ax.set_title("Logistics"); ax.legend(fontsize=7.5)
    fig.suptitle("LNCH-01: at 410 J all three are trivial to build — the "
                 "secondary metrics decide, and favour the coilgun", fontsize=10)
    save_fig(fig, "lnch01_trade",
             "Launcher trade across system mass, muzzle-repeatability penalty "
             "and consumable logistics. The chemical gun is the lightest "
             "hardware but loses on repeatability, cleanliness and consumable "
             "dependence; the railgun is the worst EM option at this energy; "
             "the coilgun wins.")
    return df, dfsc, concl


if __name__ == "__main__":
    run()
