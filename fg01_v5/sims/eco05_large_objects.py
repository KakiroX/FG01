"""
ECO-5 -- Is the kick actually EFFICIENT against large objects, or only cheap
against an expensive benchmark?

ECO-4 reported architecture B (large tracked objects) as "competitive" on cost
per object against a demonstrated-class ADR tug. That comparison was against tug
PROGRAMME costs, which are dominated by development, rendezvous and capture
hardware. It did not ask the prior question: how efficiently does a momentum
kick deliver impulse at all?

The answer is a single number. A rocket ejects propellant at its exhaust
velocity; FG01 ejects tungsten at the intercept velocity. So the kick has an
effective exhaust velocity

    v_eff = beta * v_rel = 1.2 * 619 = 743 m/s   ->   I_sp_equivalent = 75.7 s

which is worse than cold gas. The projectile mass needed is therefore

    m_p / m_target = gamma * dv / (beta * v_rel)

a FIXED FRACTION of the target's mass, independent of target size. At the design
point that is 14%. For a 1 cm fragment 14% of 1.4 g is 0.2 g and nobody cares.
For an 8.3 t upper stage it is 1.17 t of tungsten that must be launched.

This module quantifies where that crossover bites.
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from fg01.constants import (GAMMA, G0, C_LAUNCH_RANGE, P_W_RANGE, R_E)
from fg01.orbital import dv_shift, v_circ
from fg01.relmotion import nodal_rate
from fg01.io_utils import save_table, save_fig, save_json, banner, CB

BETA = 1.2007
V_REL = 619.0
DV_REQ = float(dv_shift(865e3, 200e3))
DAY = 86400.0
C_LAUNCH = C_LAUNCH_RANGE[0]        # $3000/kg, favourable end
P_W = P_W_RANGE[1]                  # $100/kg tungsten, conservative end


def dv_per_engagement(h_km=865.0, inc_deg=98.8, mission_yr=10.0):
    a = R_E + h_km * 1e3
    om = abs(float(nodal_rate(a, np.radians(inc_deg))))
    return np.pi * float(v_circ(h_km * 1e3)) / (3.5 * om * mission_yr * 365.25 * DAY)


def run():
    banner("ECO-5  Is the kick efficient against large objects?")

    v_eff = BETA * V_REL
    isp_eq = v_eff / G0
    frac = GAMMA * DV_REQ / (BETA * V_REL)
    print(f"  effective exhaust velocity  v_eff = beta*v_rel = {v_eff:.0f} m/s")
    print(f"  equivalent specific impulse       = {isp_eq:.1f} s")
    print(f"  projectile mass as fraction of target mass = {frac:.3f} "
          f"({frac*100:.1f}%)")
    print(f"  (cold gas ~70 s, hydrazine ~230 s, bipropellant ~320 s, "
          f"Hall ~1600 s, gridded ion ~3000 s)")

    # ---- head-to-head impulse delivery -------------------------------------
    rows = []
    for m_t, lab in ((1.41e-3, "1 cm fragment (1.4 g)"),
                     (0.176, "5 cm fragment (176 g)"),
                     (10.0, "10 kg object"),
                     (100.0, "100 kg object"),
                     (900.0, "defunct smallsat (0.9 t)"),
                     (1400.0, "SL-8 upper stage (1.4 t)"),
                     (8300.0, "SL-16 upper stage (8.3 t)")):
        m_fg01 = GAMMA * DV_REQ * m_t / (BETA * V_REL)
        for isp, plab, m_veh in ((230.0, "hydrazine tug", 300.0),
                                 (320.0, "bipropellant tug", 400.0),
                                 (1600.0, "Hall-thruster tug", 500.0),
                                 (3000.0, "gridded-ion tug", 500.0)):
            m_prop = (m_t + m_veh) * (np.exp(DV_REQ / (isp * G0)) - 1)
            rows.append(dict(
                target=lab, m_target_kg=m_t, competitor=plab, isp_s=isp,
                fg01_projectile_kg=m_fg01,
                tug_propellant_kg=m_prop,
                fg01_over_tug=m_fg01 / m_prop,
                fg01_launch_cost_usd=m_fg01 * (P_W + C_LAUNCH),
                tug_propellant_launch_cost_usd=m_prop * C_LAUNCH,
                consumable_delta_usd=m_fg01 * (P_W + C_LAUNCH) - m_prop * C_LAUNCH))
    df = pd.DataFrame(rows)
    save_table(df, "eco05_impulse_efficiency",
               "Mass that must be launched to deliver the same 51.8 m/s kick: "
               "FG01 projectile versus tug propellant. The tug figure includes "
               "accelerating the tug's own mass.",
               tags={"v_eff": "DERIVED", "Isp values": "SOURCED (standard)"})
    print("\n  mass to deliver the same 51.8 m/s kick (vs gridded-ion tug):")
    print(df[df.competitor == "gridded-ion tug"][
        ["target", "fg01_projectile_kg", "tug_propellant_kg", "fg01_over_tug",
         "consumable_delta_usd"]].to_string(index=False))

    # ---- full mission accounting, platform scaled to the magazine ---------
    # THIS is what ECO-4 did not do: the platform must carry the magazine.
    mission = []
    dv_tour_per = dv_per_engagement()
    for m_t, lab in ((900.0, "smallsat 0.9 t"), (1400.0, "SL-8 1.4 t"),
                     (8300.0, "SL-16 8.3 t")):
        m_shot = GAMMA * DV_REQ * m_t / (BETA * V_REL)
        for n_obj in (3, 5, 10):
            magazine = n_obj * m_shot
            # platform dry mass: fixed bus + structure to carry the magazine
            m_bus = 400.0
            m_struct = 0.25 * magazine            # tankage/structure fraction
            m_dry = m_bus + m_struct
            m_wet0 = m_dry + magazine
            # tour delta-v must move the whole stack; recoil must be cancelled
            dv_tour = n_obj * dv_tour_per
            recoil_impulse = magazine * V_REL
            dv_recoil = recoil_impulse / max((m_dry + magazine / 2), 1.0)
            dv_total = dv_tour + dv_recoil
            isp = 3000.0
            m_prop = m_wet0 * (np.exp(dv_total / (isp * G0)) - 1)
            m_launch_total = m_wet0 + m_prop
            # costs
            c_plat = 50e6
            c_launch_stack = m_launch_total * C_LAUNCH
            c_material = magazine * P_W
            total = c_plat + c_launch_stack + c_material
            mission.append(dict(
                target=lab, m_target_kg=m_t, n_objects=n_obj,
                m_shot_kg=m_shot, magazine_kg=magazine,
                m_dry_kg=m_dry, propellant_kg=m_prop,
                dv_tour_m_s=dv_tour, dv_recoil_m_s=dv_recoil,
                launched_stack_kg=m_launch_total,
                launch_cost_musd=c_launch_stack / 1e6,
                total_cost_musd=total / 1e6,
                cost_per_object_musd=total / n_obj / 1e6,
                cost_per_kg_removed=total / (n_obj * m_t)))
    dfm = pd.DataFrame(mission)
    save_table(dfm, "eco05_mission_scaled",
               "Full mission accounting with the platform scaled to carry its "
               "magazine -- the step ECO-4 omitted. Recoil delta-v is included: "
               "firing the magazine retrograde pushes the platform prograde.",
               tags={"structure fraction": "UNVALIDATED (0.25)"})
    print("\n  mission accounting with the platform scaled to the magazine:")
    print(dfm[["target", "n_objects", "magazine_kg", "propellant_kg",
               "dv_recoil_m_s", "launched_stack_kg", "cost_per_object_musd",
               "cost_per_kg_removed"]].to_string(index=False))

    # ---- crossover ---------------------------------------------------------
    cross = []
    for m_t in np.geomspace(1.0, 2e4, 60):
        m_shot = GAMMA * DV_REQ * m_t / (BETA * V_REL)
        n_obj = 10
        magazine = n_obj * m_shot
        m_dry = 400.0 + 0.25 * magazine
        m_wet0 = m_dry + magazine
        dv_total = n_obj * dv_tour_per + magazine * V_REL / max(m_dry + magazine / 2, 1)
        m_prop = m_wet0 * (np.exp(dv_total / (3000 * G0)) - 1)
        m_stack = m_wet0 + m_prop
        fg01_cost = 50e6 + m_stack * C_LAUNCH + magazine * P_W
        # electric tug servicing the same 10 objects
        m_prop_tug = sum((m_t + 500.0) * (np.exp(DV_REQ / (3000 * G0)) - 1)
                         for _ in range(n_obj))
        m_tug_stack = 500.0 + m_prop_tug + 500.0 * (
            np.exp(n_obj * dv_tour_per / (3000 * G0)) - 1)
        tug_cost = 50e6 + m_tug_stack * C_LAUNCH
        cross.append(dict(
            m_target_kg=m_t,
            fg01_cost_per_object_musd=fg01_cost / n_obj / 1e6,
            tug_cost_per_object_musd=tug_cost / n_obj / 1e6,
            fg01_penalty_musd=(fg01_cost - tug_cost) / n_obj / 1e6,
            fg01_launched_stack_kg=m_stack,
            tug_launched_stack_kg=m_tug_stack))
    dfx = pd.DataFrame(cross)
    save_table(dfx, "eco05_crossover",
               "Cost per object for FG01 against an electric tug of equal "
               "platform cost servicing the same 10 objects, as target mass "
               "varies. Both carry the same tour Delta-v.")
    # where does the extra cost exceed a plausible value of not docking?
    DOCKING_VALUE_MUSD = 2.0
    over = dfx[dfx.fg01_penalty_musd > DOCKING_VALUE_MUSD]
    m_cross = float(over.m_target_kg.min()) if len(over) else np.inf

    def at_mass(m):
        i = int(np.argmin(np.abs(dfx.m_target_kg.to_numpy() - m)))
        return dfx.iloc[i]

    pen_1400 = float(at_mass(1400).fg01_penalty_musd)
    pen_8300 = float(at_mass(8300).fg01_penalty_musd)
    # where is FG01 actually cheaper than the tug?
    cheaper = dfx[dfx.fg01_penalty_musd < 0]
    m_cheap_max = float(cheaper.m_target_kg.max()) if len(cheaper) else 0.0
    print(f"\n  FG01 vs an equal-platform electric tug, per object:")
    print(f"    at 100 kg : {float(at_mass(100).fg01_penalty_musd):+.2f} M USD")
    print(f"    at 1.4 t  : {pen_1400:+.2f} M USD")
    print(f"    at 8.3 t  : {pen_8300:+.2f} M USD")
    print(f"  FG01 is CHEAPER than the tug below ~{m_cheap_max:,.0f} kg "
          f"(the tug has to accelerate its own mass)")
    print(f"  crossover (penalty exceeds a ${DOCKING_VALUE_MUSD:.0f}M value for "
          f"not having to dock): target mass ~{m_cross:,.0f} kg")

    concl = {
        "effective_exhaust_velocity_m_s": float(v_eff),
        "equivalent_isp_s": float(isp_eq),
        "projectile_mass_fraction_of_target": float(frac),
        "headline": (
            f"A momentum kick has an effective exhaust velocity of beta*v_rel = "
            f"{v_eff:.0f} m/s, i.e. an equivalent specific impulse of "
            f"{isp_eq:.0f} seconds -- below cold gas. It is therefore an "
            "intrinsically INEFFICIENT way to deliver impulse. The projectile "
            f"mass is a fixed {frac*100:.0f}% of the target's mass regardless of "
            "target size, so the penalty is invisible for a 1.4 g fragment "
            "(0.2 g of tungsten) and dominant for an 8.3 t stage (1.17 t)."),
        "correction_to_ECO4": (
            "ECO-4 reported architecture B as competitive at $5.6M-20M per "
            "object. That comparison was against ADR tug PROGRAMME costs and "
            "did not scale the platform to carry its own magazine. Doing both "
            "-- charging structure for the magazine, propellant to move the "
            "loaded stack, and the recoil delta-v from firing it -- raises the "
            "cost and, more importantly, shows the comparison was against the "
            "wrong quantity. Against an electric tug of equal platform cost "
            "servicing the same objects, FG01 is more expensive at every target "
            "mass, and the gap widens linearly with it."),
        "the_real_trade": (
            "FG01 trades propulsive efficiency for not having to dock. A tug "
            "must rendezvous, match rates with a tumbling derelict, grapple it "
            "and control the mated stack -- the hard, expensive, still largely "
            "undemonstrated part of ADR. FG01 stands off at 10 m and shoots. "
            "That is a real advantage, and the honest question is whether it is "
            f"worth the consumable penalty, which is ${pen_1400:.1f}M per object "
            f"at 1.4 t and ${pen_8300:.1f}M at 8.3 t, growing linearly with "
            "target mass."),
        "the_symmetry": (
            "There is a genuine interior optimum in target mass, from two "
            "opposing effects. BELOW about "
            f"{m_cheap_max:,.0f} kg FG01 needs LESS launched mass than a tug, "
            "because a tug must accelerate its own ~500 kg to deliver the kick "
            "while the projectile carries no dead mass -- delivering 51.8 m/s "
            "to a 1.4 g fragment costs 0.2 mg of tungsten against 0.89 kg of "
            "xenon. ABOVE it the 14% projectile fraction dominates and the tug "
            "wins by 50-75x. The kick is mass-efficient exactly where the "
            "target is small -- and that is precisely where the tour Delta-v "
            "makes it uneconomic anyway."),
        "recoil_finding": (
            "Firing the magazine retrograde pushes the platform prograde. For "
            "10 SL-16 engagements the magazine is 11.7 t and the recoil "
            "Delta-v is 789 m/s -- larger than the entire tour budget, and a "
            "term neither the v6 plan nor ECO-4 accounted for. It costs 976 kg "
            "of xenon to cancel. Recoil scales with magazine mass, so it is "
            "another penalty that grows with target mass."),
        "crossover_kg": float(m_cross),
        "fg01_cheaper_below_kg": float(m_cheap_max),
        "verdict": (
            "EFFICIENT: no. Competitive within a bounded range: plausibly, for "
            f"targets below roughly {m_cross:,.0f} kg, IF avoiding docking is "
            "worth a few million dollars per object -- which is a programmatic "
            "judgement this study cannot settle. Above that mass the 14% "
            "projectile fraction dominates and the case fails. The sweet spot "
            "is NOT the multi-tonne upper stages ECO-4 highlighted; it is "
            "intermediate objects of tens to hundreds of kilograms, which are "
            "too big for a laser, too small to justify a bespoke tug, and "
            "tumbling enough to make docking unattractive."),
        "why_small_debris_still_fails": (
            "Note this does NOT rescue small debris. There the projectile "
            "fraction is irrelevant (0.2 g) and the killer is the 107 m/s of "
            "tour Delta-v per engagement, which is charged per OBJECT and is "
            "identical whatever the object masses. That is why cost per object "
            "is roughly flat while cost per kilogram falls by nine orders of "
            "magnitude across the size range."),
    }
    save_json(concl, "eco05_conclusions")

    # ---- figures -----------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.3))
    ax = axes[0]
    ax.loglog(dfx.m_target_kg, dfx.fg01_cost_per_object_musd, color=CB[0], lw=2,
              label="FG01 momentum kick")
    ax.loglog(dfx.m_target_kg, dfx.tug_cost_per_object_musd, color=CB[1], lw=2,
              label="electric tug, equal platform cost")
    ax.axvline(m_cross, color="0.4", ls=":", lw=1.4)
    ax.text(m_cross * 1.1, 6, f"~{m_cross:,.0f} kg", fontsize=7.5, rotation=90)
    ax.set_xlabel("target mass (kg)")
    ax.set_ylabel("cost per object (M USD)")
    ax.set_title("Both service 10 objects with the same tour Δv")
    ax.legend(fontsize=8)
    ax = axes[1]
    isps = np.array([75.7, 70, 230, 320, 1600, 3000])
    labs = ["FG01 kick\n(β·v_rel)", "cold gas", "hydrazine", "bipropellant",
            "Hall", "gridded ion"]
    cols = [CB[1]] + [CB[0]] * 5
    ax.barh(np.arange(len(isps)), isps, color=cols)
    ax.set_yticks(np.arange(len(isps)))
    ax.set_yticklabels(labs, fontsize=8)
    ax.set_xscale("log")
    ax.set_xlabel("specific impulse (s)")
    ax.set_title("The kick as a propulsion system")
    fig.suptitle("ECO-5: a momentum kick is an inefficient way to move mass — "
                 "the question is whether not docking is worth it", fontsize=10)
    save_fig(fig, "eco05_efficiency",
             "Left: cost per object against target mass for FG01 versus an "
             "electric tug of equal platform cost. Right: the kick's effective "
             "specific impulse against real propulsion systems. FG01 must "
             "launch 14% of the target's mass as tungsten, which is negligible "
             "for grams and decisive for tonnes.")
    return df, dfm, dfx, concl


if __name__ == "__main__":
    run()
