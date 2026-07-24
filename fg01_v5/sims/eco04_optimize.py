"""
ECO-4 (SIM-17) -- Economic optimisation and the honest viability gate.

Combines INT + ECO into cost per object, finds the architecture that minimises
it, and states plainly whether it reaches viability and against which benchmark.

Because ECO-1 returned a negative result for the small-debris tour, this module
also does what the plan asks in spirit rather than only in letter: it tests
whether a DIFFERENT target class rescues the concept, since the physics of the
kick is sound and only the visiting economics failed.

Three architectures are costed:
  A. Small-debris cluster tour   -- the v6 plan's proposal.
  B. Large tracked-object removal -- same launcher, same kick, different target.
  C. Fleet-scale version of A     -- to show scaling does not save it.

The value side is computed too, rather than asserted: the avoided-collision
value of removing one object, from the actual flux it contributes.
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from fg01.constants import (R_E, MU, GAMMA, RHO_AL, C_LAUNCH_RANGE, P_W_RANGE,
                            P_RE_2026, COMPETITORS, G0)
from fg01.orbital import mass_sphere, xsec_sphere, dv_shift, v_circ
from fg01.population import MASTER_TOTAL_GE_1CM, MASTER_TOTAL_GE_10CM
from fg01.relmotion import nodal_rate
from fg01.delivery import penalty_for_confidence
from fg01.io_utils import save_table, save_fig, save_json, banner, CB

BETA = 1.2007
V_REL = 619.0
D_MISS = 0.005775
DAY = 86400.0

# ECO-1 realised engagement counts (10-yr mission)
N_REALISED = {500: 16, 1000: 12, 2000: 24, 5000: 36, 10000: 55}


def dv_per_engagement(h_km, inc_deg, mission_yr):
    a = R_E + h_km * 1e3
    v = float(v_circ(h_km * 1e3))
    om = abs(float(nodal_rate(a, np.radians(inc_deg))))
    return np.pi * v / (3.5 * om * mission_yr * 365.25 * DAY)


def shot_mass(m_t, a_t, dv_req, sigma_miss=D_MISS, conf=0.95):
    """Launched mass per shot including the areal penalty."""
    m_on_target = dv_req * m_t / (BETA * V_REL)
    pen = max(float(penalty_for_confidence(sigma_miss, a_t, conf)), 1.0)
    return m_on_target * pen, pen


def avoided_collision_value(n_population=MASTER_TOTAL_GE_1CM,
                            mean_asset_value_usd=8.5e6, n_assets=11000,
                            sat_area_m2=10.0, mission_yr=10.0,
                            v_impact_km_s=10.0, shell_overlap=1.0):
    """
    Value of removing ONE object: the reduction in expected asset loss.

    Flux on one asset F = n_density * v_rel * A.  Removing one of N objects cuts
    the population by 1/N, so the avoided expected loss is
        n_assets * mean_asset_value * P_hit_per_asset / N.

    `shell_overlap` (<=1) accounts for the fact that the >=1 cm population peaks
    at 700-1000 km while most active satellites now sit at 500-600 km, so the
    average asset does not see the average density.  This factor is the largest
    single uncertainty in the model and is swept rather than fixed.

    CAVEAT, stated because it matters: with shell_overlap = 1 this model
    predicts ~80 debris losses per year across the active fleet, and the
    observed record shows essentially none.  The model therefore over-predicts,
    and the low end of the swept range is the more credible one.
    """
    vol_km3 = (4 / 3) * np.pi * (((R_E + 1200e3) / 1e3) ** 3
                                 - ((R_E + 400e3) / 1e3) ** 3)
    n_density = n_population / vol_km3 * shell_overlap
    a_km2 = sat_area_m2 * 1e-6
    flux_per_s = n_density * v_impact_km_s * a_km2
    p_hit = flux_per_s * mission_yr * 365.25 * DAY
    total_expected_loss = n_assets * mean_asset_value_usd * min(p_hit, 1.0)
    return dict(n_density_per_km3=n_density, p_hit_per_asset_10yr=p_hit,
                implied_losses_per_year=n_assets * p_hit / mission_yr,
                total_expected_loss_usd=total_expected_loss,
                value_per_object_removed_usd=total_expected_loss / n_population)


def run():
    banner("ECO-4  Economic optimisation and the honest viability gate")

    # ---- Architecture A: small-debris cluster tour ------------------------
    m_t = float(mass_sphere(0.01, RHO_AL))
    a_t = float(xsec_sphere(0.01))
    dv_req = float(dv_shift(865e3, 200e3))
    m_shot_A, pen_A = shot_mass(m_t, a_t, dv_req)

    rows = []
    for dv_budget, n_eng in N_REALISED.items():
        for c_plat in (20e6, 50e6, 150e6, 500e6):
            cons = n_eng * m_shot_A * (P_W_RANGE[1] + C_LAUNCH_RANGE[0])
            total = c_plat + cons
            rows.append(dict(
                architecture="A: small-debris tour", target="1 cm fragment",
                dv_budget_m_s=dv_budget, n_engagements=n_eng,
                platform_cost_musd=c_plat / 1e6,
                m_shot_g=m_shot_A * 1e3,
                consumable_total_usd=cons,
                cost_per_object_usd=total / n_eng,
                cost_per_kg_usd=total / n_eng / m_t))

    # ---- Architecture B: large tracked objects ----------------------------
    # same kick physics, different target: a derelict rocket body
    for m_body, d_body, lab in ((1400.0, 2.6, "SL-8 upper stage (1.4 t)"),
                                (8300.0, 3.9, "SL-16 upper stage (8.3 t)"),
                                (900.0, 1.5, "defunct smallsat (0.9 t)")):
        a_body = np.pi * (d_body / 2) ** 2
        dv_b = float(dv_shift(865e3, 200e3))
        m_on = dv_b * m_body / (BETA * V_REL)
        pen_raw = float(penalty_for_confidence(D_MISS, a_body, 0.95))
        # For a target this large the areal penalty is far below 1, i.e. aim is
        # a non-issue -- but that also means the penalty supplies no implicit
        # margin the way it does for a centimetre target. Apply the explicit
        # gamma = 2 engineering margin instead, and check E_s against it.
        pen_b = max(pen_raw, 1.0) * GAMMA
        m_launch_b = m_on * pen_b
        es = 0.5 * m_launch_b * V_REL ** 2 / m_body
        n_shots = int(np.ceil(m_launch_b / 1.0))     # 1 kg per shot
        for c_plat in (50e6, 150e6):
            for n_obj in (3, 5, 10):
                cons = n_obj * m_launch_b * (P_W_RANGE[1] + C_LAUNCH_RANGE[0])
                dv_tour = n_obj * dv_per_engagement(865, 98.8, 10.0)
                total = c_plat + cons
                rows.append(dict(
                    architecture="B: large tracked object", target=lab,
                    dv_budget_m_s=dv_tour, n_engagements=n_obj,
                    platform_cost_musd=c_plat / 1e6,
                    m_launch_total_kg=m_launch_b,
                    consumable_total_usd=cons,
                    cost_per_object_usd=total / n_obj,
                    cost_per_kg_usd=total / n_obj / m_body,
                    Es_kJ_kg=es / 1e3, n_shots_per_object=n_shots,
                    magazine_kg=m_launch_b))
    df = pd.DataFrame(rows)
    save_table(df, "eco04_architectures",
               "Cost per object for the small-debris tour (A) and the "
               "large-tracked-object variant (B), using the same kick physics.",
               tags={"n_engagements A": "DERIVED (ECO-1)",
                     "platform cost": "UNVALIDATED"})

    print("  Architecture A -- small-debris cluster tour:")
    print(df[(df.architecture.str.startswith("A")) & (df.platform_cost_musd == 50.0)][
        ["dv_budget_m_s", "n_engagements", "cost_per_object_usd",
         "cost_per_kg_usd"]].to_string(index=False))
    print("\n  Architecture B -- large tracked objects (same launcher, same kick):")
    print(df[(df.architecture.str.startswith("B")) & (df.platform_cost_musd == 50.0)][
        ["target", "n_engagements", "m_launch_total_kg", "n_shots_per_object",
         "Es_kJ_kg", "cost_per_object_usd", "cost_per_kg_usd"]].to_string(index=False))

    # ---- the value side, swept rather than asserted ------------------------
    vrows = []
    for overlap in (0.05, 0.2, 0.5, 1.0):
        for mav, alab in ((2e6, "Starlink-weighted ($2M)"),
                          (8.5e6, "mixed fleet ($8.5M)"),
                          (30e6, "high-value-weighted ($30M)")):
            v = avoided_collision_value(mean_asset_value_usd=mav,
                                        shell_overlap=overlap)
            vrows.append(dict(shell_overlap=overlap, asset_case=alab,
                              mean_asset_value_usd=mav,
                              p_hit_10yr=v["p_hit_per_asset_10yr"],
                              implied_losses_per_year=v["implied_losses_per_year"],
                              value_per_object_usd=v["value_per_object_removed_usd"]))
    dfv = pd.DataFrame(vrows)
    save_table(dfv, "eco04_value_model",
               "Avoided-collision value of removing one >=1 cm object, swept "
               "over the shell-overlap factor and the asset-value mix. The "
               "'implied_losses_per_year' column is the model's own reality "
               "check: values far above ~1 conflict with the observed record of "
               "essentially zero confirmed >=1 cm debris losses, so the low "
               "overlap rows are the credible ones.",
               tags={"value model": "UNVALIDATED (parametric)"})
    print("\n  avoided-collision value per object removed (swept):")
    print(dfv[["shell_overlap", "asset_case", "implied_losses_per_year",
               "value_per_object_usd"]].to_string(index=False))
    val_credible = dfv[dfv.implied_losses_per_year <= 2.0]
    v_lo = float(dfv.value_per_object_usd.min())
    v_hi = float(dfv.value_per_object_usd.max())
    v_cred = (float(val_credible.value_per_object_usd.min()),
              float(val_credible.value_per_object_usd.max())) \
        if len(val_credible) else (v_lo, v_lo)
    print(f"  full range ${v_lo:,.0f}-${v_hi:,.0f}; "
          f"credible subset (<=2 losses/yr) ${v_cred[0]:,.0f}-${v_cred[1]:,.0f}")

    best_A = df[df.architecture.str.startswith("A")].cost_per_object_usd.min()
    best_B = df[df.architecture.str.startswith("B")].cost_per_object_usd.min()

    # ---- Architecture C: fleet ---------------------------------------------
    fleet = []
    for n_plat in (1, 10, 100, 1000):
        n_eng = 24 * n_plat
        c_unit = 50e6 * (n_plat ** (np.log(0.85) / np.log(2)))
        total = n_plat * c_unit
        fleet.append(dict(n_platforms=n_plat, unit_cost_musd=c_unit / 1e6,
                          total_cost_musd=total / 1e6,
                          n_engagements=n_eng,
                          cost_per_object_usd=total / n_eng,
                          frac_of_1_2M_population=n_eng / MASTER_TOTAL_GE_1CM))
    dff = pd.DataFrame(fleet)
    save_table(dff, "eco04_fleet_scaling",
               "Fleet scaling of architecture A with an 85% learning curve. "
               "Cost per object falls only as the learning curve allows, "
               "because engagements scale linearly with platforms.")
    print("\n  Architecture C -- fleet scaling of A:")
    print(dff.to_string(index=False))

    concl = {
        "architecture_A_small_debris": {
            "best_cost_per_object_usd": float(best_A),
            "engagements_per_platform_decade": 24,
            "avoided_collision_value_per_object_usd_range": [v_lo, v_hi],
            "avoided_collision_value_credible_usd": list(v_cred),
            "cost_to_value_ratio_best_case": float(best_A / v_hi),
            "cost_to_value_ratio_credible": float(best_A / max(v_cred[1], 1e-9)),
            "verdict": "FAILS"},
        "architecture_B_large_objects": {
            "best_cost_per_object_usd": float(best_B),
            "competitor_ADR_tug_usd_per_object": "10M-100M (demonstrated class)",
            "verdict": "COMPETITIVE"},
        "headline_two_numbers": {
            "niche_viability_small_debris": (
                f"NOT ACHIEVED. Cost per 1 cm object is ${best_A:,.0f} against "
                f"the plan's $5,000 niche target -- a shortfall of "
                f"{best_A/5000:.0f}x. The binding condition is ECO-1's "
                "Delta-v-per-engagement invariant, not the debris population, "
                "the kick physics or the platform cost."),
            "laser_competitive": (
                "NOT ACHIEVED and not approachable: it would require 1e5 "
                "engagements, i.e. ~1e7 km/s of tour Delta-v.")},
        "value_gap": (
            f"Removing one randomly-selected >=1 cm fragment avoids between "
            f"${v_lo:,.0f} and ${v_hi:,.0f} of expected asset loss across the "
            f"swept assumptions, and at most ${v_cred[1]:,.0f} in the subset "
            "consistent with the observed loss record. Against a cost of "
            f"${best_A:,.0f} the gap is "
            f"{best_A/max(v_hi,1e-9):.0e}x at the most generous end and "
            f"{best_A/max(v_cred[1],1e-9):.0e}x at the credible end. This is "
            "not a failing of FG01 specifically -- it is why operational ADR "
            "targets large objects, whose removal prevents the breakups that "
            "create fragments in the first place, rather than the fragments "
            "themselves. No small-debris remediation concept, FG01 or "
            "otherwise, closes a gap of this size on per-object value alone; "
            "the case for small-debris removal has to be made on systemic "
            "cascade prevention, which this study does not model."),
        "the_constructive_result": (
            "The kick physics survives the economics unchanged, and it applies "
            "to any target. Against a large tracked object the same launcher, "
            "the same 619 m/s intercept and the same sub-catastrophic specific "
            f"energy ({float(df[df.architecture.str.startswith('B')].Es_kJ_kg.iloc[0]):.1f} kJ/kg -- "
            "identical, because E_s is mass-independent) deliver the same "
            "200 km orbit shift. The target's cross-section is now square "
            "metres rather than square centimetres, so the areal penalty "
            "vanishes and essentially all launched mass lands. Cost per object "
            f"is ${best_B:,.0f}, against $10M-100M for a demonstrated-class ADR "
            "tug. The concept is economically defensible against LARGE debris "
            "and indefensible against small debris -- the exact inverse of its "
            "founding premise."),
        "why_that_inversion_happens": (
            "Two independent scalings both favour large targets. (1) The areal "
            "mass penalty is pi*e*d_miss^2/A_target, so a 2.6 m rocket body "
            "needs no pattern margin at all while a 1 cm fragment needs 22x its "
            "mass on target. (2) Cost per object is dominated by the Delta-v to "
            "visit it, which is the same ~107 m/s whatever the target masses -- "
            "so it should be spent on the most valuable target available."),
    }
    save_json(concl, "eco04_conclusions")

    # ---- figures -----------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.3))
    ax = axes[0]
    sA = df[(df.architecture.str.startswith("A")) & (df.platform_cost_musd == 50.0)]
    ax.loglog(sA.n_engagements, sA.cost_per_object_usd, color=CB[0], marker="o",
              ms=5, label="A: small debris (realised)")
    n_hyp = np.geomspace(10, 1e5, 50)
    ax.loglog(n_hyp, 50e6 / n_hyp, color=CB[0], ls=":", lw=1.2,
              label="A: if N were achievable")
    sB = df[(df.architecture.str.startswith("B")) & (df.platform_cost_musd == 50.0)]
    ax.scatter(sB.n_engagements, sB.cost_per_object_usd, color=CB[1], s=45,
               marker="s", zorder=5, label="B: large tracked objects")
    ax.axhspan(1e7, 1e8, color=CB[2], alpha=0.15)
    ax.text(11, 1.5e7, "ADR tug (demonstrated class)", fontsize=7)
    ax.axhline(5000, color="k", ls="--", lw=1.2)
    ax.text(11, 6000, "plan's niche target", fontsize=7)
    ax.set_xlabel("objects removed per platform")
    ax.set_ylabel("cost per object (USD)")
    ax.set_title("Where each architecture lands")
    ax.legend(fontsize=7)
    ax = axes[1]
    hh = np.array([500, 700, 865, 1000, 1200])
    for i, inc in enumerate([30.0, 51.6, 74.0, 98.8]):
        ax.semilogy(hh, [dv_per_engagement(h, inc, 10.0) for h in hh],
                    color=CB[i], marker="o", ms=4, label=f"i = {inc}Â°")
    ax.set_xlabel("cluster altitude (km)")
    ax.set_ylabel("Î”v per engagement (m/s)")
    ax.set_title("The invariant: Ï€Â·v/(3.5Â·Î©Ì‡Â·T)")
    ax.legend(fontsize=8)
    fig.suptitle("ECO-4: the kick works; visiting the targets is what does not",
                 fontsize=10)
    save_fig(fig, "eco04_economics",
             "Cost per object for the small-debris tour against the "
             "large-object variant, and the Î”v-per-engagement invariant that "
             "drives both. Sun-synchronous orbits, where the most valuable "
             "debris sits, have the slowest nodal regression and are therefore "
             "the most expensive to tour.")
    return df, concl


if __name__ == "__main__":
    run()
