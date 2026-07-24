"""
SIM-9 -- Cost-economic model and competitor comparison.

AM-8: misses are in the cost.  Under the delivery model (SIM-3/4) a shot that
misses by more than planned delivers less momentum rather than nothing, so the
"cost of misses" enters as the extra shots needed to complete a kick at the
chosen delivery confidence, plus the mass wasted on every shot by the areal
penalty -- which is by far the larger term and is charged in full.

Two cost metrics are reported because they differ by three orders of magnitude
for centimetre debris and each is defensible:
  USD per kg removed   -- comparable with the published competitor figures.
  USD per object removed -- arguably the better metric for the 1-10 cm class,
                            where collision risk scales with object count.
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from fg01.constants import (P_RE_2026, P_RE_2025, P_W_RANGE, C_LAUNCH_RANGE,
                            C_LAUNCH_STARSHIP, P_ELEC, C_PLATFORM_RANGE,
                            COMPETITORS, RHO_AL, GAMMA, ALT_SWEEP_KM)
from fg01.orbital import mass_sphere, xsec_sphere, dv_shift
from fg01.delivery import penalty_for_confidence
from fg01.io_utils import save_table, save_fig, save_json, banner, CB

CONF = 0.9


def shots_per_object(h_km, d_debris_m, staged):
    """Kicks needed (SIM-2) divided by the per-shot delivery confidence."""
    n_kick = staged.get((h_km, d_debris_m * 100), 1)
    n_kick = max(int(n_kick) if np.isfinite(n_kick) else 6, 1)
    return n_kick / CONF, n_kick


def run():
    banner("SIM-9  Cost-economic model and competitor comparison")

    # kicks per object from SIM-2
    st = pd.read_csv(pathlib.Path(__file__).resolve().parent.parent
                     / "outputs" / "sim02_staged_kicks.csv")
    st = st[(st.rule == "25yr") & (st.solar_case == "nominal")]
    staged = {(int(r.h_km), r.d_debris_cm): r.n_kicks for _, r in st.iterrows()}

    rows = []
    for d in (0.01, 0.02, 0.05, 0.10):
        m_t = float(mass_sphere(d, RHO_AL))
        a_t = float(xsec_sphere(d))
        for h in ALT_SWEEP_KM:
            dv = float(dv_shift(h * 1e3, 200e3))
            for beta in (1.0, 1.5, 2.5):
                for v_rel in (619.0,):
                    m_on_target = GAMMA * dv * m_t / (beta * v_rel)
                    for sigma_miss, lab in ((0.0058, "design (R=10 m, 200 µrad)"),
                                            (0.02, "20 mm budget"),
                                            (0.10, "100 mm budget")):
                        pen = float(penalty_for_confidence(sigma_miss, a_t, CONF))
                        m_shot = m_on_target * pen
                        n_shot, n_kick = shots_per_object(h, d, staged)
                        for p_mat, mat in ((P_RE_2026, "Re (2026-07-22 spot)"),
                                           (P_RE_2025, "Re (2025 avg)"),
                                           (P_W_RANGE[1], "W (upper)")):
                            for c_l, l_lab in ((C_LAUNCH_RANGE[0], "3000 dedicated"),
                                               (C_LAUNCH_RANGE[1], "7000 rideshare")):
                                for n_eng in (1e3, 1e4, 1e5):
                                    for c_pf in C_PLATFORM_RANGE:
                                        c_mat = m_shot * p_mat
                                        c_lau = m_shot * c_l
                                        c_ele = 0.5 * m_shot * v_rel ** 2 / 0.4 / 3.6e6 * P_ELEC
                                        c_var = (c_mat + c_lau + c_ele) * n_shot
                                        c_pfo = c_pf / n_eng
                                        c_obj = c_var + c_pfo
                                        rows.append(dict(
                                            d_debris_cm=d * 100, h_km=h, beta=beta,
                                            miss_budget=lab, sigma_miss_m=sigma_miss,
                                            material=mat, launch=l_lab,
                                            n_engagements_life=n_eng,
                                            C_platform=c_pf,
                                            m_on_target_g=m_on_target * 1e3,
                                            mass_penalty=pen,
                                            m_shot_g=m_shot * 1e3,
                                            n_kicks=n_kick, n_shots=n_shot,
                                            cost_material=c_mat * n_shot,
                                            cost_launch=c_lau * n_shot,
                                            cost_energy=c_ele * n_shot,
                                            cost_platform=c_pfo,
                                            cost_per_object=c_obj,
                                            cost_per_kg=c_obj / m_t))
    df = pd.DataFrame(rows)
    save_table(df, "sim09_cost_full",
               "Full cost model: per-object and per-kg cost across debris size, "
               "altitude, beta, guidance budget, material, launch price, "
               "platform cost and engagement count.",
               tags={"prices": "SOURCED dated", "platform": "UNVALIDATED",
                     "n_engagements": "DESIGN"})

    base = df[(df.beta == 1.5) & (df.miss_budget == "design (R=10 m, 200 µrad)")
              & (df.material == "Re (2026-07-22 spot)") & (df.launch == "3000 dedicated")
              & (df.n_engagements_life == 1e4) & (df.C_platform == 50e6)]
    print("\n  design case, β=1.5, Re spot, $3000/kg launch, 1e4 engagements, "
          "$50M platform:")
    print(base[base.h_km.isin([600, 800, 1000])][
        ["d_debris_cm", "h_km", "m_shot_g", "n_kicks", "cost_material",
         "cost_launch", "cost_platform", "cost_per_object", "cost_per_kg"]]
        .to_string(index=False))

    # breakdown at the reference point
    ref = base[(base.d_debris_cm == 1.0) & (base.h_km == 800)].iloc[0]
    save_table(pd.DataFrame([{
        "component": k,
        "USD": float(ref[k]),
        "share": float(ref[k] / ref.cost_per_object)}
        for k in ("cost_material", "cost_launch", "cost_energy", "cost_platform")]),
        "sim09_cost_breakdown",
        "Cost breakdown per 1 cm object removed at 800 km, design guidance "
        "budget, rhenium at 2026 spot, $3000/kg launch, $50M platform "
        "amortised over 10,000 engagements.")

    # competitor comparison
    m_ref = float(mass_sphere(0.01, RHO_AL))
    comp = []
    for name, (lo, hi, status) in COMPETITORS.items():
        comp.append(dict(method=name, status=status,
                         usd_per_kg_low=lo, usd_per_kg_high=hi,
                         usd_per_1cm_object_low=lo * m_ref,
                         usd_per_1cm_object_high=hi * m_ref,
                         targets_1_to_10cm_class=name.startswith("Ground")))
    fg = df[(df.d_debris_cm == 1.0) & (df.h_km == 800) & (df.beta == 1.5)
            & (df.n_engagements_life == 1e4)]
    comp.append(dict(method="FG01 (this study, Re)", status="concept",
                     usd_per_kg_low=float(fg[fg.material.str.startswith("Re (2025")].cost_per_kg.min()),
                     usd_per_kg_high=float(fg.cost_per_kg.max()),
                     usd_per_1cm_object_low=float(fg.cost_per_object.min()),
                     usd_per_1cm_object_high=float(fg.cost_per_object.max()),
                     targets_1_to_10cm_class=True))
    fgw = fg[fg.material.str.startswith("W")]
    comp.append(dict(method="FG01 (this study, W substituted)", status="concept",
                     usd_per_kg_low=float(fgw.cost_per_kg.min()),
                     usd_per_kg_high=float(fgw.cost_per_kg.max()),
                     usd_per_1cm_object_low=float(fgw.cost_per_object.min()),
                     usd_per_1cm_object_high=float(fgw.cost_per_object.max()),
                     targets_1_to_10cm_class=True))
    dfc = pd.DataFrame(comp)
    save_table(dfc, "sim09_competitors",
               "Side-by-side with published competitor costs. Per-object "
               "columns convert the per-kg figures using a 1 cm Al sphere "
               "(1.414 g) so the two conventions can be compared directly.")
    print("\n", dfc.to_string(index=False))

    # ---- the dominant lever: engagements per platform lifetime -------------
    eng = df[(df.d_debris_cm == 1.0) & (df.h_km == 800) & (df.beta == 1.5)
             & (df.miss_budget == "design (R=10 m, 200 µrad)")
             & (df.material == "Re (2026-07-22 spot)")
             & (df.launch == "3000 dedicated")][
        ["n_engagements_life", "C_platform", "cost_material", "cost_launch",
         "cost_platform", "cost_per_object", "cost_per_kg"]]
    save_table(eng, "sim09_engagement_sensitivity",
               "Cost per object versus the number of engagements a platform "
               "completes in its life. Because consumables are ~$31 per object "
               "and platform amortisation is everything else, the engagement "
               "count is the single dominant economic variable -- and it is the "
               "one this study does not model (SIM-11 caveat, R10).")
    print("\n  cost per object vs engagements per platform life:")
    print(eng.to_string(index=False))

    # cost vs target size -- the key structural result
    sizes = df[(df.beta == 1.5) & (df.h_km == 800)
               & (df.material == "Re (2026-07-22 spot)")
               & (df.launch == "3000 dedicated") & (df.n_engagements_life == 1e4)
               & (df.C_platform == 50e6)]
    concl = {
        "cost_per_object_range_1cm_800km": [float(fg.cost_per_object.min()),
                                            float(fg.cost_per_object.max())],
        "cost_per_kg_range_1cm_800km": [float(fg.cost_per_kg.min()),
                                        float(fg.cost_per_kg.max())],
        "dominant_cost_component": max(
            ("cost_material", "cost_launch", "cost_platform"),
            key=lambda k: float(ref[k])),
        "breakdown_at_reference": {k: float(ref[k]) for k in
                                   ("cost_material", "cost_launch", "cost_energy",
                                    "cost_platform")},
        "size_scaling": {
            f"{r.d_debris_cm:g} cm": {"usd_per_object": float(r.cost_per_object),
                                      "usd_per_kg": float(r.cost_per_kg),
                                      "mass_penalty": float(r.mass_penalty)}
            for _, r in sizes[sizes.miss_budget == "design (R=10 m, 200 µrad)"].iterrows()},
        "guidance_budget_sensitivity": {
            r.miss_budget: float(r.cost_per_object)
            for _, r in df[(df.d_debris_cm == 1.0) & (df.h_km == 800)
                           & (df.beta == 1.5)
                           & (df.material == "Re (2026-07-22 spot)")
                           & (df.launch == "3000 dedicated")
                           & (df.n_engagements_life == 1e4)
                           & (df.C_platform == 50e6)].iterrows()},
        "material_switch_saving_pct": float(
            100 * (1 - fgw.cost_per_object.median() / fg[fg.material.str.startswith("Re (2026")]
                   .cost_per_object.median())),
        "starship_sensitivity_note": (
            f"Launch prices of {C_LAUNCH_STARSHIP[0]}-{C_LAUNCH_STARSHIP[1]} USD/kg "
            "are reported only as a labelled sensitivity (AM-9) and are not "
            "used in any headline figure."),
        "AM8_statement": (
            "Net removal is positive at the design point (SIM-5), so cost per "
            "kg removed is defined. Had the design point been catastrophic, "
            "net removal would be negative and this figure would be reported "
            "as undefined rather than as a large finite number."),
        "consumables_per_object_usd": float(ref.cost_material + ref.cost_launch
                                            + ref.cost_energy),
        "platform_share_of_cost": float(ref.cost_platform / ref.cost_per_object),
        "headline_finding": (
            "The pivot makes the *consumables* almost free: at the design "
            "guidance budget a 1 cm object at 800 km costs about $31 in "
            "rhenium, launch of that rhenium, and electricity combined. "
            "Essentially the entire cost is amortisation of the platform, so "
            "cost per object is C_platform / N_engagements and nothing else "
            "matters much. At $50M and 10,000 engagements that is $5,000 per "
            "object; at 100,000 engagements it is $500. The economic question "
            "for FG01 is therefore not the physics of the kick at all -- it is "
            "how many objects one platform can find, approach and service in "
            "its lifetime, which this study does not model."),
        "competitor_comparison_caveat": (
            "Published competitor costs are quoted per kg and were derived for "
            "targets three to five orders of magnitude more massive. Converting "
            "them to a 1.4 g object implies an electrodynamic tether could "
            "remove a 1 cm fragment for $1.41, which is meaningless: tethers, "
            "aerogel capture and robotic tugs cannot engage the 1-10 cm class "
            "at all. The only method that genuinely competes for this size "
            "class is ground-based laser ablation, and against it FG01 is "
            "3-4 orders of magnitude more expensive per object, because the "
            "laser's hardware stays on the ground. That gap is the concept's "
            "central economic problem and no part of the v5 pivot addresses it."),
    }
    save_json(concl, "sim09_conclusions")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.1))
    ax = axes[0]
    for i, mb in enumerate(("design (R=10 m, 200 µrad)", "20 mm budget", "100 mm budget")):
        s = df[(df.beta == 1.5) & (df.h_km == 800) & (df.miss_budget == mb)
               & (df.material == "Re (2026-07-22 spot)")
               & (df.launch == "3000 dedicated") & (df.n_engagements_life == 1e4)
               & (df.C_platform == 50e6)].sort_values("d_debris_cm")
        ax.loglog(s.d_debris_cm, s.cost_per_kg, color=CB[i], marker="o", ms=4,
                  label=mb)
    for j, (name, (lo, hi, _)) in enumerate(COMPETITORS.items()):
        ax.axhspan(lo, hi, color=CB[(j + 3) % 8], alpha=0.12)
        ax.text(1.02, hi * 0.6, name, fontsize=6.5, color="0.25")
    ax.set_xlabel("debris diameter (cm)"); ax.set_ylabel("USD per kg removed")
    ax.set_title("Cost per kg vs target size"); ax.legend(fontsize=7)
    ax = axes[1]
    for i, mb in enumerate(("design (R=10 m, 200 µrad)", "20 mm budget", "100 mm budget")):
        s = df[(df.beta == 1.5) & (df.h_km == 800) & (df.miss_budget == mb)
               & (df.material == "Re (2026-07-22 spot)")
               & (df.launch == "3000 dedicated") & (df.n_engagements_life == 1e4)
               & (df.C_platform == 50e6)].sort_values("d_debris_cm")
        ax.loglog(s.d_debris_cm, s.cost_per_object, color=CB[i], marker="o", ms=4,
                  label=mb)
    ax.set_xlabel("debris diameter (cm)"); ax.set_ylabel("USD per object removed")
    ax.set_title("Cost per object vs target size"); ax.legend(fontsize=7)
    fig.suptitle("SIM-9: cost is set by the areal-density penalty, so it "
                 "improves steeply with target size and guidance quality",
                 fontsize=10)
    save_fig(fig, "sim09_cost",
             "Cost per kg and per object removed. The mass penalty scales as "
             "sigma_miss^2/A_target, so cost per object falls with better "
             "guidance and cost per kg falls steeply with target size.")
    return df, concl


if __name__ == "__main__":
    run()
