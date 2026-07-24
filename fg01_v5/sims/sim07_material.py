"""
SIM-7 -- Material trade study: rhenium versus tungsten and alternatives.

The project is named for rhenium, so the honest answer matters more here than
anywhere else. Momentum p = mv is independent of density at fixed mass, so the
concept's founding justification for rhenium is physically void (R3). The one
axis on which density can legitimately matter is the orbital persistence of
missed material: at fixed grain mass, a denser grain is smaller and has a lower
area-to-mass ratio -- which makes it decay *slower*, not faster. That is the
opposite of the advantage claimed for it, and it is quantified here.
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from fg01.constants import (RHO_RE, RHO_W, RHO_TA, RHO_STEEL, P_RE_2026,
                            P_RE_2025, P_W_RANGE, P_TA, P_STEEL,
                            RE_PRODUCTION_TPY, W_PRODUCTION_TPY,
                            H_ABL_RE, H_ABL_W, H_ABL_AL, C_LAUNCH_RANGE)
from fg01.orbital import area_to_mass_sphere
from fg01.io_utils import save_table, save_fig, save_json, banner, CB

MATERIALS = {
    "Rhenium":  dict(rho=RHO_RE, price=(P_RE_2025, P_RE_2026), prod_tpy=RE_PRODUCTION_TPY,
                     oxide="Re2O7, volatile (m.p. 297 C, sublimes)", h_abl=H_ABL_RE,
                     tox="Re2O7 is volatile and moderately toxic; disperses as vapour"),
    "Tungsten": dict(rho=RHO_W, price=P_W_RANGE, prod_tpy=(W_PRODUCTION_TPY, W_PRODUCTION_TPY),
                     oxide="WO3, refractory", h_abl=H_ABL_W,
                     tox="WO3 low volatility, low acute toxicity"),
    "Tantalum": dict(rho=RHO_TA, price=(P_TA, P_TA * 2), prod_tpy=(2000.0, 2000.0),
                     oxide="Ta2O5, inert", h_abl=7.0e6, tox="Ta2O5 inert"),
    "Steel":    dict(rho=RHO_STEEL, price=(P_STEEL, 3.0), prod_tpy=(1.9e9, 1.9e9),
                     oxide="Fe oxides", h_abl=6.0e6, tox="benign"),
}


def run():
    banner("SIM-7  Material trade study (rhenium vs tungsten)")
    m_grain = 1e-6            # 1 mg reference grain
    rows = []
    for name, m in MATERIALS.items():
        d = 2 * (3 * m_grain / (4 * np.pi * m["rho"])) ** (1 / 3)
        aom = float(area_to_mass_sphere(d, m["rho"]))
        # momentum per kg is identical by construction; show it explicitly
        p_per_kg = 619.0
        c_lo = m["price"][0] + C_LAUNCH_RANGE[0]
        c_hi = m["price"][1] + C_LAUNCH_RANGE[1]
        rows.append(dict(
            material=name, density_kg_m3=m["rho"],
            momentum_per_kg_at_619_Ns=p_per_kg,
            grain_diameter_1mg_mm=d * 1e3,
            A_over_m_1mg_m2_kg=aom,
            price_usd_kg_low=m["price"][0], price_usd_kg_high=m["price"][1],
            production_t_yr_low=m["prod_tpy"][0], production_t_yr_high=m["prod_tpy"][1],
            delivered_cost_usd_kg_low=c_lo, delivered_cost_usd_kg_high=c_hi,
            H_ablation_MJ_kg=m["h_abl"] / 1e6,
            oxide=m["oxide"], toxicity=m["tox"]))
    df = pd.DataFrame(rows)
    save_table(df, "sim07_material_comparison",
               "Material comparison on the metrics that actually differentiate: "
               "momentum per unit mass (identical), delivered cost, supply, "
               "missed-grain area-to-mass, ablation energy, oxide behaviour.",
               tags={"prices": "SOURCED, dated", "H_abl": "SOURCED range"})
    print(df[["material", "density_kg_m3", "grain_diameter_1mg_mm",
              "A_over_m_1mg_m2_kg", "price_usd_kg_low", "price_usd_kg_high",
              "production_t_yr_high"]].to_string(index=False))

    # supply-ceiling analysis: how many engagements per year can each material
    # support at the SIM-4 design shot mass?
    sup = []
    for shot_g in (2.69, 175.0):
        for name, m in MATERIALS.items():
            prod = m["prod_tpy"][1] * 1000.0        # kg/yr
            # assume at most 1% of world production is available to the programme
            avail = 0.01 * prod
            sup.append(dict(shot_mass_g=shot_g, material=name,
                            world_production_kg_yr=prod,
                            available_at_1pct_kg_yr=avail,
                            shots_per_year=avail / (shot_g * 1e-3),
                            objects_per_year=avail / (shot_g * 1e-3)))
    dfs = pd.DataFrame(sup)
    save_table(dfs, "sim07_supply_ceiling",
               "Engagements per year supportable at 1% of world production.")

    weights = {"delivered_cost": 0.35, "supply": 0.25, "momentum": 0.0,
               "missed_persistence": 0.20, "environmental": 0.20}
    score_rows = []
    for _, r in df.iterrows():
        cost_s = 1.0 / np.log10(r.delivered_cost_usd_kg_high)
        sup_s = np.log10(r.production_t_yr_high) / 10.0
        # lower A/m = slower decay of missed grains = worse
        pers_s = r.A_over_m_1mg_m2_kg / df.A_over_m_1mg_m2_kg.max()
        env_s = {"Rhenium": 0.5, "Tungsten": 0.8, "Tantalum": 0.9, "Steel": 1.0}[r.material]
        total = (weights["delivered_cost"] * cost_s / 0.25
                 + weights["supply"] * sup_s
                 + weights["missed_persistence"] * pers_s
                 + weights["environmental"] * env_s)
        score_rows.append(dict(material=r.material, cost_score=cost_s,
                               supply_score=sup_s, persistence_score=pers_s,
                               environmental_score=env_s, weighted_total=total))
    dsc = pd.DataFrame(score_rows).sort_values("weighted_total", ascending=False)
    save_table(dsc, "sim07_weighted_ranking",
               "Weighted ranking. Weights stated explicitly: delivered cost "
               "0.35, supply 0.25, missed-grain persistence 0.20, "
               "environmental 0.20, momentum 0.00 (identical across materials).")
    print("\n", dsc.to_string(index=False))

    re_ = df[df.material == "Rhenium"].iloc[0]
    w_ = df[df.material == "Tungsten"].iloc[0]
    concl = {
        "momentum_finding": (
            "At fixed projectile mass, momentum delivered is identical for all "
            "materials: p = mv contains no density term. The concept's founding "
            "justification for rhenium -- that it is the fourth-densest stable "
            "element -- has no bearing on the momentum budget."),
        "density_effect_on_missed_grains": {
            "grain_diameter_1mg_mm": {"Re": float(re_.grain_diameter_1mg_mm),
                                      "W": float(w_.grain_diameter_1mg_mm)},
            "A_over_m_1mg": {"Re": float(re_.A_over_m_1mg_m2_kg),
                             "W": float(w_.A_over_m_1mg_m2_kg)},
            "finding": (
                "A denser grain of the same mass is smaller and has a LOWER "
                "area-to-mass ratio, so it decays SLOWER, not faster. Rhenium's "
                f"A/m is {(1 - re_.A_over_m_1mg_m2_kg / w_.A_over_m_1mg_m2_kg)*100:.1f}% "
                "below tungsten's at equal grain mass. The one "
                "axis on which density could have favoured rhenium therefore "
                "runs slightly against it. In this design it is moot: SIM-2 "
                "shows missed mass is on an immediate-reentry trajectory "
                "regardless of material, so persistence is zero either way.")},
        "cost_ratio_Re_over_W_material": float(P_RE_2026 / P_W_RANGE[1]),
        "supply_ratio_W_over_Re": float(W_PRODUCTION_TPY / RE_PRODUCTION_TPY[1]),
        "delivered_cost_ratio_Re_over_W": float(
            re_.delivered_cost_usd_kg_high / w_.delivered_cost_usd_kg_high),
        "verdict": (
            "Tungsten dominates rhenium on every metric that differentiates "
            "them: 73-243x cheaper as raw material, 1037x more available, "
            "chemically better behaved on reentry, and marginally better on "
            "missed-grain decay. Rhenium's only remaining advantages are a 9% "
            "density edge, which buys a 3% smaller grain and nothing else, and "
            "a higher melting point that is irrelevant at a 619 m/s impact. "
            "The honest recommendation is tungsten. Because launch cost "
            "dominates delivered cost at the small per-shot masses v5 needs, "
            "the switch changes total cost by a factor of ~1.6, not ~100 -- so "
            "the material choice is not what decides viability, but rhenium "
            "cannot be defended on physics."),
        "if_rhenium_retained": (
            "If the programme retains rhenium for non-technical reasons, the "
            "supply ceiling is the binding issue: at 1% of world production "
            "(0.81 t/yr) and the 2.69 g design shot, rhenium supports ~300,000 "
            "engagements per year, which is not limiting; at the 175 g shot it "
            "supports ~4,600, which is."),
    }
    save_json(concl, "sim07_conclusions")

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.0))
    ax = axes[0]
    x = np.arange(len(df))
    ax.bar(x - 0.2, df.price_usd_kg_high, 0.4, color=CB[0], label="material")
    ax.bar(x + 0.2, df.delivered_cost_usd_kg_high, 0.4, color=CB[1],
           label="material + launch")
    ax.set_yscale("log"); ax.set_xticks(x); ax.set_xticklabels(df.material, fontsize=8)
    ax.set_ylabel("USD/kg (upper bound)"); ax.legend(fontsize=8)
    ax.set_title("Delivered cost")
    ax = axes[1]
    ax.bar(x, df.production_t_yr_high, 0.5, color=CB[2])
    ax.set_yscale("log"); ax.set_xticks(x); ax.set_xticklabels(df.material, fontsize=8)
    ax.set_ylabel("world production (t/yr)")
    ax.set_title("Supply")
    fig.suptitle("SIM-7: rhenium is dominated by tungsten on both axes that "
                 "differentiate materials", fontsize=10)
    save_fig(fig, "sim07_material_trade",
             "Material cost and supply. Because momentum transfer is "
             "density-independent, these are the axes that actually "
             "discriminate, and rhenium loses both.")
    return df, concl


if __name__ == "__main__":
    run()
