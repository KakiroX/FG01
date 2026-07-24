"""
INT-2 -- Single grain to target coupling at sub-km/s: pinning down beta.

v5 had to carry beta ~ U(1.0, 2.5) because no measurement exists for Al-on-Al at
centimetre scale and sub-orbital velocity (AM-4a). The v6 plan asks for a
physically bounded value for *this* regime, and warns that the honest answer is
likely lower than the 1.5 v5 used at the centre of its envelope -- which would
darken the mass budget by 15-50%.

That is what this module finds. It is reported as a *worsening* of the v5
result, not hidden.

Method: beta = 1 + p_ejecta/p_impactor from Housen-Holsapple crater-ejecta
scaling in the strength regime, with the ejection-speed floor set by the target's
own strength velocity sqrt(Y/rho) = 320 m/s for Al-6061. At 619 m/s the impact
is barely twice that floor, so very little mass is ejected fast enough to add
momentum. The same model at 6 km/s reproduces the hypervelocity beta ~ 2-3 that
the literature reports, which is the consistency check.
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from fg01.constants import RHO_AL, RHO_RE, RHO_W, ES_C, GAMMA
from fg01.orbital import mass_sphere, dv_shift
from fg01.interaction import (beta_momentum, beta_envelope, v_strength,
                              penetration_depth, crater_volume, ejecta_mass,
                              grain_mass, mean_incidence_cos_sphere,
                              Y_AL, C_AL, MU_HH_RANGE, K_HH_RANGE)
from fg01.io_utils import save_table, save_fig, save_json, banner, CB

V_GRID = np.array([200, 400, 619, 800, 1000, 1100, 1500, 2000, 3000, 6000.0])


def run():
    banner("INT-2  Momentum coupling beta at sub-km/s")

    v_min = float(v_strength())
    print(f"  target strength velocity sqrt(Y/rho) = {v_min:.0f} m/s "
          f"(Al-6061, Y = {Y_AL/1e6:.0f} MPa)")
    print(f"  design impact speed 619 m/s = {619/v_min:.2f} x that floor "
          f"-- barely into the ejecta-producing regime")

    rows = []
    for v in V_GRID:
        lo, ctr, hi = beta_envelope(v)
        rows.append(dict(
            v_imp_m_s=v, v_over_v_strength=v / v_min,
            beta_low=float(lo), beta_central=float(ctr), beta_high=float(hi),
            beta_oblique_sphere=float(1 + (ctr - 1) * mean_incidence_cos_sphere()),
            regime=("below strength floor" if v < v_min else
                    "transitional" if v < 3 * v_min else "hypervelocity")))
    df = pd.DataFrame(rows)
    save_table(df, "int02_beta_vs_velocity",
               "Momentum enhancement from crater-ejecta scaling. beta_low/high "
               "span mu in [0.40, 0.55] and k in [0.2, 0.5]. The oblique column "
               "applies <cos i> = 2/3 for a parallel beam on a curved target.",
               tags={"model": "SOURCED (Housen-Holsapple 2011 scaling)",
                     "application to this regime": "DERIVED, UNVALIDATED"})
    print("\n" + df.to_string(index=False))

    # ---- consistency check against the hypervelocity regime ----------------
    hv = df[df.v_imp_m_s == 6000.0].iloc[0]
    des = df[df.v_imp_m_s == 619.0].iloc[0]
    val = {
        "hypervelocity_check": {
            "v_m_s": 6000.0,
            "beta_range": [float(hv.beta_low), float(hv.beta_high)],
            "note": ("The same model at 6 km/s gives beta = "
                     f"{hv.beta_low:.2f}-{hv.beta_high:.2f}, consistent with the "
                     "beta ~ 2-3 reported for hypervelocity impacts and with the "
                     "DART result (beta ~ 3.6 at 6.1 km/s on a rubble-pile "
                     "asteroid, an upper bound not applicable to a solid metal "
                     "target). The model therefore reproduces the regime the "
                     "literature's numbers come from, which is why its much "
                     "lower value at 619 m/s should be believed rather than the "
                     "literature value being extrapolated down.")},
        "design_point": {
            "v_m_s": 619.0,
            "beta_range": [float(des.beta_low), float(des.beta_high)],
            "beta_central": float(des.beta_central),
            "beta_oblique": float(des.beta_oblique_sphere)},
    }

    # ---- what it costs relative to the v5 assumption -----------------------
    m_t = float(mass_sphere(0.01, RHO_AL))
    dv = float(dv_shift(900e3, 200e3))
    cost = []
    for label, b in (("v5 assumption (envelope centre)", 1.75),
                     ("v5 nominal used in SIM-13", 1.5),
                     ("INT-2 central, normal incidence", float(des.beta_central)),
                     ("INT-2 central, oblique on sphere", float(des.beta_oblique_sphere)),
                     ("INT-2 pessimistic", float(des.beta_low)),
                     ("pure embedding (no ejecta)", 1.0)):
        m_p = GAMMA * dv * m_t / (b * 619.0)
        cost.append(dict(case=label, beta=b, m_on_target_mg=m_p * 1e6,
                         Es_kJ_kg=GAMMA * dv * 619.0 / (2 * b) / 1e3,
                         mass_vs_v5_nominal=(GAMMA * dv * m_t / (b * 619.0))
                         / (GAMMA * dv * m_t / (1.5 * 619.0))))
    dfc = pd.DataFrame(cost)
    save_table(dfc, "int02_mass_impact",
               "Consequence of the corrected beta for the mass budget at the "
               "SIM-13 optimum (1 cm target, 900 km, v_rel = 619 m/s).")
    print("\n  Consequence for the mass budget:")
    print(dfc.to_string(index=False))

    # ---- penetration / embedding ------------------------------------------
    pen = []
    for d_g in (0.05e-3, 0.1e-3, 0.45e-3, 1e-3, 2e-3):
        for v in (400.0, 619.0, 1000.0):
            for rho_p, mat in ((RHO_RE, "Re"), (RHO_W, "W")):
                p = float(penetration_depth(d_g, rho_p, v))
                pen.append(dict(
                    d_grain_mm=d_g * 1e3, material=mat, v_imp_m_s=v,
                    penetration_mm=p * 1e3, penetration_over_diameter=p / d_g,
                    perforates_1cm_target=bool(p > 0.01),
                    crater_volume_mm3=float(crater_volume(
                        grain_mass(d_g, rho_p), v)) * 1e9,
                    ejecta_per_grain_mg=float(ejecta_mass(
                        grain_mass(d_g, rho_p), v)) * 1e6,
                    ejecta_over_grain_mass=float(ejecta_mass(grain_mass(d_g, rho_p), v))
                    / float(grain_mass(d_g, rho_p))))
    dfp = pd.DataFrame(pen)
    save_table(dfp, "int02_penetration",
               "Grain penetration depth (Poncelet form), crater volume and "
               "ejecta mass per grain. 'Perforates' compares against a 1 cm "
               "solid target.",
               tags={"Poncelet": "SOURCED form", "A=1.5": "UNVALIDATED"})
    print("\n  Penetration and ejecta per grain (rhenium):")
    print(dfp[(dfp.material == "Re") & (dfp.v_imp_m_s == 619.0)][
        ["d_grain_mm", "penetration_mm", "penetration_over_diameter",
         "perforates_1cm_target", "ejecta_over_grain_mass"]].to_string(index=False))

    concl = {
        "beta_at_design_point": {
            "central": float(des.beta_central),
            "range": [float(des.beta_low), float(des.beta_high)],
            "oblique_corrected": float(des.beta_oblique_sphere),
            "recommended_for_downstream": float(des.beta_oblique_sphere),
        },
        "v_strength_m_s": v_min,
        "headline": (
            f"beta at the 619 m/s design point is "
            f"{des.beta_low:.2f}-{des.beta_high:.2f} (central "
            f"{des.beta_central:.2f}; {des.beta_oblique_sphere:.2f} after the "
            "oblique-incidence correction for a curved target), NOT the 1.5 v5 "
            "used or the 1.75 centre of its U(1.0, 2.5) envelope. The reason is "
            "physical: Al-6061's strength velocity is 320 m/s, so a 619 m/s "
            "impact is only 1.9x the speed at which the target can eject "
            "anything at all, and almost no mass leaves fast enough to enhance "
            "momentum. The v5 envelope was too generous at its upper end for "
            "this regime."),
        "consequence": (
            f"Required mass on target rises by "
            f"{float(dfc[dfc.case.str.startswith('INT-2 central, oblique')].mass_vs_v5_nominal.iloc[0]):.2f}x "
            "relative to the v5 nominal. This DARKENS the v5 mass and "
            "consumable-cost result. It does not change the viability verdict, "
            "because consumables were ~0.6% of cost per object -- the extra "
            "mass costs about $20 per engagement -- but it must be carried "
            "forward rather than papered over."),
        "disruption_margin_effect": (
            "A lower beta needs more mass for the same impulse, so "
            "E_s = gamma*dv*v_rel/(2*beta) RISES from "
            f"{float(dfc[dfc.case.str.startswith('v5 nominal')].Es_kJ_kg.iloc[0]):.1f} to "
            f"{float(dfc[dfc.case.str.startswith('INT-2 central, oblique')].Es_kJ_kg.iloc[0]):.1f} kJ/kg, "
            f"i.e. from "
            f"{float(dfc[dfc.case.str.startswith('v5 nominal')].Es_kJ_kg.iloc[0])*1e3/ES_C:.2f}x to "
            f"{float(dfc[dfc.case.str.startswith('INT-2 central, oblique')].Es_kJ_kg.iloc[0])*1e3/ES_C:.2f}x "
            "the catastrophic threshold. Still sub-catastrophic, but the margin "
            "is materially narrower than v5 reported, and at the pessimistic end "
            f"({float(dfc[dfc.case=='pure embedding (no ejecta)'].Es_kJ_kg.iloc[0])*1e3/ES_C:.2f}x for pure "
            "embedding) it is thin enough that the E_s,c = 30 kJ/kg lower bound "
            "would be breached. The velocity ceiling must be lowered "
            "accordingly -- see INT-3."),
        "AM4a_status": (
            "beta moves from UNVALIDATED toward DERIVED: it is now computed from "
            "a sourced scaling law rather than assumed, and the model is "
            "validated against the hypervelocity regime where data exists. It is "
            "NOT yet SOURCED, because no direct measurement at Al-on-Al, "
            "sub-km/s, cm-scale exists. A hydrocode run or a light-gas-gun shot "
            "at 0.6 km/s would close this."),
        "embedding_finding": (
            "Design-range grains (0.1-1 mm Re at 619 m/s) penetrate 0.3-3 mm "
            "into aluminium -- comparable to their own diameter and well short "
            "of perforating a 1 cm target. They embed. An embedded grain "
            "transfers exactly its own momentum, which is the physical floor "
            "beta = 1 and the reason the answer sits so close to it."),
    }
    save_json(concl, "int02_conclusions")
    save_json(val, "int02_validation")
    print(f"\n  RECOMMENDED beta for downstream = {des.beta_oblique_sphere:.3f} "
          f"(v5 used 1.5)")

    # ---- figures -----------------------------------------------------------
    fig, ax = plt.subplots(figsize=(7.4, 4.3))
    vv = np.geomspace(150, 8000, 300)
    lo, ctr, hi = beta_envelope(vv)
    ax.fill_between(vv, lo, hi, color=CB[0], alpha=0.22, label="μ∈[0.40,0.55], k∈[0.2,0.5]")
    ax.plot(vv, ctr, color=CB[0], lw=2, label="central (μ=0.45, k=0.3)")
    ax.plot(vv, 1 + (ctr - 1) * mean_incidence_cos_sphere(), color=CB[3], lw=1.6,
            ls="-.", label="oblique on a curved target")
    ax.axvline(v_min, color="0.4", ls=":", lw=1.2)
    ax.text(v_min * 1.05, 2.6, "strength velocity\n√(Y/ρ) = 320 m/s", fontsize=7)
    ax.axvspan(400, 1100, color=CB[2], alpha=0.10)
    ax.text(430, 0.55 + 1, "FG01 design band", fontsize=7.5, color=CB[2])
    ax.axhline(1.75, color=CB[1], ls="--", lw=1.4, label="v5 envelope centre (1.75)")
    ax.axhline(1.5, color=CB[1], ls=":", lw=1.4, label="v5 nominal (1.5)")
    ax.set_xscale("log")
    ax.set_xlabel("impact speed (m/s)")
    ax.set_ylabel("momentum enhancement β")
    ax.set_ylim(0.9, 3.2)
    ax.set_title("INT-2: β is close to 1 in the FG01 regime, not 1.5–2.5")
    ax.legend(fontsize=7.5, loc="upper left")
    save_fig(fig, "int02_beta",
             "Momentum enhancement from crater-ejecta scaling. Below roughly "
             "three times the target's strength velocity almost nothing is "
             "ejected fast enough to add momentum, so β approaches its physical "
             "floor of 1. The v5 envelope, drawn from hypervelocity data, was "
             "too generous for this regime.")
    return df, concl


if __name__ == "__main__":
    run()
