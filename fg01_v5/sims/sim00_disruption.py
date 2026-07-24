"""
SIM-0 -- Catastrophic-disruption boundary for the 200 km-shift kick.

Redraws the v4 falsification gate with Delta-v_shift instead of
Delta-v_deorbit and shows where the design point sits relative to the NASA SBM
catastrophic threshold.  Also establishes a constraint the plan did not carry:
the *self-disposal* condition on intercept velocity.
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from fg01.constants import (ALT_SWEEP_KM, DH_NOMINAL, DH_RANGE, GAMMA,
                            ES_C, ES_C_RANGE, BETA_RANGE, MU_KM, R_E_KM)
from fg01.orbital import dv_shift, dv_direct_reentry, dv_hohmann_two_burn, v_circ
from fg01.io_utils import save_table, save_fig, save_json, banner, CB


def design_es(dv, v_rel, beta, gamma=GAMMA):
    """E_s,design = gamma * dv * v_rel / (2 beta)   [J/kg]"""
    return gamma * dv * v_rel / (2.0 * beta)


def v_max_safe(dv, beta, es_c=ES_C):
    """Margin-free safe intercept velocity: 2 E_s,c beta / dv   [m/s]"""
    return 2.0 * es_c * beta / dv


def run():
    banner("SIM-0  Catastrophic-disruption boundary (200 km-shift kick)")
    alts = np.array(ALT_SWEEP_KM, float) * 1e3
    dhs = [150e3, 200e3, 250e3]
    betas = [1.0, 1.5, 2.0, 2.5]
    escs = [30e3, 40e3, 50e3]

    # ---------------- validation first (AM: show it, don't assert it) -------
    val = {}
    # (a) dv_shift vs the first burn of a two-burn Hohmann to the same perigee
    err = []
    for h in alts:
        for dh in dhs:
            dv1, _, _ = dv_hohmann_two_burn(h, h - dh)
            err.append(abs(dv1 - dv_shift(h, dh)) / dv1)
    val["hohmann_first_burn_max_rel_err"] = float(np.max(err))

    # (b) algebraic identity  E_s = gamma*dv*v_rel/(2 beta)  from
    #     m_p = gamma*dv*m_t/(beta*v_rel) and E_s = 0.5*m_p*v_rel^2/m_t
    rng = np.random.default_rng(1)
    m_t = rng.uniform(1e-3, 1.0, 2000)
    dv_ = rng.uniform(20, 300, 2000)
    v_ = rng.uniform(200, 3000, 2000)
    b_ = rng.uniform(1.0, 2.5, 2000)
    m_p = GAMMA * dv_ * m_t / (b_ * v_)
    es_direct = 0.5 * m_p * v_ ** 2 / m_t
    es_formula = design_es(dv_, v_, b_)
    val["es_substitution_max_rel_err"] = float(
        np.max(np.abs(es_direct - es_formula) / es_formula))

    # (c) vis-viva reproduces the plan's Part 1.4 delta-v column
    plan_dv = {400: 57.7, 600: 55.2, 800: 52.9, 1000: 50.7, 1200: 48.7}
    val["plan_table_dv_shift_check"] = {
        str(k): {"plan_m_s": v, "computed_m_s": round(float(dv_shift(k * 1e3, 200e3)), 2)}
        for k, v in plan_dv.items()}
    plan_direct = {400: 87.4, 600: 142.0, 800: 193.8, 1000: 243.1, 1200: 290.0}
    val["plan_table_dv_direct_check"] = {
        str(k): {"plan_m_s": v, "computed_m_s": round(float(dv_direct_reentry(k * 1e3)), 2)}
        for k, v in plan_direct.items()}
    save_json(val, "sim00_validation")
    print("  validation:", {k: v for k, v in val.items() if "max_rel_err" in k})

    # ---------------- main sweep -------------------------------------------
    rows = []
    for h in alts:
        dvd = float(dv_direct_reentry(h))
        for dh in dhs:
            dv = float(dv_shift(h, dh))
            for beta in betas:
                for esc in escs:
                    vms = float(v_max_safe(dv, beta, esc))
                    rows.append(dict(
                        h_km=h / 1e3, dh_km=dh / 1e3, beta=beta,
                        Es_c_kJ_kg=esc / 1e3,
                        v_circ_m_s=float(v_circ(h)),
                        dv_shift_m_s=dv,
                        dv_direct_m_s=dvd,
                        dv_ratio_direct_over_shift=dvd / dv,
                        Es_design_619_kJ_kg=design_es(dv, 619.0, beta) / 1e3,
                        Es_design_direct_619_kJ_kg=design_es(dvd, 619.0, beta) / 1e3,
                        v_max_safe_m_s=vms,
                        v_design_max_m_s=vms / GAMMA,
                        v_self_disposal_min_m_s=dvd,
                        self_disposal_ok=bool(vms / GAMMA >= dvd),
                    ))
    df = pd.DataFrame(rows)
    save_table(df, "sim00_boundary_full",
               "SIM-0 disruption boundary: dv_shift, design specific energy, "
               "safe intercept-velocity ceiling and the self-disposal floor, "
               "for every (h, dh, beta, Es_c).",
               tags={"dv_shift": "DERIVED", "Es_c": "SOURCED(johnson2001)",
                     "beta": "UNVALIDATED", "gamma": "DERIVED"})

    # headline table reproducing the plan's Part 1.4 layout
    head = df[(df.dh_km == 200) & (df.beta == 1.5) & (df.Es_c_kJ_kg == 40)]
    head = head[["h_km", "dv_direct_m_s", "dv_shift_m_s",
                 "Es_design_direct_619_kJ_kg", "Es_design_619_kJ_kg",
                 "v_max_safe_m_s", "v_design_max_m_s", "v_self_disposal_min_m_s"]]
    save_table(head, "sim00_headline",
               "Design point at beta=1.5, gamma=2, Es_c=40 kJ/kg, dh=200 km, "
               "v_rel=619 m/s (plan Part 1.4 layout, recomputed).")
    print(head.to_string(index=False))

    # ---------------- conclusion numbers ------------------------------------
    sel = df[(df.dh_km == 200) & (df.Es_c_kJ_kg == 40)]
    concl = {
        "Es_design_kJ_kg_range_at_619": [float(sel.Es_design_619_kJ_kg.min()),
                                         float(sel.Es_design_619_kJ_kg.max())],
        "Es_design_v4_direct_kJ_kg_range_at_619": [
            float(sel.Es_design_direct_619_kJ_kg.min()),
            float(sel.Es_design_direct_619_kJ_kg.max())],
        "v_max_safe_range_m_s": [float(sel.v_max_safe_m_s.min()),
                                 float(sel.v_max_safe_m_s.max())],
        "v_design_max_range_m_s": [float(sel.v_design_max_m_s.min()),
                                   float(sel.v_design_max_m_s.max())],
        "fraction_of_threshold_at_619": [
            float(sel.Es_design_619_kJ_kg.min() / 40.0),
            float(sel.Es_design_619_kJ_kg.max() / 40.0)],
        "self_disposal_satisfied_everywhere_beta1_esc30": bool(
            df[(df.beta == 1.0) & (df.Es_c_kJ_kg == 30) & (df.dh_km == 200)]
            .self_disposal_ok.all()),
        "worst_case_v_design_max_vs_self_disposal": {
            str(int(r.h_km)): {"v_design_max_m_s": round(r.v_design_max_m_s, 1),
                               "v_self_disposal_min_m_s": round(r.v_self_disposal_min_m_s, 1)}
            for _, r in df[(df.beta == 1.0) & (df.Es_c_kJ_kg == 30)
                           & (df.dh_km == 200)].iterrows()},
    }
    save_json(concl, "sim00_conclusions")
    print("  Es_design (beta 1-2.5, 619 m/s):",
          [round(x, 1) for x in concl["Es_design_kJ_kg_range_at_619"]], "kJ/kg")
    print("  v4 direct-reentry equivalent:",
          [round(x, 1) for x in concl["Es_design_v4_direct_kJ_kg_range_at_619"]], "kJ/kg")

    # ---------------- figures ----------------------------------------------
    # Fig 1: v_max,safe(h) with beta as parameter, one panel per Es_c
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.9), sharey=True)
    for ax, esc in zip(axes, escs):
        for i, beta in enumerate(betas):
            s = df[(df.dh_km == 200) & (df.beta == beta) & (df.Es_c_kJ_kg == esc)]
            ax.plot(s.h_km, s.v_max_safe_m_s, color=CB[i], label=f"β={beta}")
            ax.plot(s.h_km, s.v_design_max_m_s, color=CB[i], ls="--", lw=1.0)
        s0 = df[(df.dh_km == 200) & (df.beta == 1.0) & (df.Es_c_kJ_kg == esc)]
        ax.plot(s0.h_km, s0.v_self_disposal_min_m_s, color="0.35", ls=":",
                lw=1.6, label="self-disposal floor")
        ax.axhspan(0, 200, color="0.9", zorder=0)
        ax.set_title(f"$E_{{s,c}}$ = {esc/1e3:.0f} kJ/kg")
        ax.set_xlabel("altitude (km)")
    axes[0].set_ylabel("intercept velocity $v_{rel}$ (m/s)")
    axes[0].set_ylim(0, 6200)
    axes[0].legend(fontsize=7.5, ncol=1)
    fig.suptitle("SIM-0: safe intercept-velocity ceiling for a 200 km shift "
                 "(solid = margin-free, dashed = with γ=2)", fontsize=10)
    save_fig(fig, "sim00_vmax_safe",
             "Safe intercept-velocity ceiling vs altitude. Everything below a "
             "solid line is sub-catastrophic; below the dashed line it stays "
             "sub-catastrophic with the γ=2 mass margin. The dotted line is "
             "the minimum v_rel that guarantees missed projectile mass is on "
             "an immediate-reentry trajectory.")

    # Fig 2: the headline pivot figure
    fig, ax = plt.subplots()
    for i, beta in enumerate([1.0, 1.5, 2.5]):
        s = df[(df.dh_km == 200) & (df.beta == beta) & (df.Es_c_kJ_kg == 40)]
        ax.plot(s.h_km, s.Es_design_direct_619_kJ_kg, color=CB[i], ls="--",
                label=f"v4 direct reentry, β={beta}")
        ax.plot(s.h_km, s.Es_design_619_kJ_kg, color=CB[i], ls="-",
                label=f"v5 200 km shift, β={beta}")
    ax.axhline(40, color="k", lw=1.4)
    ax.axhspan(40, 200, color="#D55E00", alpha=0.10)
    ax.text(410, 42, "catastrophic (NASA SBM, 40 kJ/kg)", fontsize=8)
    ax.set_xlabel("altitude (km)")
    ax.set_ylabel("design specific energy $E_s$ (kJ/kg)")
    ax.set_ylim(0, 130)
    ax.legend(fontsize=7.5, ncol=2)
    ax.set_title("SIM-0: the v4→v5 pivot, at the common design point "
                 "$v_{rel}$=619 m/s, γ=2")
    save_fig(fig, "sim00_pivot_headline",
             "Design-point specific energy for direct reentry (v4) versus the "
             "200 km accelerated-reentry shift (v5) at a common intercept "
             "velocity. The v5 design point is sub-catastrophic at every "
             "altitude and for every β in the envelope.")
    return df, concl


if __name__ == "__main__":
    run()
