"""
SIM-1 -- Momentum transfer and mass scaling f(h), no assumed efficiency (AM-1).

m_p,design = gamma * dv_shift(h) * m_t / (beta * v_rel)
f(h)       = m_p,design / m_t = gamma * dv_shift(h) / (beta * v_rel)   [target-mass free]

This is the *ideal* projectile mass -- the mass that must physically strike the
target.  The mass that must be *launched* is larger by the cloud/target area
ratio and is computed in SIM-4; the two must never be confused, so this module
reports m_p,design as "mass on target" throughout.
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from fg01.constants import (ALT_SWEEP_KM, GAMMA, ES_C, RHO_AL, RHO_RE, RHO_W,
                            DEBRIS_DIAM_M)
from fg01.orbital import dv_shift, dv_direct_reentry, mass_sphere
from fg01.io_utils import save_table, save_fig, save_json, banner, CB
from fg01.sbm import specific_energy

LIMITATION = ("AM-15/AM-16: treating N micro-grain impacts as one macroscopic "
              "impactor of equal total mass is UNVALIDATED at these scales. "
              "Every mass figure here is a momentum-balance result and is "
              "independent of how the mass is divided; the *disruption* "
              "consequence of that division is bounded both ways in SIM-5.")


def f_of_h(h, v_rel, beta, dh=200e3, gamma=GAMMA):
    return gamma * dv_shift(h, dh) / (beta * v_rel)


def run():
    banner("SIM-1  Momentum transfer and mass scaling f(h)")
    alts = np.array(ALT_SWEEP_KM, float) * 1e3
    betas = [1.0, 1.5, 2.0, 2.5]
    v_rels = [400.0, 619.0, 800.0, 1000.0]

    rows = []
    for h in alts:
        dv = float(dv_shift(h, 200e3))
        dvd = float(dv_direct_reentry(h))
        for beta in betas:
            for v in v_rels:
                f = gamma_f = GAMMA * dv / (beta * v)
                f_min = dv / (beta * v)
                f_v4 = GAMMA * dvd / (beta * v)
                for d in DEBRIS_DIAM_M:
                    m_t = float(mass_sphere(d, RHO_AL))
                    m_p_des = f * m_t
                    rows.append(dict(
                        h_km=h / 1e3, beta=beta, v_rel_m_s=v,
                        d_debris_cm=d * 100, m_t_g=m_t * 1e3,
                        m_p_min_g=f_min * m_t * 1e3,
                        m_p_design_g=m_p_des * 1e3,
                        f_h_g_per_g=f,
                        f_h_v4_direct_g_per_g=f_v4,
                        f_ratio_v4_over_v5=f_v4 / f,
                        Es_kJ_kg=float(specific_energy(m_p_des, v, m_t)) / 1e3,
                        Es_below_threshold=bool(
                            specific_energy(m_p_des, v, m_t) < ES_C),
                        d_grain_Re_mm_if_single=1e3 * 2 * (
                            3 * m_p_des / (4 * np.pi * RHO_RE)) ** (1 / 3),
                        d_grain_W_mm_if_single=1e3 * 2 * (
                            3 * m_p_des / (4 * np.pi * RHO_W)) ** (1 / 3),
                    ))
    df = pd.DataFrame(rows)
    save_table(df, "sim01_mass_scaling_full",
               "Momentum-balance projectile mass on target, f(h), and the "
               "consistency check E_s < E_s,c, over (h, beta, v_rel, d_debris).",
               tags={"m_p": "DERIVED", "beta": "UNVALIDATED",
                     "rho_Al": "SOURCED", "limitation": "AM-15/AM-16"})

    # design-point table at the plan's reference conditions
    ref = df[(df.v_rel_m_s == 619.0) & (df.d_debris_cm == 1.0)][
        ["h_km", "beta", "m_t_g", "m_p_min_g", "m_p_design_g", "f_h_g_per_g",
         "f_ratio_v4_over_v5", "Es_kJ_kg", "Es_below_threshold"]]
    save_table(ref, "sim01_design_point",
               "1 cm Al-6061 sphere (AM-2), v_rel = 619 m/s, dh = 200 km, gamma = 2.")
    print(ref[ref.beta.isin([1.0, 1.5, 2.5])].to_string(index=False))

    sel = df[(df.v_rel_m_s == 619.0) & (df.h_km == 800) & (df.d_debris_cm == 1.0)]
    concl = {
        "limitation": LIMITATION,
        "f_h_range_g_per_g_800km_619ms": [float(sel.f_h_g_per_g.min()),
                                          float(sel.f_h_g_per_g.max())],
        "Re_on_target_g_per_kg_Al_800km_619ms": [
            float(sel.f_h_g_per_g.min() * 1000), float(sel.f_h_g_per_g.max() * 1000)],
        "v5_mass_advantage_vs_v4_direct": {
            str(int(h)): round(float(dv_direct_reentry(h * 1e3) / dv_shift(h * 1e3, 200e3)), 2)
            for h in ALT_SWEEP_KM},
        "all_design_points_subcatastrophic_at_619": bool(
            df[df.v_rel_m_s == 619.0].Es_below_threshold.all()),
        "target_mass_independence_check": float(
            np.ptp(df[(df.h_km == 800) & (df.beta == 1.5)
                      & (df.v_rel_m_s == 619)].f_h_g_per_g)),
    }
    save_json(concl, "sim01_conclusions")
    print("  f(h) at 800 km, 619 m/s, β∈[1,2.5]:",
          [round(x, 4) for x in concl["f_h_range_g_per_g_800km_619ms"]], "g Re per g Al")
    print("  v5 needs this many times less mass than v4:",
          concl["v5_mass_advantage_vs_v4_direct"])

    # ---- figures -----------------------------------------------------------
    fig, ax = plt.subplots()
    hk = np.array(ALT_SWEEP_KM, float)
    for i, v in enumerate([400.0, 619.0, 1000.0]):
        lo = np.array([f_of_h(h * 1e3, v, 2.5) for h in hk])
        hi = np.array([f_of_h(h * 1e3, v, 1.0) for h in hk])
        ax.fill_between(hk, lo, hi, color=CB[i], alpha=0.22)
        ax.plot(hk, np.array([f_of_h(h * 1e3, v, 1.5) for h in hk]),
                color=CB[i], label=f"v5, 200 km shift, $v_{{rel}}$={v:.0f} m/s")
    lo4 = np.array([GAMMA * dv_direct_reentry(h * 1e3) / (2.5 * 619.0) for h in hk])
    hi4 = np.array([GAMMA * dv_direct_reentry(h * 1e3) / (1.0 * 619.0) for h in hk])
    ax.fill_between(hk, lo4, hi4, color="0.5", alpha=0.25)
    ax.plot(hk, np.array([GAMMA * dv_direct_reentry(h * 1e3) / (1.5 * 619.0) for h in hk]),
            color="0.3", ls="--", label="v4, direct reentry, 619 m/s")
    ax.set_yscale("log")
    ax.set_xlabel("altitude (km)")
    ax.set_ylabel("f(h) = mass on target per unit debris mass (g/g)")
    ax.set_title("SIM-1: rhenium mass on target, bands span β ∈ [1.0, 2.5]")
    ax.legend(fontsize=7.5)
    save_fig(fig, "sim01_f_of_h",
             "Required projectile mass on target per unit debris mass. Bands "
             "span the unvalidated β envelope. v5 needs 1.5–6× less mass than "
             "v4 at the same intercept velocity, and the advantage grows with "
             "altitude because dv_shift falls while dv_direct rises.")

    fig, ax = plt.subplots()
    for i, beta in enumerate([1.0, 1.5, 2.5]):
        s = df[(df.beta == beta) & (df.v_rel_m_s == 619) & (df.d_debris_cm == 1.0)]
        ax.plot(s.h_km, s.d_grain_Re_mm_if_single, color=CB[i],
                label=f"rhenium, β={beta}")
        ax.plot(s.h_km, s.d_grain_W_mm_if_single, color=CB[i], ls="--",
                label=f"tungsten, β={beta}")
    ax.set_xlabel("altitude (km)")
    ax.set_ylabel("equivalent single-slug diameter (mm)")
    ax.set_title("SIM-1: equivalent compact-projectile size, 1 cm Al target, "
                 "$v_{rel}$=619 m/s")
    ax.legend(fontsize=7.5, ncol=2)
    save_fig(fig, "sim01_grain_size",
             "Diameter of a single compact projectile carrying the full design "
             "mass. Rhenium's 9% density advantage over tungsten changes the "
             "diameter by only 3%.")
    return df, concl


if __name__ == "__main__":
    run()
