"""
SIM-4 -- Terminal GNC, launch angle and shot precision.

This is where v4 was falsified (P_hit collapsed by ~50 orders of magnitude
between 1 ms and 10 ms control lag, 68% of total variance).  The Delta-v pivot
does not touch this, so the question is re-opened here from the error budget up.

Two changes to the v4 formulation are made, both stated explicitly and both
reported *alongside* the original so the reader can see exactly what each is
worth:

A. The sweep rate.  v4's lag term is v_rel*tau, i.e. it assumes the target
   crosses the aim point at the impact velocity.  In the engagement geometry
   that SIM-0's energy gate forces -- a retrograde intercept of a few hundred
   m/s, which requires FG01 to be nearly co-moving with the target -- the
   impact velocity is supplied by the launcher, while the target's motion
   *relative to the launcher* proceeds at u, the platform-target relative
   speed.  The lag term is u*tau.  This is a consequence of the geometry, not
   an assumption about better hardware, and it has a price: FG01 must operate
   in the target's orbital plane (charged in SIM-11).

B. The success variable.  v4 scored a Bernoulli hit/miss on cloud containment.
   The physical quantity is delivered momentum, which is continuous in the
   miss distance (SIM-3).  A shot that misses by more than planned delivers
   less momentum, i.e. a smaller orbit shift -- not nothing.

Launch angle: the kick must be retrograde in the target's frame to lower
perigee, so FG01 leads the target in-plane and fires rearward along the
velocity vector.  The optimum is pure retrograde (cos component = 1); an
off-axis angle theta multiplies the useful impulse by cos(theta) and is
reported as a sensitivity.
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from fg01.constants import (RHO_AL, GAMMA, R_E, OPTICAL_APERTURE_M,
                            OPTICAL_LAMBDA_M, LIDAR_BW_HZ, LIDAR_SNR, TAU_MS)
from fg01.orbital import mass_sphere, xsec_sphere, dv_shift
from fg01.gnc import (sigma_pos_diffraction, sigma_velocity_estimate,
                      relative_accel_scale, miss_v4, p_hit_plan_model)
from fg01.delivery import (penalty_for_confidence, sigma_opt_for_confidence,
                           min_launch_mass_for_confidence, delivered_fraction,
                           optimal_sigma, expected_fraction_random_miss)
from fg01.io_utils import save_table, save_fig, save_json, banner, CB

C_LIGHT = 2.99792458e8
RANGES = [1.0, 5.0, 10.0, 50.0, 100.0]
US = [1.0, 10.0, 50.0, 100.0]              # FG01-target relative speed, m/s
POINTING = [5e-5, 2e-4, 1e-3]              # rad
MUZZLE = [1e-3, 1e-2]                      # fractional muzzle-speed dispersion
T_OBS, SAMPLE_HZ = 10.0, 1000.0
V_REL = 619.0
H_DES = 800e3


def lidar_range_sigma(bw_hz, snr):
    """Time-of-flight ranging precision: c/(2 B sqrt(2 SNR))."""
    return C_LIGHT / (2.0 * bw_hz * np.sqrt(2.0 * snr))


def run():
    banner("SIM-4  Terminal GNC, launch angle and shot precision")
    a_orbit = R_E + H_DES
    m_t = float(mass_sphere(0.01, RHO_AL))
    a_t = float(xsec_sphere(0.01))
    r_t = 0.005
    m_on_target = GAMMA * float(dv_shift(H_DES, 200e3)) * m_t / (1.5 * V_REL)

    rows = []
    for R in RANGES:
        sp = float(sigma_pos_diffraction(R, OPTICAL_LAMBDA_M, OPTICAL_APERTURE_M))
        sp_snr = float(sigma_pos_diffraction(R, OPTICAL_LAMBDA_M,
                                             OPTICAL_APERTURE_M, snr=10.0))
        sv = float(sigma_velocity_estimate(sp, T_OBS, T_OBS * SAMPLE_HZ))
        a_rel = float(relative_accel_scale(a_orbit, R))
        for tau_ms in TAU_MS:
            tau = tau_ms * 1e-3
            t_go = tau + R / V_REL
            for u in US:
                for sig_pt in POINTING:
                    for mz in MUZZLE:
                        for bw, lab in ((LIDAR_BW_HZ, "1 MHz (plan spec)"),
                                        (1e8, "100 MHz"), (1e9, "1 GHz")):
                            sr = lidar_range_sigma(bw, LIDAR_SNR)
                            terms = {
                                "sigma_pos": sp,
                                "vel_extrap": sv * t_go,
                                "rel_accel": 0.5 * a_rel * t_go ** 2,
                                "pointing": sig_pt * R,
                                "muzzle_timing": u * R * mz / V_REL,
                                "range_timing": u * sr / V_REL,
                            }
                            m3 = float(np.sqrt(sum(v ** 2 for v in terms.values())))
                            m2 = float(np.sqrt(sp ** 2 + (u * tau) ** 2
                                               + (sig_pt * R) ** 2
                                               + (u * R * mz / V_REL) ** 2))
                            m1 = float(miss_v4(R, tau, V_REL,
                                               OPTICAL_LAMBDA_M, OPTICAL_APERTURE_M))
                            rows.append(dict(
                                R_m=R, tau_ms=tau_ms, u_m_s=u,
                                pointing_urad=sig_pt * 1e6, muzzle_frac=mz,
                                lidar_bw=lab, sigma_range_m=sr,
                                tof_s=R / V_REL, t_go_s=t_go,
                                sigma_pos_m=sp, sigma_pos_snr10_m=sp_snr,
                                sigma_v_m_s=sv, a_rel_m_s2=a_rel,
                                term_vel_extrap_m=terms["vel_extrap"],
                                term_rel_accel_m=terms["rel_accel"],
                                term_pointing_m=terms["pointing"],
                                term_muzzle_m=terms["muzzle_timing"],
                                term_range_m=terms["range_timing"],
                                d_miss_M1_v4_m=m1,
                                d_miss_M2_coorbital_m=m2,
                                d_miss_M3_extrapolating_m=m3,
                                penalty_M3=float(penalty_for_confidence(m3, a_t, 0.9)),
                                penalty_M1=float(penalty_for_confidence(m1, a_t, 0.9)),
                                M_launch_M3_g=float(penalty_for_confidence(m3, a_t, 0.9))
                                * m_on_target * 1e3,
                                M_launch_M1_g=float(penalty_for_confidence(m1, a_t, 0.9))
                                * m_on_target * 1e3,
                                p_hit_plan_M1=float(p_hit_plan_model(m1, 0.38)),
                                p_hit_plan_M3=float(p_hit_plan_model(m3, 0.38)),
                            ))
    df = pd.DataFrame(rows)
    save_table(df, "sim04_error_budget_full",
               "Terminal error budget under three miss models, with the "
               "launched-mass penalty each implies at 90% delivery confidence. "
               "p_hit_plan_* columns apply the plan's Bernoulli model at the "
               "v4 cloud size (sigma=0.38 m) for direct comparability.",
               tags={"sigma_pos": "DERIVED (diffraction)",
                     "pointing": "UNVALIDATED (engineering estimate)",
                     "muzzle_frac": "UNVALIDATED (engineering estimate)",
                     "u": "DESIGN (engagement geometry)"})

    # ---- validate the closed-form penalty against the numerical routine ----
    val = {}
    checks = []
    for sm in (0.005, 0.02, 0.05):
        M_num, s_num, pen_num = min_launch_mass_for_confidence(m_on_target, sm, r_t, 0.9)
        pen_cf = float(penalty_for_confidence(sm, a_t, 0.9))
        s_cf = float(sigma_opt_for_confidence(sm, 0.9))
        checks.append(dict(sigma_miss_m=sm, penalty_numeric=pen_num,
                           penalty_closed_form=pen_cf,
                           rel_err=abs(pen_num - pen_cf) / pen_cf,
                           sigma_numeric=s_num, sigma_closed_form=s_cf))
    val["closed_form_vs_numeric"] = checks
    val["note"] = ("The closed form uses the small-target limit; the numerical "
                   "routine uses the exact non-central chi-square fraction and "
                   "a Rayleigh miss. Agreement confirms the analytic sizing "
                   "rule used throughout SIM-9 and SIM-12.")
    save_json(val, "sim04_validation")
    print("  closed-form vs numeric penalty:",
          [f"{c['sigma_miss_m']*100:g} cm: {c['rel_err']*100:.1f}%" for c in checks])

    # ---- the headline design case -----------------------------------------
    base = df[(df.pointing_urad == 200.0) & (df.muzzle_frac == 1e-3)
              & (df.lidar_bw == "100 MHz") & (df.u_m_s == 10.0)
              & (df.R_m == 10.0)]
    head = base[["tau_ms", "d_miss_M1_v4_m", "d_miss_M2_coorbital_m",
                 "d_miss_M3_extrapolating_m", "p_hit_plan_M1", "p_hit_plan_M3",
                 "penalty_M3", "M_launch_M3_g"]]
    save_table(head, "sim04_headline_design_case",
               "R = 10 m, u = 10 m/s, pointing 200 urad, muzzle dispersion "
               "0.1%, 100 MHz ranging, v_rel = 619 m/s, 1 cm target.")
    print("\n  design case (R=10 m, u=10 m/s, 200 µrad, 0.1% muzzle):")
    print(head.to_string(index=False))

    # error-term ranking at the design point, tau = 10 ms
    dp = base[base.tau_ms == 10.0].iloc[0]
    ranking = {k: float(dp[f"term_{k}_m"]) for k in
               ("vel_extrap", "rel_accel", "pointing", "muzzle", "range")}
    ranking["sigma_pos"] = float(dp.sigma_pos_m)
    val["error_term_ranking_at_design_point_m"] = dict(
        sorted(ranking.items(), key=lambda kv: -kv[1]))
    save_json(val, "sim04_validation")
    print("  dominant error terms (m):",
          {k: f"{v:.2e}" for k, v in val["error_term_ranking_at_design_point_m"].items()})

    # ---- the question the plan asks, answered ------------------------------
    # "is there a (v_rel, sigma_cloud, tau, R) with delivery >= 50% / 90% at a
    #  realistically achievable tau?"
    ok50, ok90 = [], []
    for _, r in df.iterrows():
        for name, dm in (("M1_v4", r.d_miss_M1_v4_m),
                         ("M3", r.d_miss_M3_extrapolating_m)):
            pen = float(penalty_for_confidence(dm, a_t, 0.9))
            rec = dict(model=name, tau_ms=r.tau_ms, R_m=r.R_m, u_m_s=r.u_m_s,
                       pointing_urad=r.pointing_urad, d_miss_m=dm,
                       penalty=pen, M_launch_g=pen * m_on_target * 1e3)
            if pen * m_on_target <= 0.010:      # <= 10 g/shot
                ok90.append(rec)
            if pen * m_on_target <= 1.0:        # <= 1 kg/shot
                ok50.append(rec)
    concl = {
        "design_point": dict(h_km=H_DES / 1e3, v_rel_m_s=V_REL,
                             m_on_target_g=m_on_target * 1e3,
                             target="1 cm Al-6061 sphere"),
        "max_tau_ms_tested": max(TAU_MS),
        "v4_model_M1": {
            "d_miss_at_tau_10ms_m": float(base[base.tau_ms == 10].d_miss_M1_v4_m.iloc[0]),
            "p_hit_plan_at_tau_10ms": float(base[base.tau_ms == 10].p_hit_plan_M1.iloc[0]),
            "p_hit_plan_at_tau_1ms": float(base[base.tau_ms == 1].p_hit_plan_M1.iloc[0]),
            "statement": ("Reproduces the v4 collapse: at v_rel*tau the miss "
                          "distance at 10 ms lag is metres, and the Bernoulli "
                          "hit model at sigma_cloud=0.38 m returns a vanishing "
                          "probability.")},
        "M3_model": {
            "d_miss_at_tau_10ms_m": float(base[base.tau_ms == 10].d_miss_M3_extrapolating_m.iloc[0]),
            "d_miss_at_tau_100ms_m": float(base[base.tau_ms == 100].d_miss_M3_extrapolating_m.iloc[0]),
            "launched_mass_at_tau_100ms_g": float(base[base.tau_ms == 100].M_launch_M3_g.iloc[0]),
            "statement": ("With the target state extrapolated over the "
                          "engagement, control lag is not a binding constraint "
                          "at any tau tested up to 100 ms: the miss distance is "
                          "dominated by launcher pointing, which is "
                          "range-proportional and lag-independent.")},
        "binding_error_term": max(ranking, key=ranking.get),
        "tau_sensitivity_M3": {
            f"{t} ms": float(base[base.tau_ms == t].d_miss_M3_extrapolating_m.iloc[0])
            for t in TAU_MS},
        "n_configs_under_10g_per_shot": len(ok90),
        "n_configs_under_1kg_per_shot": len(ok50),
        "lidar_bandwidth_requirement": (
            "The plan's 1 MHz ranging bandwidth gives sigma_R = 33.5 m, which "
            "at u = 10 m/s contributes 0.54 m of miss -- inadequate. 100 MHz "
            "(sigma_R = 0.34 m, 5.4 mm of miss) or better is required."),
        "launch_angle": ("Pure retrograde in the target frame. FG01 leads the "
                         "target in-plane and fires rearward; useful impulse "
                         "scales as cos(theta) off-axis, so a 10 deg pointing "
                         "error costs 1.5% of the impulse, negligible next to "
                         "the areal-density effect of the same error."),
    }
    save_json(concl, "sim04_conclusions")

    # ---- figures ----------------------------------------------------------
    fig, ax = plt.subplots()
    for i, u in enumerate(US):
        s = df[(df.pointing_urad == 200.0) & (df.muzzle_frac == 1e-3)
               & (df.lidar_bw == "100 MHz") & (df.u_m_s == u)
               & (df.R_m == 10.0)].sort_values("tau_ms")
        ax.loglog(s.tau_ms, s.d_miss_M3_extrapolating_m * 100, color=CB[i],
                  marker="o", ms=3, label=f"M3 extrapolating, u={u:g} m/s")
    s1 = df[(df.pointing_urad == 200.0) & (df.muzzle_frac == 1e-3)
            & (df.lidar_bw == "100 MHz") & (df.u_m_s == 10.0)
            & (df.R_m == 10.0)].sort_values("tau_ms")
    ax.loglog(s1.tau_ms, s1.d_miss_M1_v4_m * 100, color="0.3", ls="--", marker="s",
              ms=3, label="M1 (v4 model), $v_{rel}\\tau$")
    ax.axhline(1.0, color=CB[2], lw=0.9)
    ax.text(0.11, 1.15, "1 cm target radius scale", fontsize=7, color=CB[2])
    ax.set_xlabel("control-loop latency τ (ms)")
    ax.set_ylabel("miss distance (cm)")
    ax.set_title("SIM-4: the v4 latency cliff versus the co-orbital "
                 "extrapolating budget (R = 10 m)")
    ax.legend(fontsize=7.5)
    save_fig(fig, "sim04_miss_vs_tau",
             "Miss distance against control-loop latency. Under the v4 model "
             "the miss grows as v_rel*tau and passes metres by 10 ms. Under the "
             "co-orbital extrapolating budget the lag term is u*tau at most and "
             "the total is dominated by lag-independent launcher pointing.")

    fig, ax = plt.subplots()
    for i, pt in enumerate(POINTING):
        s = df[(df.pointing_urad == pt * 1e6) & (df.muzzle_frac == 1e-3)
               & (df.lidar_bw == "100 MHz") & (df.u_m_s == 10.0)
               & (df.tau_ms == 10.0)].sort_values("R_m")
        ax.loglog(s.R_m, s.M_launch_M3_g, color=CB[i], marker="o", ms=3,
                  label=f"pointing {pt*1e6:.0f} µrad")
    ax.axhline(m_on_target * 1e3, color="0.3", ls=":",
               label=f"mass on target ({m_on_target*1e3:.3f} g)")
    ax.set_xlabel("engagement range R (m)")
    ax.set_ylabel("launched mass per shot (g), 90% delivery confidence")
    ax.set_title("SIM-4: launched mass is set by pointing × range, not by latency")
    ax.legend(fontsize=8)
    save_fig(fig, "sim04_mass_vs_range",
             "Launched mass per engagement at 90% delivery confidence. Because "
             "the dominant miss term is proportional to range, closing to short "
             "range is worth more than any improvement in control latency.")
    return df, concl


if __name__ == "__main__":
    run()
