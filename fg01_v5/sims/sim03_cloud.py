"""
SIM-3 -- Cloud dispersion, evolution, and the cross-section SIM-4 consumes.

Three things are established here.

(1) Clohessy-Wiltshire relative motion is evaluated exactly (AM-14: no assumed
    cone angle).  At engagement time-of-flight (R/v_rel < 1 s for R <= 100 m)
    the CW solution is numerically indistinguishable from straight-line
    ballistic spreading: the orbital-curvature correction is O((n t)^2) ~ 1e-7.
    So the cloud cross-section is NOT an orbital-mechanics outcome at these
    scales -- it is set by the launcher's transverse velocity dispersion, i.e.
    it is a *design variable*.  CW matters only for release-to-intercept times
    of minutes or more, which this engagement does not use.

(2) The mass consequence.  A target immersed in the cloud intercepts only the
    mass incident on its own cross-section, so widening the cloud to tolerate a
    miss costs launched mass as sigma^2.  The Pareto curve between tolerance
    and mass is computed here and consumed by SIM-4.

(3) Mode and grain size.  Grain mass is bounded above by the requirement that
    enough grains strike the target for the delivered momentum to be
    statistically stable (Poisson), and bounded below by manufacturability and
    by reentry-survival considerations (SIM-8).
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from fg01.constants import RHO_RE, RHO_W, RHO_AL, GAMMA, MU, R_E
from fg01.orbital import mean_motion, mass_sphere, xsec_sphere, dv_shift
from fg01.delivery import (delivered_fraction, delivered_fraction_smalltarget,
                           optimal_sigma, launch_mass_deterministic,
                           grain_statistics)
from fg01.io_utils import save_table, save_fig, save_json, banner, CB

RANGES = [1.0, 5.0, 10.0, 50.0, 100.0]
V_RELS = [400.0, 619.0, 1000.0]
DVT = [1e-3, 1e-2, 0.1, 0.3, 1.0, 3.0]        # transverse velocity dispersion, m/s


def cw_spread(dv, n, t):
    """Exact CW displacement magnitudes for a unit-direction dv kick."""
    dz = (dv / n) * np.sin(n * t)
    dx = (2 * dv / n) * (1 - np.cos(n * t))
    dy = (3 * dv / n) * (n * t - np.sin(n * t))
    return dx, dy, dz


def run():
    banner("SIM-3  Cloud dispersion and the exported cross-section")

    # ---------- (1) CW vs ballistic at engagement timescales ---------------
    rows = []
    for h in (400e3, 600e3, 800e3, 1000e3):
        n = float(mean_motion(R_E + h))
        for v_rel in V_RELS:
            for R in RANGES:
                tof = R / v_rel
                for dv in DVT:
                    dx, dy, dz = cw_spread(dv, n, tof)
                    cw_mag = np.hypot(dx, dz)          # in-plane transverse
                    ballistic = dv * tof
                    rows.append(dict(
                        h_km=h / 1e3, v_rel_m_s=v_rel, R_m=R, tof_s=tof,
                        dv_T_m_s=dv, n_rad_s=n, n_times_tof=n * tof,
                        sigma_cw_m=cw_mag, sigma_ballistic_m=ballistic,
                        cw_over_ballistic=cw_mag / ballistic,
                        along_track_cw_m=dy))
    dfc = pd.DataFrame(rows)
    save_table(dfc, "sim03_cw_vs_ballistic",
               "Exact Clohessy-Wiltshire cloud spreading versus straight-line "
               "ballistic spreading at engagement times-of-flight.",
               tags={"CW": "DERIVED", "dv_T": "DESIGN (launcher dispersion)"})
    dev = float(np.max(np.abs(dfc.cw_over_ballistic - 1.0)))
    print(f"  max |CW/ballistic - 1| over all engagement cases: {dev:.3e}")
    print(f"  max n*TOF: {dfc.n_times_tof.max():.3e}  (curvature is O((n t)^2))")

    # long-timescale CW, where it does matter (context)
    long_rows = []
    n800 = float(mean_motion(R_E + 800e3))
    for t in (1.0, 10.0, 60.0, 300.0, 1800.0, 5400.0):
        dx, dy, dz = cw_spread(0.01, n800, t)
        long_rows.append(dict(t_s=t, dx_m=dx, dy_m=dy, dz_m=dz,
                              ballistic_m=0.01 * t,
                              along_track_over_ballistic=dy / (0.01 * t)))
    save_table(pd.DataFrame(long_rows), "sim03_cw_long_timescale",
               "CW evolution of a 10 mm/s dispersion at 800 km over release-to-"
               "intercept times from 1 s to one orbit (context: shows where CW "
               "departs from ballistic, and that the engagement does not "
               "operate there).")

    # ---------- (2) the mandated sigma_cloud export ------------------------
    exp_rows = []
    for v_rel in V_RELS:
        for R in RANGES:
            tof = R / v_rel
            for dv in DVT:
                sig = dv * tof
                for d_t in (0.01, 0.05, 0.10):
                    r_t = d_t / 2
                    f = float(delivered_fraction(sig, 0.0, r_t))
                    f_small = float(delivered_fraction_smalltarget(
                        sig, 0.0, xsec_sphere(d_t)))
                    exp_rows.append(dict(
                        v_rel_m_s=v_rel, R_m=R, tof_s=tof, dv_T_m_s=dv,
                        sigma_cloud_m=sig, d_target_cm=d_t * 100,
                        delivered_fraction_centred=f,
                        small_target_approx=f_small,
                        mass_penalty_centred=1.0 / f if f > 0 else np.inf))
    dfe = pd.DataFrame(exp_rows)
    save_table(dfe, "sim03_sigma_cloud_export",
               "MANDATED EXPORT to SIM-4: cloud cross-section sigma_cloud at "
               "the exact times-of-flight used for R = 1, 5, 10, 50, 100 m, "
               "with the delivered-mass fraction for a perfectly centred shot.",
               tags={"sigma_cloud": "DERIVED from dv_T (DESIGN)"})

    # ---------- (3) the Pareto the plan asks for ---------------------------
    par_rows = []
    d_t = 0.01
    r_t = d_t / 2
    for d_miss in np.geomspace(0.002, 2.0, 40):
        s_opt = optimal_sigma(d_miss, r_t)
        M_over_m, s_used, penalty = launch_mass_deterministic(1.0, d_miss, r_t)
        par_rows.append(dict(d_miss_m=d_miss, sigma_opt_m=s_opt,
                             sigma_opt_over_d_miss=s_opt / d_miss,
                             mass_penalty=penalty,
                             closed_form_penalty=np.pi * np.e * d_miss ** 2
                             / xsec_sphere(d_t)))
    dfp = pd.DataFrame(par_rows)
    save_table(dfp, "sim03_pareto_mass_penalty",
               "Launched-mass penalty versus tolerated miss distance for a "
               "1 cm target, with the optimal cloud sigma. The closed-form "
               "column is the small-target limit  pi*e*d_miss^2/A_t.")
    print(f"  optimal sigma/d_miss (asymptote 1/sqrt2 = 0.707): "
          f"{dfp.sigma_opt_over_d_miss.iloc[-1]:.4f}")
    print("  mass penalty at d_miss = 1, 3, 10, 30 cm:",
          [round(float(np.interp(x, dfp.d_miss_m, dfp.mass_penalty)), 1)
           for x in (0.01, 0.03, 0.10, 0.30)])

    # ---------- grain-size / mode recommendation ---------------------------
    # design mass on target at 800 km, beta = 1.5, v_rel = 619 m/s
    m_t = float(mass_sphere(0.01, RHO_AL))
    m_p_design = GAMMA * float(dv_shift(800e3, 200e3)) * m_t / (1.5 * 619.0)
    gr_rows = []
    for d_miss in (0.01, 0.03, 0.10, 0.30):
        s = optimal_sigma(d_miss, r_t)
        M, _, penalty = launch_mass_deterministic(m_p_design, d_miss, r_t)
        for m_grain_mg in (0.001, 0.01, 0.1, 1.0, 10.0, 100.0):
            m_grain = m_grain_mg * 1e-6
            n_tot, n_hit, cv = grain_statistics(M, m_grain, s, d_miss, r_t)
            d_re = 2 * (3 * m_grain / (4 * np.pi * RHO_RE)) ** (1 / 3)
            gr_rows.append(dict(
                d_miss_m=d_miss, sigma_cloud_m=s, M_launch_g=M * 1e3,
                mass_penalty=penalty, m_grain_mg=m_grain_mg,
                d_grain_Re_mm=d_re * 1e3,
                n_grains_total=n_tot, n_grains_on_target=n_hit,
                poisson_cv=cv,
                stable_delivery=bool(n_hit >= 100),
                Es_per_impact_J_kg=0.5 * m_grain * 619.0 ** 2 / m_t))
    dfg = pd.DataFrame(gr_rows)
    save_table(dfg, "sim03_grain_mode",
               "Grain-size selection: number of grains striking a 1 cm target, "
               "Poisson variability of the delivered momentum, and per-impact "
               "specific energy. Design mass on target 0.161 g "
               "(800 km, beta=1.5, v_rel=619 m/s).",
               tags={"per_impact_Es": "DERIVED", "swarm equivalence": "AM-15/16 UNVALIDATED"})
    print("\n  grain mode at d_miss = 10 cm:")
    print(dfg[dfg.d_miss_m == 0.10][
        ["m_grain_mg", "d_grain_Re_mm", "M_launch_g", "n_grains_total",
         "n_grains_on_target", "poisson_cv", "Es_per_impact_J_kg"]].to_string(index=False))

    ok = dfg[(dfg.d_miss_m == 0.10) & dfg.stable_delivery]
    concl = {
        "cw_equals_ballistic_at_engagement": {
            "max_abs_dev": dev,
            "statement": ("Over engagement times-of-flight (<=0.25 s for "
                          "R<=100 m) the exact CW solution differs from "
                          "straight-line ballistic spreading by <1e-6 relative. "
                          "The cloud cross-section is therefore set by launcher "
                          "transverse velocity dispersion, a design variable, "
                          "not by orbital dynamics.")},
        "sigma_opt_rule": "sigma* = d_miss / sqrt(2) (small-target limit)",
        "mass_penalty_closed_form": "M_launch/m_on_target = pi*e*d_miss^2/A_target",
        "mass_penalty_1cm_target": {
            f"{x*100:g} cm miss": round(float(np.interp(x, dfp.d_miss_m, dfp.mass_penalty)), 1)
            for x in (0.01, 0.03, 0.10, 0.30)},
        "recommended_grain_mass_mg_max": float(ok.m_grain_mg.max()) if len(ok) else None,
        "recommended_grain_diameter_Re_mm_max": float(ok.d_grain_Re_mm.max()) if len(ok) else None,
        "mode_recommendation": (
            "Discrete grains, not dust and not a single slug. Grain mass must "
            "be small enough that >=100 grains strike the target (Poisson CV "
            "<=10%); at a 10 cm miss tolerance that caps grain mass near 1 mg "
            "(~0.45 mm rhenium). Dust (0.01 mg) also satisfies the statistics "
            "but is harder to launch coherently and is more strongly perturbed "
            "by residual gas and charging; 0.1-1 mg is the recommended band."),
    }
    save_json(concl, "sim03_conclusions")

    # ---------- figures ----------------------------------------------------
    fig, ax = plt.subplots()
    ax.loglog(dfp.d_miss_m * 100, dfp.mass_penalty, color=CB[0], lw=2,
              label="exact (non-central χ², 1 cm target)")
    ax.loglog(dfp.d_miss_m * 100, dfp.closed_form_penalty, color=CB[1], ls="--",
              label=r"closed form  $\pi e\, d_{miss}^2/A_t$")
    for d_t2, c in ((0.05, CB[2]), (0.10, CB[3])):
        pen = [launch_mass_deterministic(1.0, dm, d_t2 / 2)[2] for dm in dfp.d_miss_m]
        ax.loglog(dfp.d_miss_m * 100, pen, color=c, label=f"{d_t2*100:g} cm target")
    ax.axhline(1, color="0.4", lw=0.8)
    ax.set_xlabel("tolerated miss distance (cm)")
    ax.set_ylabel("launched mass / mass on target")
    ax.set_title("SIM-3: the mass price of guidance tolerance")
    ax.legend(fontsize=8)
    save_fig(fig, "sim03_mass_penalty",
             "Launched-mass penalty as a function of the miss distance the shot "
             "pattern must tolerate. The penalty is quadratic in miss distance "
             "and inversely proportional to target cross-section, which is why "
             "cost per object rises steeply for the smallest debris.")

    fig, ax = plt.subplots()
    for i, v in enumerate(V_RELS):
        for j, dv in enumerate((0.01, 0.1, 1.0)):
            s = dfe[(dfe.v_rel_m_s == v) & (dfe.dv_T_m_s == dv)
                    & (dfe.d_target_cm == 1.0)]
            ax.loglog(s.R_m, s.sigma_cloud_m * 100, color=CB[i],
                      ls=["-", "--", ":"][j],
                      label=f"$v_{{rel}}$={v:.0f}, $\\delta v_T$={dv:g} m/s")
    ax.set_xlabel("engagement range R (m)")
    ax.set_ylabel("cloud σ at intercept (cm)")
    ax.set_title("SIM-3: shot-pattern size, σ = δ$v_T$·R/$v_{rel}$")
    ax.legend(fontsize=7, ncol=2)
    save_fig(fig, "sim03_sigma_cloud",
             "Cloud cross-section at intercept. Because time-of-flight is "
             "R/v_rel, the pattern size is set by the ratio of transverse "
             "dispersion to intercept velocity times range — a launcher design "
             "choice, tunable to match the achievable miss distance.")
    return dfc, dfe, dfp, dfg, concl


if __name__ == "__main__":
    run()
