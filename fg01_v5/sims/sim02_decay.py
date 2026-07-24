"""
SIM-2 -- Orbital decay and the accelerated-reentry success gate (NRLMSISE-00).

The primary physical success criterion of v5: does a Delta-h shift cut the
natural orbital lifetime below the regulatory bar, and at which altitudes?

Also settles the fate of *missed* projectile mass, which v4 assumed remains on
a near-circular orbit.  In the engagement geometry the energy gate forces
(retrograde intercept at several hundred m/s), that assumption is wrong: the
missed mass is left on a trajectory whose perigee is far below the atmosphere.
Both cases are computed and reported.
"""
import sys, pathlib, time
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from fg01.constants import (ALT_SWEEP_KM, RHO_AL, RHO_RE, RHO_W, CD_SPHERE,
                            F107_CASES, NRLMSISE_IMPL, ODMSP_YR, FCC_DISPOSAL_YR,
                            GRAIN_DIAM_M)
from fg01.orbital import (v_circ, dv_shift, dv_direct_reentry, r_of_h,
                          elements_from_apo_peri, area_to_mass_sphere,
                          mass_sphere, period, vis_viva)
from fg01.decay import lifetime, lifetime_years, propagate_cartesian, YEAR, H_REENTRY
from fg01.atmosphere import density
from fg01.io_utils import save_table, save_fig, save_json, banner, CB
from fg01.constants import MU, R_E

DEBRIS_D = [0.005, 0.01, 0.02, 0.05, 0.10]
CASES = ["solar_min", "nominal", "solar_max"]
INC_DEG = 85.0          # representative debris-belt inclination; sensitivity run below


def _lt(a, e, aom, case, t_max_yr=200.0):
    t, _ = lifetime(a, e, aom, inc_deg=INC_DEG, cd=CD_SPHERE, case=case,
                    t_max=t_max_yr * YEAR)
    return t / YEAR


def validation():
    """Starshine 1 cross-check + averaged-vs-Cartesian propagator agreement."""
    out = {"nrlmsise_implementation": NRLMSISE_IMPL}

    # --- Starshine 1 (NORAD 25769): 0.48 m diameter, 39.5 kg, deployed from
    # STS-96 on 1999-06-05 at ~387 km circular, i = 51.6 deg, reentered
    # 2000-02-18 -> 258 days observed.  Rising phase of solar cycle 23.
    d, m = 0.48, 39.5
    aom = (np.pi * (d / 2) ** 2) / m
    a0 = R_E + 387e3
    res = {}
    for case in CASES:
        t, _ = lifetime(a0, 0.0, aom, inc_deg=51.6, cd=CD_SPHERE, case=case,
                        t_max=20 * YEAR)
        res[case] = t / 86400.0
    obs = 258.0
    out["starshine1"] = {
        "A_over_m_m2_per_kg": aom, "observed_days": obs,
        "modelled_days": res,
        "rel_error_vs_observed": {k: (v - obs) / obs for k, v in res.items()},
        "note": ("Observed lifetime falls between the nominal and solar-max "
                 "model cases, consistent with the rising phase of cycle 23 "
                 "(F10.7 ~ 130-200 over 1999-2000). Fixed-F10.7 cases bracket "
                 "the observation; they are not expected to reproduce it "
                 "exactly.")}

    # --- averaged vs full 3-D Cartesian, same drag model, 30-day window
    a, e = elements_from_apo_peri(800e3, 600e3)
    aom_d = float(area_to_mass_sphere(0.01, RHO_AL))
    inc = np.deg2rad(INC_DEG)
    # start at apogee
    r0 = np.array([a * (1 + e), 0.0, 0.0])
    vmag = vis_viva(a * (1 + e), a)
    v0 = np.array([0.0, vmag * np.cos(inc), vmag * np.sin(inc)])
    # exactly an integer number of revolutions, J2 off, so that the comparison
    # isolates the orbit-averaging error (see propagate_cartesian docstring)
    n_rev = 400
    t_end = n_rev * float(period(a))
    cart = {}
    for dt_step in (20.0, 10.0):
        _, a_c, e_c, _, _ = propagate_cartesian(r0, v0, aom_d, case="nominal",
                                                t_end=t_end, dt=dt_step,
                                                use_j2=False)
        cart[dt_step] = (a_c, e_c)
    a_c, e_c = cart[10.0]
    from fg01.decay import _rates
    a_a, e_a, t = a, e, 0.0
    while t < t_end:
        step = min(3 * 3600.0, t_end - t)
        da, de = _rates(a_a, e_a, inc, aom_d, CD_SPHERE, "nominal")
        a_m, e_m = a_a + 0.5 * step * da, max(e_a + 0.5 * step * de, 0.0)
        da2, de2 = _rates(a_m, e_m, inc, aom_d, CD_SPHERE, "nominal")
        a_a += step * da2
        e_a = max(e_a + step * de2, 0.0)
        t += step
    out["averaged_vs_cartesian"] = {
        "window_revolutions": n_rev,
        "window_days": t_end / 86400.0,
        "delta_a_cartesian_m": float(a_c - a), "delta_a_averaged_m": float(a_a - a),
        "rel_diff": float(abs((a_c - a) - (a_a - a)) / abs(a_c - a)),
        "delta_e_cartesian": float(e_c - e), "delta_e_averaged": float(e_a - e),
        "cartesian_stepsize_convergence_rel": float(
            abs(cart[20.0][0] - cart[10.0][0]) / abs(cart[10.0][0] - a)),
        "note": ("J2 disabled in the Cartesian run and the window set to an "
                 "integer number of revolutions, so this isolates the "
                 "orbit-averaging error rather than J2 osculating "
                 "short-period terms.")}

    # --- fast-decay case propagated fully in Cartesian to reentry
    a2 = R_E + 250e3
    r0 = np.array([a2, 0.0, 0.0])
    vmag = np.sqrt(MU / a2)
    v0 = np.array([0.0, vmag * np.cos(inc), vmag * np.sin(inc)])
    t_c2, *_ , re2 = propagate_cartesian(r0, v0, aom_d, case="nominal",
                                         t_end=400 * 86400.0, dt=10.0,
                                         use_j2=False)
    t_a2 = _lt(a2, 0.0, aom_d, "nominal")
    out["fast_decay_250km_1cm"] = {
        "cartesian_days": t_c2 / 86400.0, "cartesian_reentered": bool(re2),
        "averaged_days": t_a2 * 365.25,
        "rel_diff": float(abs(t_c2 / 86400.0 - t_a2 * 365.25) / (t_c2 / 86400.0))}
    return out


def run():
    banner("SIM-2  Orbital decay and the accelerated-reentry gate")
    t0 = time.time()

    val = validation()
    save_json(val, "sim02_validation")
    print("  Starshine 1: observed 258 d, model",
          {k: round(v, 1) for k, v in val["starshine1"]["modelled_days"].items()})
    v_ac = val["averaged_vs_cartesian"]
    print(f"  averaged vs Cartesian ({v_ac['window_revolutions']} rev, "
          f"{v_ac['window_days']:.0f} d, Δa): "
          f"{v_ac['delta_a_cartesian_m']:.1f} m vs "
          f"{v_ac['delta_a_averaged_m']:.1f} m "
          f"({v_ac['rel_diff']*100:.2f}%); Cartesian step convergence "
          f"{v_ac['cartesian_stepsize_convergence_rel']*100:.3f}%")
    print("  250 km 1 cm full-Cartesian to reentry:",
          f"{val['fast_decay_250km_1cm']['cartesian_days']:.1f} d vs averaged "
          f"{val['fast_decay_250km_1cm']['averaged_days']:.1f} d "
          f"({val['fast_decay_250km_1cm']['rel_diff']*100:.1f}%)")

    # ---------------- debris lifetimes, pre- and post-kick -------------------
    rows = []
    for d in DEBRIS_D:
        aom = float(area_to_mass_sphere(d, RHO_AL))
        for h in ALT_SWEEP_KM:
            hm = h * 1e3
            for case in CASES:
                t_pre = _lt(R_E + hm, 0.0, aom, case)
                for dh in (150e3, 200e3, 250e3):
                    a, e = elements_from_apo_peri(hm, hm - dh)
                    t_post = _lt(a, e, aom, case)
                    rows.append(dict(
                        d_debris_cm=d * 100, h_km=h, dh_km=dh / 1e3,
                        solar_case=case, f107=F107_CASES[case],
                        A_over_m=aom,
                        T_pre_yr=t_pre, T_post_yr=t_post,
                        acceleration_factor=(t_pre / t_post
                                             if np.isfinite(t_pre) and t_post > 0
                                             else np.inf),
                        pass_25yr=bool(t_post < ODMSP_YR),
                        pass_5yr=bool(t_post < FCC_DISPOSAL_YR)))
    df = pd.DataFrame(rows)
    save_table(df, "sim02_decay_full",
               "Pre- and post-kick orbital lifetime for solid Al-6061 spheres, "
               "orbit-averaged Gauss propagation on real NRLMSISE-00 density.",
               tags={"rho": "SOURCED(NRLMSISE-00)", "C_D": "SOURCED",
                     "inclination": f"{INC_DEG} deg (representative)"})

    # ---------------- altitude ceilings -------------------------------------
    ceil = []
    for d in DEBRIS_D:
        for dh in (150.0, 200.0, 250.0):
            for case in CASES:
                s = df[(df.d_debris_cm == d * 100) & (df.dh_km == dh)
                       & (df.solar_case == case)].sort_values("h_km")
                def ceiling(col):
                    ok = s[s[col]]
                    return float(ok.h_km.max()) if len(ok) else np.nan
                ceil.append(dict(d_debris_cm=d * 100, dh_km=dh, solar_case=case,
                                 h_ceiling_25yr_km=ceiling("pass_25yr"),
                                 h_ceiling_5yr_km=ceiling("pass_5yr")))
    dfc = pd.DataFrame(ceil)
    save_table(dfc, "sim02_altitude_ceilings",
               "Highest starting altitude at which a single Delta-h kick puts "
               "the object inside the 25-yr and 5-yr disposal rules.")
    print("\n  Altitude ceilings (nominal solar, dh=200 km):")
    print(dfc[(dfc.dh_km == 200) & (dfc.solar_case == "nominal")].to_string(index=False))

    # ---------------- staged (multi-kick) deorbit ---------------------------
    # A single kick is not the only option: repeating the engagement lowers the
    # object 200 km at a time.  Each kick is individually sub-catastrophic, so
    # the altitude ceiling becomes an economic limit, not a physical one.
    stage_rows = []
    for d in DEBRIS_D:
        aom = float(area_to_mass_sphere(d, RHO_AL))
        for case in CASES:
            for h in ALT_SWEEP_KM:
                for rule, bar in (("25yr", ODMSP_YR), ("5yr", FCC_DISPOSAL_YR)):
                    n = 0
                    h_cur = float(h)
                    t_final = np.inf
                    while n <= 12:
                        if n == 0:
                            a, e = R_E + h_cur * 1e3, 0.0
                        else:
                            a, e = elements_from_apo_peri(h_cur * 1e3,
                                                          (h_cur - 200.0) * 1e3)
                        t_f = _lt(a, e, aom, case)
                        if t_f < bar:
                            t_final = t_f
                            break
                        n += 1
                        h_cur -= 200.0
                        if h_cur < 200.0:
                            break
                    stage_rows.append(dict(d_debris_cm=d * 100, h_km=h,
                                           solar_case=case, rule=rule,
                                           n_kicks=n if np.isfinite(t_final) else np.nan,
                                           T_final_yr=t_final))
    dfs = pd.DataFrame(stage_rows)
    save_table(dfs, "sim02_staged_kicks",
               "Number of successive 200 km kicks needed to bring an object "
               "inside the 25-yr / 5-yr rule (each kick sub-catastrophic).")
    print("\n  Staged kicks to meet 25-yr rule (nominal solar):")
    piv = dfs[(dfs.rule == "25yr") & (dfs.solar_case == "nominal")].pivot(
        index="h_km", columns="d_debris_cm", values="n_kicks")
    print(piv.to_string())

    # ---------------- fate of missed projectile mass ------------------------
    miss_rows = []
    for h in ALT_SWEEP_KM:
        hm = h * 1e3
        r = float(r_of_h(hm))
        vc = float(v_circ(hm))
        dvd = float(dv_direct_reentry(hm))
        for v_rel in (200.0, 400.0, 619.0, 800.0, 1000.0):
            v_new = vc - v_rel                      # retrograde launch
            energy = 0.5 * v_new ** 2 - MU / r
            a_new = -MU / (2 * energy) if energy < 0 else np.inf
            if np.isfinite(a_new) and a_new > 0:
                e_new = 1.0 - (min(r, 2 * a_new - r)) / a_new if a_new > r / 2 else 1.0
                rp = 2 * a_new - r if 2 * a_new - r < r else r
                hp = rp - R_E
                immediate = hp < H_REENTRY
                t_to_perigee = period(a_new) / 2.0 / 60.0 if immediate else np.nan
            else:
                hp, immediate, t_to_perigee, a_new = -np.inf, True, np.nan, np.inf
            for dgr in GRAIN_DIAM_M:
                aom_re = float(area_to_mass_sphere(dgr, RHO_RE))
                aom_w = float(area_to_mass_sphere(dgr, RHO_W))
                if immediate:
                    t_re = t_w = 0.0
                else:
                    ee = (r - rp) / (r + rp)
                    t_re = _lt((r + rp) / 2, ee, aom_re, "nominal")
                    t_w = _lt((r + rp) / 2, ee, aom_w, "nominal")
                # counterfactual: grain left on the original circular orbit
                t_cf_re = _lt(r, 0.0, aom_re, "nominal")
                miss_rows.append(dict(
                    h_km=h, v_rel_m_s=v_rel, dv_direct_m_s=dvd,
                    d_grain_mm=dgr * 1e3,
                    perigee_alt_km=hp / 1e3 if np.isfinite(hp) else -np.inf,
                    immediate_reentry=bool(immediate),
                    minutes_to_reentry=t_to_perigee,
                    T_decay_Re_yr=t_re, T_decay_W_yr=t_w,
                    T_decay_if_left_circular_yr=t_cf_re,
                    fcc_5yr_pass=bool((t_re if not immediate else 0.0) < FCC_DISPOSAL_YR)))
    dfm = pd.DataFrame(miss_rows)
    save_table(dfm, "sim02_missed_mass_fate",
               "Fate of projectile mass that misses. A retrograde intercept at "
               "v_rel > dv_direct leaves the missed mass on a sub-orbital "
               "trajectory that reenters on its first perigee pass; the "
               "counterfactual column is the v4 assumption of a grain left on "
               "the original circular orbit.")
    imm = dfm.groupby(["h_km", "v_rel_m_s"]).immediate_reentry.first().reset_index()
    print("\n  Missed-mass immediate reentry (True = reenters within one orbit):")
    print(imm.pivot(index="h_km", columns="v_rel_m_s",
                    values="immediate_reentry").to_string())

    # ---------------- conclusions -------------------------------------------
    n = df[(df.solar_case == "nominal") & (df.dh_km == 200)]
    concl = {
        "propagator": "orbit-averaged Gauss variational equations, 96-node "
                      "Gauss-Legendre quadrature in eccentric anomaly, vector "
                      "atmospheric co-rotation, real NRLMSISE-00 density",
        "nrlmsise_implementation": NRLMSISE_IMPL,
        "inclination_deg": INC_DEG,
        "ceilings_nominal_dh200": dfc[(dfc.dh_km == 200) &
                                      (dfc.solar_case == "nominal")]
            .set_index("d_debris_cm")[["h_ceiling_25yr_km", "h_ceiling_5yr_km"]]
            .to_dict(orient="index"),
        "acceleration_factor_1cm_nominal": {
            str(int(r.h_km)): (None if not np.isfinite(r.acceleration_factor)
                               else round(float(r.acceleration_factor), 2))
            for _, r in n[n.d_debris_cm == 1.0].iterrows()},
        "missed_mass_immediate_reentry_at_619": bool(
            dfm[(dfm.v_rel_m_s == 619.0)].immediate_reentry.all()),
        "min_v_rel_for_self_disposal_m_s": {
            str(int(h)): round(float(dv_direct_reentry(h * 1e3)), 1)
            for h in ALT_SWEEP_KM},
        "runtime_s": round(time.time() - t0, 1),
    }
    save_json(concl, "sim02_conclusions")

    # ---------------- figures ------------------------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.0), sharey=True)
    for ax, case in zip(axes, CASES):
        for i, d in enumerate([0.5, 1.0, 2.0, 5.0, 10.0]):
            s = df[(df.d_debris_cm == d) & (df.dh_km == 200)
                   & (df.solar_case == case)].sort_values("h_km")
            ax.plot(s.h_km, np.clip(s.T_post_yr, 1e-3, 500), color=CB[i],
                    label=f"{d:g} cm, post-kick")
            ax.plot(s.h_km, np.clip(s.T_pre_yr, 1e-3, 500), color=CB[i],
                    ls=":", lw=1.0)
        ax.axhline(25, color="k", lw=1.2)
        ax.axhline(5, color="k", lw=1.2, ls="--")
        ax.set_yscale("log")
        ax.set_title(f"{case} (F10.7={F107_CASES[case]:.0f})")
        ax.set_xlabel("pre-kick altitude (km)")
    axes[0].set_ylabel("orbital lifetime (yr)")
    axes[0].text(405, 27, "25 yr", fontsize=7.5)
    axes[0].text(405, 5.4, "5 yr", fontsize=7.5)
    axes[0].legend(fontsize=7, ncol=1)
    fig.suptitle("SIM-2: lifetime after a 200 km kick (solid) vs. undisturbed "
                 "(dotted); values at 500 yr are capped", fontsize=10)
    save_fig(fig, "sim02_lifetime_vs_altitude",
             "Post-kick orbital lifetime against pre-kick altitude for five "
             "debris sizes and three solar activity cases. Crossings of the "
             "25-yr and 5-yr lines are the altitude ceilings.")

    fig, ax = plt.subplots()
    for i, d in enumerate([0.5, 1.0, 2.0, 5.0, 10.0]):
        s = dfs[(dfs.rule == "25yr") & (dfs.solar_case == "nominal")
                & (dfs.d_debris_cm == d)].sort_values("h_km")
        ax.step(s.h_km, s.n_kicks, where="mid", color=CB[i], label=f"{d:g} cm")
    ax.set_xlabel("pre-kick altitude (km)")
    ax.set_ylabel("successive 200 km kicks needed")
    ax.set_title("SIM-2: staged deorbit — kicks required to meet the 25-yr rule "
                 "(nominal solar)")
    ax.legend(fontsize=8, ncol=3)
    save_fig(fig, "sim02_staged_kicks",
             "Number of successive sub-catastrophic 200 km kicks needed to "
             "bring an object inside the 25-yr rule. The single-kick altitude "
             "ceiling is where this first exceeds one.")
    print(f"\n  SIM-2 runtime {time.time()-t0:.0f} s")
    return df, dfc, dfs, dfm, concl


if __name__ == "__main__":
    run()
