"""
INT-0 -- Engagement geometry and the co-orbital premise (the validation gate).

The v6 plan calls premise #3 "the most attackable claim in the whole project"
and asks for it to be tested independently, including the worst-case
near-transverse approach, with the honest boundary stated.

This module answers it in three parts:

  A. KINEMATICS (decides the architecture on its own).  Only the retrograde
     component of the impact relative velocity lowers the target's orbit.  For a
     crossing engagement that component is sin(theta/2) of the total, while the
     specific energy scales with the total squared.  Useful-delta-v per unit
     specific energy is 2*beta*cos(psi)/|w|, which evaluates to beta/v_orbital
     for *any* crossing angle and 2*beta/|w| for a co-orbital shot.  Co-orbital
     wins by 2*v_orbital/|w| ~ 24x at the design point.  This is kinematics, not
     an assumption about hardware.

  B. GUIDANCE.  Sweep approach geometry and control latency, under both a
     no-lead aim law (the v4/v5-plan form) and a lead-computing one, and report
     the cone within which the latency-independent miss survives.

  C. TIMELINE.  A crossing target is inside tracking range for a bounded time.
     The full detect-track-solve-fire sequence must fit inside it.  This, not
     the control loop, is what actually forbids fast-crossing engagements.

Full 3-D relative propagation is used as the check on the closed forms.
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from fg01.constants import MU, R_E, GAMMA, ES_C, RHO_AL, OPTICAL_LAMBDA_M, OPTICAL_APERTURE_M
from fg01.orbital import v_circ, mass_sphere, xsec_sphere, dv_shift, state_from_elements
from fg01.relmotion import (plane_angle, v_rel_crossing, retrograde_efficiency,
                            useful_dv_per_specific_energy, engagement_window,
                            miss_no_lead, miss_with_lead, sigma_v_from_track,
                            coorbital_drift_rate, dv_plane_change)
from fg01.io_utils import save_table, save_fig, save_json, banner, CB

H = 900e3
V_ORB = float(v_circ(H))
TAUS = [1e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1]
RANGES = [1.0, 10.0, 100.0, 1000.0]


def part_a_kinematics():
    """Retrograde efficiency and energy cost vs. approach geometry."""
    rows = []
    for theta_deg in [0.0, 0.01, 0.05, 0.1, 0.5, 1.0, 2.0, 5.0, 8.5, 15.0, 30.0, 45.0, 90.0]:
        th = np.radians(theta_deg)
        if theta_deg == 0.0:
            # co-orbital: |w| supplied entirely by the launcher, pure retrograde
            for w in (400.0, 619.0, 1000.0):
                rows.append(dict(
                    mode="co-orbital (launcher-supplied)", theta_deg=0.0,
                    w_m_s=w, cos_psi=1.0, useful_frac=1.0,
                    Es_kJ_kg_at_design=GAMMA * float(dv_shift(H, 200e3)) * w / (2 * 1.15) / 1e3,
                    useful_dv_per_Es=float(useful_dv_per_specific_energy(1.15, w, 1.0)),
                    v_transverse_m_s=np.nan))
        else:
            w = float(v_rel_crossing(V_ORB, th))
            cpsi = float(retrograde_efficiency(th))
            rows.append(dict(
                mode="crossing (geometry-supplied)", theta_deg=theta_deg,
                w_m_s=w, cos_psi=cpsi, useful_frac=cpsi,
                Es_kJ_kg_at_design=np.nan,
                useful_dv_per_Es=float(useful_dv_per_specific_energy(1.15, w, cpsi)),
                v_transverse_m_s=w))
    df = pd.DataFrame(rows)

    # what a crossing engagement costs in disruption energy for a fixed useful dv
    dv_needed = float(dv_shift(H, 200e3))
    m_t = float(mass_sphere(0.01, RHO_AL))
    cost = []
    for theta_deg in [0.1, 0.5, 1.0, 2.0, 5.0, 8.5, 15.0, 45.0, 90.0]:
        th = np.radians(theta_deg)
        w = float(v_rel_crossing(V_ORB, th))
        cpsi = float(retrograde_efficiency(th))
        beta = 1.15
        # m_p to deliver dv_needed retrograde
        m_p = dv_needed * m_t / (beta * w * cpsi)
        es = 0.5 * m_p * w ** 2 / m_t
        cost.append(dict(theta_deg=theta_deg, w_m_s=w, cos_psi=cpsi,
                         m_p_required_g=m_p * 1e3, Es_kJ_kg=es / 1e3,
                         catastrophic=bool(es >= ES_C),
                         times_threshold=es / ES_C))
    dfc = pd.DataFrame(cost)
    m_p_co = dv_needed * m_t / (1.15 * 619.0)
    es_co = 0.5 * m_p_co * 619.0 ** 2 / m_t
    dfc.loc[len(dfc)] = dict(theta_deg=0.0, w_m_s=619.0, cos_psi=1.0,
                             m_p_required_g=m_p_co * 1e3, Es_kJ_kg=es_co / 1e3,
                             catastrophic=bool(es_co >= ES_C),
                             times_threshold=es_co / ES_C)
    return df, dfc.sort_values("theta_deg")


def part_a2_combined():
    """
    The physically correct intermediate case: the platform is *nearly* co-orbital
    but its plane differs from the target's by theta, so it already drifts at
    v_geom = 2 v sin(theta/2) before firing.  The launcher then adds v_L
    retrograde.  The impact relative velocity is the vector sum:

        w = v_geom_vec + v_L * (-u_t)

    Only the launcher part is reliably retrograde; the geometric part is
    essentially cross-track.  So

        |w|^2 = v_geom^2 + v_L^2      (the two are close to orthogonal)
        retrograde component = v_L

    and the specific energy needed for a fixed retrograde kick grows as
    (1 + (v_geom/v_L)^2).  This gives ECO-1 a usable firing cone instead of the
    two idealised extremes.
    """
    m_t = float(mass_sphere(0.01, RHO_AL))
    dv_needed = float(dv_shift(H, 200e3))
    beta = 1.15
    rows = []
    for v_L in (400.0, 619.0, 1000.0):
        for theta_deg in [0.0, 0.05, 0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0]:
            v_geom = float(v_rel_crossing(V_ORB, np.radians(theta_deg)))
            w = np.hypot(v_geom, v_L)
            m_p = dv_needed * m_t / (beta * v_L)         # retrograde part only
            es = 0.5 * m_p * w ** 2 / m_t
            rows.append(dict(
                v_launch_m_s=v_L, theta_deg=theta_deg, v_geom_m_s=v_geom,
                v_geom_over_v_launch=v_geom / v_L,
                w_total_m_s=w, m_p_g=m_p * 1e3, Es_kJ_kg=es / 1e3,
                energy_penalty=1 + (v_geom / v_L) ** 2,
                times_threshold=es / ES_C,
                subcatastrophic=bool(es < ES_C),
                transverse_rate_m_s=v_geom))
    return pd.DataFrame(rows)


def part_b_guidance():
    """Miss distance vs. approach geometry, latency and aim law."""
    rows = []
    for theta_deg in [0.0, 0.001, 0.01, 0.1, 1.0, 8.5, 45.0]:
        th = np.radians(theta_deg)
        if theta_deg == 0.0:
            v_T_list = [(dh, float(coorbital_drift_rate(dh, H)))
                        for dh in (10.0, 100.0, 1000.0, 5000.0)]
        else:
            v_T_list = [(np.nan, float(v_rel_crossing(V_ORB, th)))]
        for dh, v_T in v_T_list:
            for R in RANGES:
                w = 619.0
                tof = R / w
                sp = (OPTICAL_LAMBDA_M / OPTICAL_APERTURE_M) * R / 10.0   # SNR=10
                # tracking arc is limited by how long the target is in view
                t_view = float(engagement_window(max(v_T, 1e-3), 10 * R))
                t_obs = min(10.0, t_view)
                sv = float(sigma_v_from_track(sp, max(t_obs, 1e-3), 1000.0))
                for tau in TAUS:
                    rows.append(dict(
                        theta_deg=theta_deg, dh_offset_m=dh, v_transverse_m_s=v_T,
                        R_m=R, tof_s=tof, tau_s=tau,
                        t_view_s=t_view, t_obs_s=t_obs,
                        sigma_pos_m=sp, sigma_v_m_s=sv,
                        miss_no_lead_m=float(miss_no_lead(v_T, tau, tof)),
                        miss_with_lead_m=float(miss_with_lead(sv, tau, tof)),
                        n_samples_in_view=t_view * 1000.0))
    return pd.DataFrame(rows)


def part_c_timeline():
    """Engagement-timeline feasibility for crossing targets."""
    rows = []
    for theta_deg in [0.0, 0.001, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 8.5, 45.0]:
        if theta_deg == 0.0:
            v_T = float(coorbital_drift_rate(1000.0, H))     # 1 km offset
        else:
            v_T = float(v_rel_crossing(V_ORB, np.radians(theta_deg)))
        for r_track in (100.0, 1000.0, 10000.0):
            t_win = float(engagement_window(v_T, r_track))
            rows.append(dict(
                theta_deg=theta_deg, v_transverse_m_s=v_T, r_track_m=r_track,
                window_s=t_win,
                # a full autonomous detect->track->solve->fire sequence
                feasible_1s_sequence=bool(t_win > 1.0),
                feasible_100ms_sequence=bool(t_win > 0.1),
                shots_possible_at_1Hz=t_win * 1.0,
                aim_slew_rate_deg_s=np.degrees(v_T / max(r_track / 10, 1.0))))
    return pd.DataFrame(rows)


def validation_3d():
    """
    Check the closed forms against explicit 3-D states: build two circular
    orbits with a given plane separation, place both at the common node, and
    measure |w| and the retrograde component directly.
    """
    out = []
    a = R_E + H
    for theta_deg in (0.1, 1.0, 8.5, 45.0, 90.0):
        th = np.radians(theta_deg)
        # orbit 1 in the equatorial plane, orbit 2 rotated by theta about x
        r1, v1, _ = state_from_elements(a, 0.0, 0.0, 0.0)
        r2, v2, _ = state_from_elements(a, 0.0, th, 0.0)
        assert np.allclose(r1, r2, rtol=1e-9), "states must share the node"
        w = v1 - v2
        u_t = v2 / np.linalg.norm(v2)
        w_mag = float(np.linalg.norm(w))
        retro = float(-np.dot(w, u_t))
        out.append(dict(
            theta_deg=theta_deg,
            w_3d_m_s=w_mag, w_closed_form=float(v_rel_crossing(V_ORB, th)),
            rel_err_w=abs(w_mag - float(v_rel_crossing(V_ORB, th)))
            / float(v_rel_crossing(V_ORB, th)),
            cos_psi_3d=retro / w_mag,
            cos_psi_closed_form=float(retrograde_efficiency(th)),
            rel_err_cos_psi=abs(retro / w_mag - float(retrograde_efficiency(th)))
            / float(retrograde_efficiency(th))))
    return pd.DataFrame(out)


def run():
    banner("INT-0  Engagement geometry and the co-orbital premise")

    val = validation_3d()
    save_table(val, "int00_validation_3d",
               "Closed-form |w| and retrograde efficiency checked against "
               "explicit 3-D inertial states at a common node.")
    print("  3-D validation, max relative error: "
          f"|w| {val.rel_err_w.max():.2e}, cos(psi) {val.rel_err_cos_psi.max():.2e}")

    dfa, dfc = part_a_kinematics()
    save_table(dfa, "int00_kinematics",
               "Retrograde efficiency cos(psi) and useful-delta-v per unit "
               "specific energy, co-orbital vs. crossing.",
               tags={"cos_psi": "DERIVED", "beta": "1.15 from INT-2"})
    save_table(dfc, "int00_crossing_energy_cost",
               "Projectile mass and specific energy required to deliver the "
               "SAME 51.8 m/s retrograde kick to a 1 cm target, as a function "
               "of approach geometry.")
    print("\n  Cost of delivering the same retrograde kick vs. approach angle:")
    print(dfc[["theta_deg", "w_m_s", "cos_psi", "m_p_required_g", "Es_kJ_kg",
               "times_threshold", "catastrophic"]].to_string(index=False))

    df2 = part_a2_combined()
    save_table(df2, "int00_combined_geometry",
               "The operational case: a nearly-co-orbital platform whose plane "
               "differs from the target's by theta. Specific energy for a fixed "
               "retrograde kick grows as 1 + (v_geom/v_launch)^2.",
               tags={"combined geometry": "DERIVED"})
    ok = df2[(df2.v_launch_m_s == 619.0) & df2.subcatastrophic]
    theta_cone = float(ok.theta_deg.max())
    ok10 = df2[(df2.v_launch_m_s == 619.0) & (df2.energy_penalty <= 1.10)]
    theta_cone10 = float(ok10.theta_deg.max())
    print("\n  Combined geometry at v_launch = 619 m/s:")
    print(df2[df2.v_launch_m_s == 619.0][
        ["theta_deg", "v_geom_m_s", "w_total_m_s", "Es_kJ_kg", "energy_penalty",
         "times_threshold", "subcatastrophic"]].to_string(index=False))
    print(f"  -> sub-catastrophic out to theta = {theta_cone:g} deg; "
          f"<=10% energy penalty out to {theta_cone10:g} deg")

    dfb = part_b_guidance()
    save_table(dfb, "int00_guidance_sweep",
               "Miss distance vs. approach geometry, engagement range and "
               "control latency, under a no-lead and a lead-computing aim law.",
               tags={"sigma_pos": "DERIVED (diffraction/SNR)"})

    dft = part_c_timeline()
    save_table(dft, "int00_timeline",
               "Engagement window: how long a target of a given transverse rate "
               "remains inside tracking range.")
    print("\n  Engagement window at 10 km tracking range:")
    print(dft[dft.r_track_m == 10000.0][
        ["theta_deg", "v_transverse_m_s", "window_s", "feasible_1s_sequence",
         "aim_slew_rate_deg_s"]].to_string(index=False))

    # ---- the honest boundary -------------------------------------------------
    # cone within which a no-lead aim law still lands inside a 1 cm target
    m_t = float(mass_sphere(0.01, RHO_AL))
    r_t = 0.005
    bound = []
    for tau in TAUS:
        for R in RANGES:
            tof = R / 619.0
            v_T_max = r_t / (tau + tof)
            th_max = 2 * np.arcsin(np.clip(v_T_max / (2 * V_ORB), 0, 1))
            bound.append(dict(tau_ms=tau * 1e3, R_m=R, tof_s=tof,
                              v_T_max_no_lead_m_s=v_T_max,
                              theta_max_no_lead_deg=np.degrees(th_max),
                              dh_offset_equivalent_m=v_T_max * 2 * (R_E + H)
                              / V_ORB))
    dfbound = pd.DataFrame(bound)
    save_table(dfbound, "int00_no_lead_cone",
               "The honest boundary: maximum transverse rate (and equivalent "
               "plane angle) for which a NON-extrapolating aim law still lands "
               "inside a 1 cm target.")
    print("\n  No-lead aim law: maximum tolerable transverse rate")
    print(dfbound[dfbound.R_m == 10.0][
        ["tau_ms", "v_T_max_no_lead_m_s", "theta_max_no_lead_deg"]].to_string(index=False))

    co = dfb[(dfb.theta_deg == 0.0) & (dfb.R_m == 10.0) & (dfb.dh_offset_m == 1000.0)]
    cross = dfb[(dfb.theta_deg == 8.5) & (dfb.R_m == 10.0)]
    concl = {
        "part_A_kinematic_verdict": {
            "retrograde_efficiency_crossing": "cos(psi) = sin(theta/2)",
            "useful_dv_per_Es_crossing": "beta / v_orbital  (independent of theta)",
            "useful_dv_per_Es_coorbital": "2 beta / |w|",
            "coorbital_advantage_factor": float(2 * V_ORB / 619.0),
            "statement": (
                f"A crossing engagement wastes all but sin(theta/2) of its "
                f"relative velocity on a cross-track impulse that changes "
                f"inclination rather than altitude, while paying the full |w|^2 "
                f"in disruption energy. Delivering the same 51.8 m/s retrograde "
                f"kick from an 8.5 deg crossing needs "
                f"{float(dfc[dfc.theta_deg == 8.5].m_p_required_g.iloc[0]):.3f} g "
                f"at {float(dfc[dfc.theta_deg == 8.5].Es_kJ_kg.iloc[0]):.0f} kJ/kg "
                f"({float(dfc[dfc.theta_deg == 8.5].times_threshold.iloc[0]):.0f}x "
                f"the catastrophic threshold) versus "
                f"{float(dfc[dfc.theta_deg == 0.0].m_p_required_g.iloc[0]):.3f} g at "
                f"{float(dfc[dfc.theta_deg == 0.0].Es_kJ_kg.iloc[0]):.1f} kJ/kg "
                f"co-orbital. The co-orbital architecture is forced by "
                f"kinematics alone, {2*V_ORB/619.0:.0f}x more efficient in "
                f"useful delta-v per unit disruption energy -- before any "
                f"guidance argument is made.")},
        "part_B_guidance_verdict": {
            "coorbital_miss_no_lead_m": float(co.miss_no_lead_m.min()),
            "coorbital_miss_with_lead_m": float(co.miss_with_lead_m.min()),
            "crossing_8p5deg_miss_no_lead_m": float(cross.miss_no_lead_m.min()),
            "crossing_8p5deg_miss_with_lead_m": float(cross.miss_with_lead_m.min()),
            "statement": (
                "The latency-independent miss is NOT a property of the guidance "
                "law -- it is a property of the geometry. Co-orbital, the "
                "transverse rate is the drift rate between neighbouring orbits "
                "(0.5 m/s at a 1 km altitude offset), so even a non-extrapolating "
                "aim law lands inside a 1 cm target at any latency tested. "
                "Crossing, the transverse rate is the full relative velocity and "
                "a non-extrapolating law misses by metres. An extrapolating law "
                "recovers the crossing case on paper, but see Part C.")},
        "part_C_timeline_verdict": {
            "statement": (
                "The binding constraint on crossing engagements is not the "
                "control loop but the engagement window. At 8.5 deg the target "
                f"crosses at {float(v_rel_crossing(V_ORB, np.radians(8.5))):.0f} m/s "
                "and is inside a 10 km tracking sphere for "
                f"{float(engagement_window(v_rel_crossing(V_ORB, np.radians(8.5)), 10000.0)):.1f} s, "
                "with the aim point slewing at "
                f"{float(np.degrees(v_rel_crossing(V_ORB, np.radians(8.5))/1000.0)):.0f} deg/s. "
                "Detection, track formation, solution and firing must all fit "
                "inside that, on a 1 cm unlit target, with one attempt per "
                "conjunction. Co-orbital, the same target is available "
                "indefinitely and the shot can be repeated.")},
        "premise_3_verdict": (
            "UPHELD, with the boundary stated: the latency-independent 5.8 mm "
            "miss holds for co-orbital engagements, which the plan's own energy "
            "gate independently forces. It does NOT hold for crossing "
            "engagements -- but crossing engagements are ruled out on kinematics "
            "(Part A) before guidance is even considered, so nothing downstream "
            "depends on the disputed claim. The correct statement is not 'lag "
            "does not matter' but 'the only geometry that transfers momentum "
            "efficiently is also the geometry in which lag does not matter'."),
        "firing_constraint_for_ECO": {
            "theta_max_subcatastrophic_deg": theta_cone,
            "theta_max_10pct_energy_penalty_deg": theta_cone10,
            "v_geom_max_m_s": float(v_rel_crossing(V_ORB, np.radians(theta_cone10))),
            "plane_match_tolerance_deg": theta_cone10,
            "note": ("ECO-1 may count a target as engageable only if the "
                     "platform's orbital plane is within this cone of the "
                     "target's, so that the launcher rather than the geometry "
                     "supplies the impact velocity. This is the binding "
                     "operational constraint the economics must respect, and it "
                     "is far tighter than the naive 'safe velocity ceiling' "
                     "reading of the v5 energy gate."),
        },
    }
    save_json(concl, "int00_conclusions")
    print("\n  co-orbital advantage factor: "
          f"{concl['part_A_kinematic_verdict']['coorbital_advantage_factor']:.1f}x")
    print("  plane-match tolerance handed to ECO-1: "
          f"{theta_cone10:g} deg (<=10% energy penalty), "
          f"{theta_cone:g} deg (sub-catastrophic limit)")

    # ---- figures -------------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2))
    ax = axes[0]
    th = np.linspace(0.05, 90, 400)
    w = v_rel_crossing(V_ORB, np.radians(th))
    ax.semilogy(th, retrograde_efficiency(np.radians(th)) * 100, color=CB[0], lw=2,
                label="retrograde fraction cos ψ")
    ax.axhline(100, color=CB[1], ls="--", lw=1.6, label="co-orbital (100%)")
    ax.set_xlabel("plane angle between platform and target (deg)")
    ax.set_ylabel("useful fraction of relative velocity (%)")
    ax.set_title("Only sin(θ/2) of a crossing shot is useful")
    ax.legend(fontsize=8)
    ax = axes[1]
    es = []
    for t in th:
        thr = np.radians(t)
        ww = float(v_rel_crossing(V_ORB, thr))
        cp = float(retrograde_efficiency(thr))
        m_t = float(mass_sphere(0.01, RHO_AL))
        m_p = float(dv_shift(H, 200e3)) * m_t / (1.15 * ww * cp)
        es.append(0.5 * m_p * ww ** 2 / m_t / 1e3)
    ax.semilogy(th, es, color=CB[0], lw=2, label="crossing engagement")
    ax.axhline(40, color="k", lw=1.3)
    ax.axhspan(40, 1e5, color="#D55E00", alpha=0.10)
    m_t = float(mass_sphere(0.01, RHO_AL))
    m_p_co = float(dv_shift(H, 200e3)) * m_t / (1.15 * 619.0)
    ax.axhline(0.5 * m_p_co * 619 ** 2 / m_t / 1e3, color=CB[2], ls="--", lw=1.6,
               label="co-orbital at 619 m/s")
    ax.text(2, 47, "catastrophic threshold", fontsize=7.5)
    ax.set_xlabel("plane angle between platform and target (deg)")
    ax.set_ylabel("specific energy for the same kick (kJ/kg)")
    ax.set_title("Crossing shots are catastrophic at any angle")
    ax.legend(fontsize=8)
    fig.suptitle("INT-0: the co-orbital architecture is forced by kinematics, "
                 "not chosen for guidance convenience", fontsize=10)
    save_fig(fig, "int00_kinematics",
             "Retrograde efficiency and disruption cost versus approach "
             "geometry. A crossing engagement wastes all but sin(θ/2) of its "
             "relative velocity on a useless cross-track impulse while paying "
             "the full |w|² in specific energy, so it exceeds the catastrophic "
             "threshold at every plane angle.")

    fig, ax = plt.subplots(figsize=(7.4, 4.3))
    for i, tau in enumerate([1e-4, 1e-3, 1e-2, 1e-1]):
        s = dfb[(dfb.tau_s == tau) & (dfb.R_m == 10.0)].copy()
        s = s.sort_values("v_transverse_m_s")
        ax.loglog(s.v_transverse_m_s, s.miss_no_lead_m, color=CB[i], marker="o",
                  ms=3, label=f"no lead, τ={tau*1e3:g} ms")
    ax.axhline(0.005, color="k", ls=":", lw=1.4)
    ax.text(1e-2, 0.0062, "1 cm target radius", fontsize=7.5)
    ax.axvspan(0.05, 5, color=CB[2], alpha=0.12)
    ax.text(0.09, 2e-5, "co-orbital\nregime", fontsize=7.5, color=CB[2])
    ax.set_xlabel("target transverse rate (m/s)")
    ax.set_ylabel("miss distance (m)")
    ax.set_title("INT-0: where the latency term does and does not bite (R = 10 m)")
    ax.legend(fontsize=7.5)
    save_fig(fig, "int00_miss_vs_transverse",
             "Miss distance under a non-extrapolating aim law against target "
             "transverse rate. In the co-orbital regime the miss stays inside "
             "the target at every latency tested; the latency term only bites "
             "for geometries that INT-0 Part A has already excluded.")
    return dfa, dfb, dft, concl


if __name__ == "__main__":
    run()
