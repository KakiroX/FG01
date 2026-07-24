"""
LAS-2 -- Accuracy requirements, head to head.

The claim under test: FG01 requires less accuracy than laser ablation, and is
therefore more robust than the cost comparison alone suggests.

The claim is TRUE, but not for the reason usually given, and one popular version
of it is FALSE. Both are reported.

WHAT IS FALSE. "The laser needs metre-level orbit prediction; FG01 needs less."
In *linear* terms FG01 needs the tighter number: it must place a ~9 mm cloud on
a 1 cm target, against the laser's ~25 cm spot radius. FG01's linear requirement
is roughly 28x tighter.

WHAT IS TRUE, and decisive. Accuracy is an *angular* problem, and the two
systems work at ranges differing by five orders of magnitude. The right metric
is the required pointing accuracy measured against what the sensor can
physically deliver -- its own diffraction limit. On that metric FG01 has ~130x
margin and the laser has ~1.6x, i.e. the laser must operate essentially at the
physical limit of its optics while FG01 has two orders of magnitude of headroom.

Four independent accuracy findings are quantified here:
  A. angular requirement vs. sensor diffraction limit (the margin metric)
  B. resolved vs. unresolved target, and the photon budget
  C. external cueing accuracy required before the engagement can begin
  D. graceful degradation vs. threshold cliff -- with the honest counterpoint
     that a laser firing many pulses averages over its own cliff
  E. the ablation-feedback problem: the laser perturbs the orbit it is trying
     to track, and must re-converge within ~10 pulses
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from fg01.constants import (R_E, RHO_AL, SEED, OPTICAL_LAMBDA_M,
                            OPTICAL_APERTURE_M)
from fg01.orbital import mass_sphere, xsec_sphere, dv_shift
from fg01.laser import (spot_diameter, atmospheric_transmission, fried_parameter,
                        ao_strehl, impulse_per_pulse, pass_geometry,
                        apparent_magnitude, photon_rate, resolved_ratio,
                        point_ahead_angle, F_TH, C_LIGHT)
from fg01.interaction import mass_efficiency, sigma_from_cloud_diameter
from fg01.io_utils import save_table, save_fig, save_json, banner, CB

H = 900e3
D_TARGET = 0.01
M_T = float(mass_sphere(D_TARGET, RHO_AL))
A_T = float(xsec_sphere(D_TARGET))
R_T = D_TARGET / 2
DV_REQ = float(dv_shift(H, 200e3))

# FG01 design point (v6)
FG_RANGE = 10.0
FG_SIGMA_CLOUD = 0.00876          # m, INT-3/4
FG_POINTING = 2e-4                # rad, design
FG_MISS = 0.005775                # m, achieved miss
FG_U = 10.0                       # m/s platform-target relative speed
FG_VREL = 619.0

# Laser design point (LAS-1 optimum)
LAS_RANGE = 1.186e6
LAS_APERTURE = 5.0
LAS_LAM = 532e-9
LAS_EPULSE = 20000.0
LAS_SPOT = 0.505
LAS_RETRO = 0.608
LAS_PRF = 10.0
LAS_PULSES = 241.0


def run():
    banner("LAS-2  Accuracy requirements, head to head")

    # ---- A. angular requirement vs. what the sensor can deliver -----------
    fg_req_ang = FG_SIGMA_CLOUD / FG_RANGE
    fg_diff = 1.22 * OPTICAL_LAMBDA_M / OPTICAL_APERTURE_M
    las_req_ang = (LAS_SPOT / 2) / LAS_RANGE
    las_diff = 1.22 * LAS_LAM / LAS_APERTURE

    rows = [
        dict(system="FG01 (co-orbital kick)", range_m=FG_RANGE,
             linear_tolerance_m=FG_SIGMA_CLOUD,
             angular_tolerance_rad=fg_req_ang,
             angular_tolerance_urad=fg_req_ang * 1e6,
             sensor_aperture_m=OPTICAL_APERTURE_M,
             diffraction_limit_urad=fg_diff * 1e6,
             margin_over_diffraction=fg_req_ang / fg_diff,
             design_pointing_urad=FG_POINTING * 1e6,
             achieved_miss_m=FG_MISS),
        dict(system="Ground laser ablation", range_m=LAS_RANGE,
             linear_tolerance_m=LAS_SPOT / 2,
             angular_tolerance_rad=las_req_ang,
             angular_tolerance_urad=las_req_ang * 1e6,
             sensor_aperture_m=LAS_APERTURE,
             diffraction_limit_urad=las_diff * 1e6,
             margin_over_diffraction=las_req_ang / las_diff,
             design_pointing_urad=las_req_ang * 1e6 / 3,
             achieved_miss_m=np.nan),
    ]
    dfa = pd.DataFrame(rows)
    save_table(dfa, "las02_angular_requirement",
               "Pointing requirement measured against the diffraction limit of "
               "each system's own sensing aperture. The margin column is the "
               "meaningful comparison: it is the headroom each system has "
               "against the physics of its own optics.",
               tags={"all": "DERIVED"})
    print(dfa[["system", "range_m", "linear_tolerance_m", "angular_tolerance_urad",
               "diffraction_limit_urad", "margin_over_diffraction"]].to_string(index=False))
    print(f"\n  ANGULAR tolerance ratio (laser is harder by): "
          f"{fg_req_ang/las_req_ang:,.0f}x")
    print(f"  LINEAR tolerance ratio (FG01 is harder by): "
          f"{(LAS_SPOT/2)/FG_SIGMA_CLOUD:,.1f}x  <-- the honest counterpoint")
    print(f"  MARGIN over own diffraction limit: FG01 {fg_req_ang/fg_diff:,.0f}x, "
          f"laser {las_req_ang/las_diff:,.2f}x")

    # ---- B. resolved vs unresolved, and the photon budget -----------------
    rows = []
    for lab, R, D, lam in (("FG01 terminal sensor", FG_RANGE,
                            OPTICAL_APERTURE_M, OPTICAL_LAMBDA_M),
                           ("Ground telescope", LAS_RANGE, LAS_APERTURE, 550e-9),
                           ("Ground telescope, 10 m", LAS_RANGE, 10.0, 550e-9)):
        rr = float(resolved_ratio(D_TARGET, R, D, lam))
        pr = float(photon_rate(D_TARGET, R, D))
        rows.append(dict(
            sensor=lab, range_m=R, aperture_m=D,
            target_angular_size_urad=D_TARGET / R * 1e6,
            diffraction_limit_urad=1.22 * lam / D * 1e6,
            resolved_ratio=rr,
            resolved=bool(rr > 1),
            apparent_magnitude=float(apparent_magnitude(D_TARGET, R)),
            photoelectrons_per_s=pr,
            photoelectrons_per_10ms=pr * 0.01,
            snr_10ms=pr * 0.01 / np.sqrt(max(pr * 0.01, 1.0) + 100.0)))
    dfb = pd.DataFrame(rows)
    save_table(dfb, "las02_detection",
               "Detection geometry and photon budget for a 1 cm, albedo-0.1 "
               "sphere. 'resolved_ratio' is target angular size divided by the "
               "diffraction limit: above 1 the target is an extended object "
               "whose centroid is directly measurable; below 1 it is an "
               "unresolved point source.",
               tags={"photometry": "DERIVED", "albedo 0.1": "UNVALIDATED"})
    print("\n  detection:")
    print(dfb[["sensor", "target_angular_size_urad", "diffraction_limit_urad",
               "resolved_ratio", "apparent_magnitude", "photoelectrons_per_s",
               "snr_10ms"]].to_string(index=False))

    # ---- C. external cueing required to begin the engagement --------------
    # FG01 must be cued well enough to close to its own sensor's acquisition
    # range; the laser must be cued to within its own spot at full range.
    # Cueing is evaluated at the range at which the cue is USED: FG01 uses it to
    # close from a rendezvous approach (~50 km) to its own sensor's acquisition
    # range, after which its onboard sensor closes the loop; the laser must use
    # it at full slant range because it has no means of improving on it.
    cue = []
    for lab, need_m, R in (
            ("FG01: close to optical acquisition range", 1000.0, 50e3),
            ("Ground laser: place spot on target", LAS_SPOT / 2, LAS_RANGE)):
        cue.append(dict(system=lab, required_position_knowledge_m=need_m,
                        cue_used_at_range_m=R,
                        angular_equivalent_urad=need_m / R * 1e6,
                        loop_closed_onboard=bool(lab.startswith("FG01"))))
    dfc = pd.DataFrame(cue)
    save_table(dfc, "las02_cueing",
               "External orbit knowledge each system needs before an engagement "
               "can start. FG01 uses external cueing only to get within range "
               "of its own sensor, which then closes the loop; the laser must "
               "be handed a solution accurate to its spot radius at 1,186 km "
               "because it has no way to improve on it.")
    print("\n  external cueing requirement:")
    print(dfc.to_string(index=False))
    print(f"  -> the laser needs the catalogue to be "
          f"{1000.0/(LAS_SPOT/2):,.0f}x more accurate than FG01 does")

    # ---- D. graceful vs cliff --------------------------------------------
    # FG01: delivered mass fraction falls smoothly with miss distance.
    # Laser: fluence at the target falls as a Gaussian; once it drops below the
    # ablation threshold the delivered impulse is exactly zero.
    R_las, elev, retro = pass_geometry(H, np.radians(-7.0))
    T_atm = float(atmospheric_transmission(elev))
    w_beam = LAS_SPOT / 2                       # 1/e^2 beam radius
    _, F_peak, _ = impulse_per_pulse(LAS_EPULSE, LAS_SPOT, A_T, T_atm)
    F_peak = float(F_peak)

    deg = []
    for x in np.geomspace(0.05, 20.0, 220):
        # FG01: miss = x * nominal miss
        d_fg = x * FG_MISS
        eps = float(mass_efficiency(FG_SIGMA_CLOUD, R_T, d_fg))
        eps0 = float(mass_efficiency(FG_SIGMA_CLOUD, R_T, 0.0))
        # laser: pointing error = x * nominal tolerance (spot radius)
        d_las = x * w_beam
        F_at = F_peak * np.exp(-2 * d_las ** 2 / w_beam ** 2)
        las_frac = (F_at / F_peak) if F_at >= F_TH else 0.0
        deg.append(dict(
            error_multiple=x,
            fg01_delivered_fraction=eps / eps0,
            fg01_miss_mm=d_fg * 1e3,
            laser_pointing_error_m=d_las,
            laser_fluence_J_cm2=F_at / 1e4,
            laser_above_threshold=bool(F_at >= F_TH),
            laser_delivered_fraction=las_frac))
    dfd = pd.DataFrame(deg)
    save_table(dfd, "las02_degradation",
               "Delivered momentum as each system's aiming error grows, "
               "normalised to its own nominal tolerance. FG01 decays smoothly; "
               "the laser holds up and then goes to exactly zero when the "
               "fluence crosses the ablation threshold.")
    fg_at2 = float(np.interp(2.0, dfd.error_multiple, dfd.fg01_delivered_fraction))
    las_cliff = float(dfd[dfd.laser_delivered_fraction > 0].error_multiple.max())
    print(f"\n  at 2x the nominal aiming error: FG01 still delivers "
          f"{fg_at2*100:.1f}% of nominal")
    print(f"  laser delivers zero beyond {las_cliff:.2f}x its nominal tolerance "
          f"(fluence falls under the {F_TH/1e4:.0f} J/cm² ablation threshold)")

    # ---- honest counterpoint: the laser averages over its own cliff -------
    rng = np.random.default_rng(SEED)
    cp = []
    for sig_mult in (0.25, 0.5, 1.0, 2.0, 4.0):
        n = 200000
        d = sig_mult * w_beam * np.sqrt(rng.chisquare(2, n))
        F_at = F_peak * np.exp(-2 * d ** 2 / w_beam ** 2)
        hit = F_at >= F_TH
        frac = np.where(hit, F_at / F_peak, 0.0)
        # FG01 with the same relative pointing error
        d_fg = sig_mult * FG_MISS * np.sqrt(rng.chisquare(2, n))
        eps = mass_efficiency(FG_SIGMA_CLOUD, R_T, d_fg)
        eps0 = float(mass_efficiency(FG_SIGMA_CLOUD, R_T, 0.0))
        cp.append(dict(
            pointing_sigma_multiple=sig_mult,
            laser_pulse_hit_rate=float(hit.mean()),
            laser_mean_delivered_fraction=float(frac.mean()),
            laser_pulses_needed=LAS_PULSES / max(frac.mean(), 1e-9),
            laser_engagement_time_s=LAS_PULSES / max(frac.mean(), 1e-9) / LAS_PRF,
            fg01_mean_delivered_fraction=float((eps / eps0).mean()),
            fg01_shots_needed=1.0 / max(float((eps / eps0).mean()), 1e-9)))
    dfe = pd.DataFrame(cp)
    save_table(dfe, "las02_repetition_counterpoint",
               "The honest counterpoint to the cliff argument: a laser firing "
               "hundreds of pulses averages over its own threshold, so the "
               "cliff degrades throughput rather than causing outright failure "
               "-- until the pointing error grows enough that the engagement "
               "no longer fits inside a pass.")
    print("\n  effect of pointing error on each system:")
    print(dfe.to_string(index=False))

    # ---- E. the ablation-feedback problem ---------------------------------
    dv_pulse = DV_REQ / LAS_PULSES
    fb = []
    for n_unmodelled in (1, 3, 10, 30, 100):
        dt = n_unmodelled / LAS_PRF
        drift = 0.5 * dv_pulse * n_unmodelled * dt      # accumulated displacement
        fb.append(dict(
            pulses_without_feedback=n_unmodelled,
            elapsed_s=dt,
            unmodelled_displacement_m=drift,
            fraction_of_spot_radius=drift / w_beam,
            target_lost=bool(drift > w_beam)))
    dff = pd.DataFrame(fb)
    save_table(dff, "las02_ablation_feedback",
               "The coupled ablation-and-estimation problem: each pulse changes "
               "the target's velocity, so an estimator that does not model its "
               "own effect walks off the target.")
    print("\n  ablation feedback (laser only):")
    print(dff.to_string(index=False))
    lost_at = dff[dff.target_lost].pulses_without_feedback.min() if dff.target_lost.any() else None

    # ---- point-ahead -------------------------------------------------------
    v_t_las = np.sqrt(3.986004418e14 / (R_E + H))
    pa_las = float(point_ahead_angle(v_t_las, R_las))
    pa_fg = FG_U * (FG_RANGE / FG_VREL) / FG_RANGE
    pah = pd.DataFrame([
        dict(system="Ground laser", point_ahead_urad=pa_las * 1e6,
             tolerance_urad=las_req_ang * 1e6,
             fractional_accuracy_needed=las_req_ang / pa_las),
        dict(system="FG01", point_ahead_urad=pa_fg * 1e6,
             tolerance_urad=fg_req_ang * 1e6,
             fractional_accuracy_needed=fg_req_ang / pa_fg),
    ])
    save_table(pah, "las02_point_ahead",
               "Point-ahead angle and the fractional accuracy with which it "
               "must be computed. The laser's lead is small but its tolerance "
               "is smaller; FG01's lead is large but its tolerance is larger.")
    print("\n  point-ahead:")
    print(pah.to_string(index=False))

    concl = {
        "headline": (
            "FG01's accuracy advantage is real and is best stated as headroom "
            "against the physical limit of its own optics: FG01's pointing "
            f"requirement sits {fg_req_ang/fg_diff:,.0f}x above the diffraction "
            "limit of its 10 cm sensor, while the ground laser's sits "
            f"{las_req_ang/las_diff:,.2f}x above the diffraction limit of a 5 m "
            "telescope. The laser must work essentially at the physical limit "
            "of its aperture; FG01 has two orders of magnitude of margin."),
        "angular_vs_linear": {
            "angular_ratio_laser_harder": float(fg_req_ang / las_req_ang),
            "linear_ratio_fg01_harder": float((LAS_SPOT / 2) / FG_SIGMA_CLOUD),
            "note": (
                "Both are true and they point opposite ways. In linear terms "
                "FG01 is the tighter requirement (9 mm against 25 cm). The "
                "angular framing is the operative one because pointing hardware, "
                "diffraction and atmospheric turbulence all act on angles, not "
                "on metres. Any argument for FG01 that quotes the linear number "
                "is wrong and will be caught.")},
        "resolved_target": {
            "fg01_resolved_ratio": float(dfb[dfb.sensor.str.startswith("FG01")].resolved_ratio.iloc[0]),
            "laser_resolved_ratio": float(dfb[dfb.sensor == "Ground telescope"].resolved_ratio.iloc[0]),
            "statement": (
                "FG01 sees the target as an extended object resolved across "
                f"{float(dfb[dfb.sensor.str.startswith('FG01')].resolved_ratio.iloc[0]):.0f} "
                "diffraction widths, with ~10^11 photoelectrons per second. The "
                "ground telescope sees an unresolved point source "
                f"{1/float(dfb[dfb.sensor == 'Ground telescope'].resolved_ratio.iloc[0]):.0f}x "
                "below its own resolution at magnitude "
                f"{float(dfb[dfb.sensor == 'Ground telescope'].apparent_magnitude.iloc[0]):.1f}. "
                "Detection is feasible -- this is not the objection sometimes "
                "made -- but every position measurement is a centroid of a "
                "point source rather than a direct measurement of an extended "
                "one.")},
        "cueing_ratio": float(1000.0 / (LAS_SPOT / 2)),
        "graceful_vs_cliff": {
            "fg01_delivered_at_2x_error": fg_at2,
            "laser_zero_beyond_multiple": las_cliff,
            "statement": (
                f"At twice its nominal aiming error FG01 still delivers "
                f"{fg_at2*100:.1f}% of nominal momentum; the laser delivers "
                f"exactly zero beyond {las_cliff:.2f}x its tolerance, because "
                "fluence crosses the ablation threshold. The failure modes "
                "differ in kind, not degree.")},
        "honest_counterpoint": (
            "The cliff is per-pulse, and the laser fires hundreds of pulses per "
            "object, so it averages over its own threshold: at 1x pointing "
            f"sigma the pulse hit rate is "
            f"{float(dfe[dfe.pointing_sigma_multiple == 1.0].laser_pulse_hit_rate.iloc[0]):.2f} "
            "and the engagement simply takes longer. The cliff therefore "
            "degrades throughput rather than causing outright mission failure, "
            "and any argument that treats a laser miss as total loss overstates "
            "the case. What does NOT average out is the tracking requirement: "
            "FG01 needs ONE aiming solution, the laser needs "
            f"{LAS_PULSES:.0f} consecutive ones."),
        "ablation_feedback": {
            "dv_per_pulse_m_s": float(dv_pulse),
            "pulses_before_target_lost": (int(lost_at) if lost_at else None),
            "statement": (
                f"Each pulse changes the target's velocity by {dv_pulse*1e3:.0f} "
                "mm/s, so an estimator that does not model the effect of its own "
                f"shots walks off the target within {lost_at} pulses "
                f"({(lost_at or 0)/LAS_PRF:.1f} s). The laser must therefore "
                "solve a coupled ablation-and-estimation problem in real time "
                "for the whole engagement. FG01 fires once and has no such "
                "coupling.")},
        "summary_of_the_asymmetry": (
            "FG01: one aiming solution, 200 urad, on a resolved target 10 m "
            "away, in vacuum, with graceful degradation and no feedback loop. "
            "Ground laser: 241 consecutive aiming solutions, 0.21 urad, on an "
            "unresolved point source 1,186 km away, through turbulent "
            "atmosphere, with a hard threshold and a self-perturbing target. "
            "That is the accuracy case for FG01, and it is a strong one -- but "
            "it is a case about engineering risk, not about cost, and LAS-3 "
            "shows why that distinction matters."),
    }
    save_json(concl, "las02_conclusions")

    # ---- figures -----------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.3))
    ax = axes[0]
    ax.semilogx(dfd.error_multiple, dfd.fg01_delivered_fraction * 100,
                color=CB[0], lw=2.2, label="FG01 — momentum delivered")
    ax.semilogx(dfd.error_multiple, dfd.laser_delivered_fraction * 100,
                color=CB[1], lw=2.2, label="Laser — momentum delivered")
    ax.axvline(1.0, color="0.5", ls=":", lw=1.2)
    ax.text(1.05, 60, "nominal\ntolerance", fontsize=7, color="0.35")
    ax.set_xlabel("aiming error / nominal tolerance")
    ax.set_ylabel("delivered momentum (% of nominal)")
    ax.set_title("Graceful decay vs. threshold cliff")
    ax.legend(fontsize=8)
    ax = axes[1]
    x = np.arange(2)
    ax.bar(x - 0.2, [fg_req_ang * 1e6, las_req_ang * 1e6], 0.4, color=CB[0],
           label="required pointing")
    ax.bar(x + 0.2, [fg_diff * 1e6, las_diff * 1e6], 0.4, color=CB[1],
           label="sensor diffraction limit")
    ax.set_yscale("log")
    ax.set_xticks(x); ax.set_xticklabels(["FG01\n(10 m)", "Ground laser\n(1,186 km)"],
                                         fontsize=8)
    ax.set_ylabel("angle (µrad)")
    for i, m in enumerate([fg_req_ang / fg_diff, las_req_ang / las_diff]):
        ax.text(i, max(fg_req_ang, las_req_ang) * 1e6 * 3, f"{m:,.0f}× margin"
                if m > 10 else f"{m:,.2f}× margin", ha="center", fontsize=8)
    ax.legend(fontsize=8)
    ax.set_title("Headroom against the physics of the optics")
    fig.suptitle("LAS-2: the accuracy asymmetry, stated in the units that matter",
                 fontsize=10)
    save_fig(fig, "las02_accuracy",
             "Left: delivered momentum against aiming error. FG01 decays "
             "smoothly because a partially-missed cloud still deposits mass; "
             "the laser goes to exactly zero once fluence drops below the "
             "ablation threshold. Right: each system's pointing requirement "
             "against the diffraction limit of its own aperture — FG01 works "
             "with two orders of magnitude of margin, the laser at the limit.")
    return dfa, dfd, concl


if __name__ == "__main__":
    run()
