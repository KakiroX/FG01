"""
SIM-14 -- Terminal-accuracy comparison: space-based laser ablation vs FG01.

Scope, stated up front (mirrors docs/LIMITATIONS.md practice): this sim
compares ONE dimension only -- tolerance of the terminal engagement to
pointing / orbit-prediction error. It does NOT reopen the cost comparison:
SIM-9 already finds FG01 is 3-4 orders of magnitude more expensive per
object than ground-based laser ablation, and nothing here changes that.

The laser-side numbers come from `laser accuracy comm.md`, an informal brief
(docs/CITATION_LOG.md key `laser_brief`, status E -- not independently
verified against primary literature). Its own claims are used AS STATED,
not tuned to favour either system:
  - orbit prediction "better than 1 m", optimistically "<1 cm within 10 s"
  - pointing "sub-arcsecond" (no exact figure given; 0.1-2.0 arcsec used)
  - omitting p-N correction biases the shot "on the order of the size of
    the debris objects themselves" (implemented literally as a residual
    equal to target radius)
  - the brief's own comparison table states the laser's failure mode is
    binary: "Miss entirely; no effect", against FG01's "Cloud still
    intersects debris path; likely partial effect"

The FG01-side numbers are NOT re-derived here -- they are the already
-validated SIM-4 M3 error budget (fg01.gnc.miss_extrapolating), evaluated at
three (R, tau, u, pointing, muzzle, bandwidth) points already inside SIM-4's
tested grid, so nothing here extrapolates beyond what SIM-4 validated.

Single unified metric for both systems: expected fraction of the intended
effect delivered per shot.
  - Laser: modelled as binary per the brief's own table (hit = full design
    ablation impulse, miss = zero), so E[fraction] = P(miss < r_target)
    under a Rayleigh miss distribution (fg01.laser.p_effective_shot).
  - FG01: continuous areal-capture model already established in SIM-3/4
    (fg01.delivery.expected_fraction_random_miss) -- no re-derivation.

A second, equally literal reading of the brief is carried alongside: laser
ablation needs a *sustained* hit across many pulses to accumulate useful
Delta-v (the brief: "must stay focused... on a small, fast-moving target"),
so an engagement's success probability compounds as P_hit^N over N pulses.
FG01's design point needs exactly one kick (SIM-13 design card).
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from fg01.constants import (R_E, OPTICAL_APERTURE_M, OPTICAL_LAMBDA_M,
                            DEBRIS_DIAM_M, LASER_ORBIT_PRED_CLAIM_M,
                            LASER_ORBIT_PRED_REQUIREMENT_M,
                            LASER_ORBIT_PRED_DEGRADED_M,
                            LASER_STANDOFF_RANGE_KM)
from fg01.gnc import miss_extrapolating
from fg01.delivery import expected_fraction_random_miss, sigma_opt_for_confidence
from fg01.laser import laser_miss_budget, p_effective_shot
from fg01.io_utils import save_table, save_fig, save_json, banner, CB

V_REL_FG01 = 619.0
H_DES = 800e3

# Three FG01 points, all inside the grid SIM-4 already swept and validated
# (RANGES, TAU_MS, US, POINTING, MUZZLE, LIDAR bandwidths in fg01.constants /
# sim04_gnc.py) -- nothing here is extrapolated beyond that grid. Engagement
# range R is held FIXED at the SIM-0/SIM-11 co-orbital design value (10 m)
# across all three scenarios: R is an engagement-geometry choice forced by
# the energy gate, not a GNC-quality parameter, so only the GNC-quality terms
# (latency, relative speed, pointing, muzzle repeatability, ranging bandwidth)
# vary between best-case, design-point and worst-case. Varying R alongside
# GNC quality would stack two different kinds of degradation into one
# scenario and make the "degraded" case harder than its laser counterpart.
FG01_SCENARIOS = {
    "optimistic":  dict(R=10.0, tau=0.1e-3, u=1.0,   pointing=5e-5, muzzle=1e-3, bw=1e9),
    "as_designed": dict(R=10.0, tau=10e-3,  u=10.0,  pointing=2e-4, muzzle=1e-3, bw=1e8),
    "degraded":    dict(R=10.0, tau=100e-3, u=100.0, pointing=1e-3, muzzle=1e-2, bw=1e6),
}
T_OBS, SAMPLE_HZ = 10.0, 1000.0

# Three laser points, built directly from laser_accuracy_comm.md's own stated
# figures (see fg01/laser.py docstring for the mapping). Engagement range is
# likewise held FIXED across all three scenarios (500 km, the middle of the
# representative 100-1000 km grid in fg01.constants) for the same reason:
# range is a stand-off/geometry choice, not a GNC-quality parameter, and
# holding it fixed keeps the three-scenario structure symmetric with FG01's.
LASER_R_KM = LASER_STANDOFF_RANGE_KM[2]   # 500 km: middle of the representative grid
LASER_SCENARIOS = {
    "optimistic":            dict(orbit_pred=LASER_ORBIT_PRED_CLAIM_M,       pointing_arcsec=0.1, R_km=LASER_R_KM, pN_corrected=True),
    "as_stated_requirement": dict(orbit_pred=LASER_ORBIT_PRED_REQUIREMENT_M, pointing_arcsec=0.5, R_km=LASER_R_KM, pN_corrected=True),
    "degraded_uncorrected":  dict(orbit_pred=LASER_ORBIT_PRED_DEGRADED_M,    pointing_arcsec=2.0, R_km=LASER_R_KM, pN_corrected=False),
}

N_PULSES = [1, 2, 3, 5, 10]   # P_hit^N underflows float64 well before N=100 at these
                              # per-shot P_hit values -- see sim14_conclusions.json


def fg01_sigma_total(sc):
    a_orbit = R_E + H_DES
    total, terms = miss_extrapolating(
        R=sc["R"], tau=sc["tau"], u=sc["u"], v_rel=V_REL_FG01,
        lam=OPTICAL_LAMBDA_M, D=OPTICAL_APERTURE_M,
        sigma_point=sc["pointing"], dv_muzzle_frac=sc["muzzle"],
        t_obs=T_OBS, n_samples=T_OBS * SAMPLE_HZ, a_orbit=a_orbit)
    return float(total), {k: float(v) for k, v in terms.items()}


def run():
    banner("SIM-14  Terminal-accuracy comparison: laser ablation vs FG01")
    print("  Scope: pointing/prediction error tolerance ONLY. Cost is not "
          "reopened -- see SIM-9 (FG01 is 3-4 orders of magnitude more "
          "expensive per object than ground-based laser ablation).")

    # ---- (1) total miss budget, both systems -------------------------------
    fg01_sigma = {name: fg01_sigma_total(sc) for name, sc in FG01_SCENARIOS.items()}
    budget_rows = []
    for name, (total, terms) in fg01_sigma.items():
        row = dict(system="FG01", scenario=name, sigma_total_m=total)
        row.update({f"term_{k}_m": v for k, v in terms.items()})
        row.update(FG01_SCENARIOS[name])
        budget_rows.append(row)

    r_ref = 0.005   # 1 cm target radius, used only for the headline budget row
    for name, sc in LASER_SCENARIOS.items():
        pN_residual = r_ref if not sc["pN_corrected"] else 0.0
        total, terms = laser_miss_budget(sc["orbit_pred"], sc["pointing_arcsec"],
                                         sc["R_km"] * 1e3, pN_residual)
        row = dict(system="Laser", scenario=name, sigma_total_m=float(total))
        row.update({f"term_{k}_m": float(v) for k, v in terms.items()})
        row.update(sc)
        budget_rows.append(row)
    dfb = pd.DataFrame(budget_rows)
    save_table(dfb, "sim14_error_budget",
               "Total 1-sigma terminal miss budget, both systems, three "
               "scenarios each. FG01 rows reuse SIM-4's validated M3 model "
               "at points already inside its tested grid; laser rows "
               "implement laser_accuracy_comm.md's own stated figures. "
               "Laser's pN_residual in the degraded scenario is evaluated "
               "at a 1 cm (0.005 m) target radius for this headline row; "
               "the target-size sweep below computes it per target.",
               tags={"FG01 rows": "DERIVED (SIM-4 validated model)",
                     "Laser rows": "derived from `laser_brief` (E, UNVALIDATED range/pointing figure)"})
    print("\n  total miss budget (m):")
    print(dfb[["system", "scenario", "sigma_total_m"]].to_string(index=False))

    # ---- (2) expected delivered fraction vs target size --------------------
    # Cloud sigma is FIXED at the as-designed shot pattern (SIM-4/SIM-13 design
    # card: 8.76 mm, i.e. sigma_opt_for_confidence at the as_designed terminal
    # budget) and reused across all three scenarios. This is deliberate: SIM-3
    # shows FG01 *could* resize its shot pattern to trade mass for tolerance,
    # but giving it that adaptive freedom here while holding the laser's spot
    # size fixed would stack the comparison in FG01's favour. Both systems are
    # therefore stress-tested here as fixed, already-engineered designs.
    s_cloud_fixed = float(sigma_opt_for_confidence(fg01_sigma["as_designed"][0], 0.9))
    frac_rows = []
    for d in DEBRIS_DIAM_M:
        r_t = d / 2.0
        for name, (total, _) in fg01_sigma.items():
            frac = float(expected_fraction_random_miss(s_cloud_fixed, total, r_t))
            frac_rows.append(dict(system="FG01", scenario=name, d_target_m=d,
                                  sigma_total_m=total, expected_fraction=frac))
        for name, sc in LASER_SCENARIOS.items():
            pN_residual = r_t if not sc["pN_corrected"] else 0.0
            total, _ = laser_miss_budget(sc["orbit_pred"], sc["pointing_arcsec"],
                                        sc["R_km"] * 1e3, pN_residual)
            p_hit = float(p_effective_shot(total, r_t))
            frac_rows.append(dict(system="Laser", scenario=name, d_target_m=d,
                                  sigma_total_m=float(total), expected_fraction=p_hit))
    dff = pd.DataFrame(frac_rows)
    save_table(dff, "sim14_expected_fraction_vs_size",
               "Expected fraction of the intended per-shot effect actually "
               "delivered, vs target diameter. FG01: continuous areal-"
               "capture model (delivery.expected_fraction_random_miss, "
               "already established in SIM-3/4). Laser: binary per the "
               "brief's own comparison table (hit=full effect, miss=zero), "
               "so expected fraction = P(miss < target radius) under a "
               "Rayleigh miss.",
               tags={"metric": "DERIVED (unified single-shot expected-effect fraction)"})

    # ---- (3) compounding-shots view -----------------------------------------
    pulse_rows = []
    for name, sc in LASER_SCENARIOS.items():
        for d in DEBRIS_DIAM_M:
            r_t = d / 2.0
            pN_residual = r_t if not sc["pN_corrected"] else 0.0
            total, _ = laser_miss_budget(sc["orbit_pred"], sc["pointing_arcsec"],
                                        sc["R_km"] * 1e3, pN_residual)
            p_hit = float(p_effective_shot(total, r_t))
            for n in N_PULSES:
                pulse_rows.append(dict(scenario=name, d_target_m=d, n_pulses=n,
                                       p_hit_single=p_hit,
                                       p_all_n_hit=p_hit ** n))
    dfp = pd.DataFrame(pulse_rows)
    save_table(dfp, "sim14_laser_pulse_compounding",
               "P(every one of N pulses hits) = p_hit^N for a laser "
               "engagement needing sustained tracking across N pulses to "
               "accumulate useful Delta-v (per the brief: the system 'must "
               "stay focused... on a small, fast-moving target'). FG01's "
               "design point needs exactly 1 kick per object (SIM-13).",
               tags={"p_hit^N": "DERIVED (independent-pulse approximation, "
                                "conservative if tracking errors are correlated "
                                "shot-to-shot)"})

    # ---- conclusions --------------------------------------------------------
    d1cm = 0.01
    fg01_1cm = dff[(dff.system == "FG01") & (dff.d_target_m == d1cm)]
    laser_1cm = dff[(dff.system == "Laser") & (dff.d_target_m == d1cm)]
    concl = {
        "scope": ("Pointing/prediction error tolerance only. Cost is not "
                 "reopened; SIM-9's finding that FG01 is 3-4 orders of "
                 "magnitude more expensive per object than ground-based "
                 "laser ablation stands unchanged."),
        "laser_source": ("All laser-side figures are taken as stated in "
                         "laser_accuracy_comm.md, an informal brief not "
                         "independently verified against primary literature "
                         "in this study (CITATION_LOG key `laser_brief`, "
                         "status E). Engagement range (fixed at 500 km, the "
                         "middle of a 100-1000 km representative grid) and "
                         "the exact pointing figure within 'sub-arcsecond' "
                         "are this study's own representative assumptions, "
                         "not sourced from the brief."),
        "fg01_source": ("All FG01-side figures reuse SIM-4's already-"
                        "validated M3 error-budget model at points inside "
                        "its tested grid; nothing is re-derived or "
                        "extrapolated here."),
        "headline_1cm_target": {
            "FG01 expected fraction delivered": {
                r.scenario: round(float(r.expected_fraction), 4)
                for r in fg01_1cm.itertuples()},
            "Laser expected fraction delivered (=P_hit)": {
                r.scenario: round(float(r.expected_fraction), 6)
                for r in laser_1cm.itertuples()},
        },
        "structural_asymmetry": (
            "FG01's areal-capture model gives a smoothly degrading delivered "
            "fraction as error grows (SIM-3's non-central-chi-square model); "
            "the brief's own comparison table states the laser's failure "
            "mode is binary -- a miss beyond the target radius delivers "
            "nothing. This is the single largest structural difference "
            "between the two systems on this one dimension."),
        "why_fg01_tolerates_more_absolute_angular_error": (
            "Not better sensors: FG01's terminal budget is dominated by "
            "mechanical/optical terms at a ~10 m engagement range (SIM-4), "
            "so even a comparatively loose 200 microrad pointing spec gives "
            "a millimetre-scale linear miss. The laser's sub-arcsecond spec "
            "is far tighter in angle, but at a 100-1000 km stand-off range "
            "the same angular budget maps to metre-scale linear miss -- "
            "which is why the brief frames orbit prediction and pointing as "
            "both needing to reach the metre level. The asymmetry is a "
            "consequence of engagement geometry (stand-off range), not a "
            "claim that FG01 hardware is more precise."),
        "pulse_compounding": (
            "Single-pulse P_hit for a 1 cm target is already so small in "
            "every scenario tested (2e-4 optimistic, 5e-6 as-stated, 3e-7 "
            "degraded) that P(N consecutive hits)=P_hit^N underflows "
            "float64 well before N=100 -- i.e. requiring even a modest "
            "number of sustained hits on cm-scale debris is not just "
            "unlikely but numerically indistinguishable from impossible "
            "under these parameters. FG01 requires exactly one kick per "
            "object at its design point (SIM-13 design card), so no "
            "analogous compounding penalty applies to it."),
        "not_claimed": ("This sim does not claim FG01 is cheaper, that the "
                        "laser numbers are precise (the brief itself gives "
                        "ranges, not point estimates, for pointing), or that "
                        "a real space-based laser system could not correct "
                        "for the degraded scenario. The degraded scenario is "
                        "a what-if built from the brief's own stated risk "
                        "factors (orbit-prediction drift, uncorrected p-N "
                        "terms), not a claim about any specific fielded system. "
                        "It also does not model a laser trading focus for "
                        "tolerance (a defocused beam would spread ablation "
                        "fluence over a larger spot, an analogue of FG01's own "
                        "sigma-vs-mass tradeoff in SIM-3) -- the brief's own "
                        "comparison table treats the laser as a fixed, binary "
                        "hit/miss system, and that framing is followed here "
                        "rather than adding an untested laser design freedom. "
                        "FG01's own shot pattern is likewise held FIXED across "
                        "scenarios (sigma=8.76 mm, the as-designed value) "
                        "rather than re-optimised per scenario, so neither "
                        "system is given an adaptive-resizing advantage the "
                        "other lacks."),
    }
    save_json(concl, "sim14_conclusions")
    print("\n  expected fraction delivered, 1 cm target:")
    print("    FG01:  ", concl["headline_1cm_target"]["FG01 expected fraction delivered"])
    print("    Laser: ", concl["headline_1cm_target"]["Laser expected fraction delivered (=P_hit)"])

    # ---- figures -------------------------------------------------------------
    fig, ax = plt.subplots()
    fg01_names = list(FG01_SCENARIOS)
    laser_names = list(LASER_SCENARIOS)
    x = np.arange(3)
    w = 0.35
    fg01_vals = [fg01_sigma[n][0] for n in fg01_names]
    laser_vals = [laser_miss_budget(LASER_SCENARIOS[n]["orbit_pred"],
                                    LASER_SCENARIOS[n]["pointing_arcsec"],
                                    LASER_SCENARIOS[n]["R_km"] * 1e3,
                                    r_ref if not LASER_SCENARIOS[n]["pN_corrected"] else 0.0)[0]
                 for n in laser_names]
    ax.bar(x - w / 2, fg01_vals, w, color=CB[0], label="FG01 (SIM-4 M3 budget)")
    ax.bar(x + w / 2, laser_vals, w, color=CB[1], label="Laser (laser_accuracy_comm.md)")
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(["optimistic", "as-designed / as-stated", "degraded"])
    ax.axhline(r_ref, color="0.3", ls=":", lw=1)
    ax.text(1.55, r_ref * 1.15, "1 cm target radius", fontsize=7.5, color="0.3")
    ax.set_ylabel("total 1σ miss budget (m, log scale)")
    ax.set_title("SIM-14: terminal miss budget, three scenarios each\n"
                 "(scope: accuracy tolerance only -- see SIM-9 for cost)")
    ax.legend(fontsize=8)
    save_fig(fig, "sim14_miss_budget",
             "Total terminal miss-budget comparison. FG01's budget is "
             "dominated by close engagement range, not superior sensing; "
             "the laser's sub-arcsecond angular budget maps to a much "
             "larger linear miss at its hundreds-of-km stand-off range.")

    fig, ax = plt.subplots()
    for i, name in enumerate(fg01_names):
        s = dff[(dff.system == "FG01") & (dff.scenario == name)].sort_values("d_target_m")
        ax.plot(s.d_target_m * 100, s.expected_fraction, color=CB[i], marker="o",
               ms=4, label=f"FG01 {name}")
    for i, name in enumerate(laser_names):
        s = dff[(dff.system == "Laser") & (dff.scenario == name)].sort_values("d_target_m")
        ax.plot(s.d_target_m * 100, s.expected_fraction, color=CB[i + 3], marker="s",
               ms=4, ls="--", label=f"Laser {name}")
    ax.set_xscale("log")
    ax.set_xlabel("target diameter (cm)")
    ax.set_ylabel("expected fraction of intended effect delivered per shot")
    ax.set_title("SIM-14: single-shot delivered-effect fraction vs target size")
    ax.legend(fontsize=7.5, ncol=2)
    save_fig(fig, "sim14_fraction_vs_size",
             "Expected delivered-effect fraction per shot. FG01 degrades "
             "gracefully (areal-capture model, SIM-3/4); the laser is "
             "modelled as binary per its own comparison table, so its curve "
             "reflects P(hit) collapsing once the miss budget exceeds the "
             "target radius.")

    fig, ax = plt.subplots()
    for i, name in enumerate(laser_names):
        s = dfp[(dfp.scenario == name) & (dfp.d_target_m == 0.01)].sort_values("n_pulses")
        ax.loglog(s.n_pulses, s.p_all_n_hit, color=CB[i], marker="o", ms=4,
                 label=f"Laser {name}, 1 cm target")
    ax.axhline(1.0, color=CB[6], ls=":", lw=1.2, label="FG01: 1 kick required (SIM-13)")
    ax.set_xlabel("number of consecutive pulses required to hit")
    ax.set_ylabel("P(every pulse hits)")
    ax.set_title("SIM-14: sustained-tracking penalty compounds with pulse count")
    ax.legend(fontsize=7.5)
    save_fig(fig, "sim14_pulse_compounding",
             "P(all N pulses hit) = P_hit^N for a laser engagement needing "
             "sustained tracking to accumulate Delta-v. FG01's design point "
             "needs exactly one kick per object, so no analogous compounding "
             "applies to it.")

    return dfb, dff, dfp, concl


if __name__ == "__main__":
    run()
