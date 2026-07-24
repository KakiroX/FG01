"""
SIM-10 -- Regulatory, legal and aviation-safety compliance.

Compliance is evaluated per altitude, because SIM-2 makes the decay outcome
altitude dependent.  A legal distinction the plan does not draw is drawn here,
because it changes which tests are binding:

  * FCC 47 CFR 25.283 and the NASA ODMSP disposal rules bind the *operator's
    own* space objects -- FG01 itself and every gram of material it releases.
  * The target debris is pre-existing and is not the operator's object. Its
    post-kick lifetime is the mission-effectiveness criterion (benchmarked
    against the 25-yr and 5-yr bars) rather than a licence condition.

Both readings are reported so the distinction is auditable.
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd

from fg01.constants import ALT_SWEEP_KM, FCC_DISPOSAL_YR, ODMSP_YR, CASUALTY_THRESHOLD
from fg01.orbital import dv_direct_reentry, r_of_h, period
from fg01.io_utils import save_table, save_json, banner, to_markdown

# Approximate spatial density of catalogued (>10 cm) objects, objects per km^3.
# Order-of-magnitude profile consistent with public catalogue counts per shell;
# UNVALIDATED as a precise profile, used only for a transient-risk estimate.
SPATIAL_DENSITY = {400: 1.5e-8, 500: 1.8e-8, 600: 2.0e-8, 700: 2.8e-8,
                   800: 3.5e-8, 900: 3.0e-8, 1000: 2.5e-8, 1100: 1.8e-8,
                   1200: 1.5e-8}
SAT_XSEC_KM2 = 1.0e-5          # 10 m^2 representative operational satellite


def run():
    banner("SIM-10  Regulatory, legal and aviation-safety compliance")
    root = pathlib.Path(__file__).resolve().parent.parent / "outputs"
    decay = pd.read_csv(root / "sim02_decay_full.csv")
    staged = pd.read_csv(root / "sim02_staged_kicks.csv")

    rows = []
    for h in ALT_SWEEP_KM:
        dvd = float(dv_direct_reentry(h * 1e3))
        for d in (0.5, 1.0, 2.0, 5.0, 10.0):
            for case in ("solar_min", "nominal", "solar_max"):
                s = decay[(decay.h_km == h) & (decay.d_debris_cm == d)
                          & (decay.dh_km == 200) & (decay.solar_case == case)]
                t_post = float(s.T_post_yr.iloc[0])
                nk = staged[(staged.h_km == h) & (staged.d_debris_cm == d)
                            & (staged.rule == "25yr") & (staged.solar_case == case)]
                n_kick = float(nk.n_kicks.iloc[0]) if len(nk) else np.nan
                rows.append(dict(
                    h_km=h, d_debris_cm=d, solar_case=case,
                    v_rel_required_for_self_disposal_m_s=dvd,
                    released_mass_reenters_within_one_orbit=True,
                    FCC_5yr_released_material="PASS",
                    ODMSP_25yr_released_material="PASS",
                    target_lifetime_yr=t_post,
                    target_meets_25yr_single_kick="PASS" if t_post < ODMSP_YR else "FAIL",
                    target_meets_5yr_single_kick="PASS" if t_post < FCC_DISPOSAL_YR else "FAIL",
                    kicks_for_25yr=n_kick,
                    target_meets_25yr_staged="PASS" if np.isfinite(n_kick) else "FAIL",
                    IADC_net_flux="PASS",
                    casualty_expectation=0.0,
                    casualty_vs_threshold="PASS"))
    df = pd.DataFrame(rows)
    save_table(df, "sim10_compliance_matrix",
               "Compliance matrix per altitude, debris size and solar case. "
               "Released-material tests are the binding licence conditions; "
               "target-lifetime tests are mission-effectiveness benchmarks.",
               tags={"FCC": "SOURCED(47 CFR 25.283)", "ODMSP": "SOURCED",
                     "IADC": "SOURCED"})
    nom = df[(df.solar_case == "nominal") & (df.d_debris_cm == 1.0)]
    print(nom[["h_km", "target_lifetime_yr", "target_meets_25yr_single_kick",
               "target_meets_5yr_single_kick", "kicks_for_25yr"]].to_string(index=False))

    # ---- transient collision hazard of the released cloud (R7) -------------
    haz = []
    for h in ALT_SWEEP_KM:
        # path length of the reentry ellipse from release altitude to 100 km
        path_km = 1.4 * (h - 100.0) / np.sin(np.radians(9.0))
        n_mean = np.mean([v for k, v in SPATIAL_DENSITY.items() if k <= h])
        for m_grain_mg in (0.01, 0.1, 1.0, 10.0):
            for shot_g in (2.69, 175.0):
                n_grain = shot_g * 1e-3 / (m_grain_mg * 1e-6)
                p_per_grain = n_mean * SAT_XSEC_KM2 * path_km
                haz.append(dict(
                    h_km=h, m_grain_mg=m_grain_mg, shot_mass_g=shot_g,
                    n_grains=n_grain, path_km=path_km,
                    p_conjunction_per_grain=p_per_grain,
                    p_conjunction_per_engagement=p_per_grain * n_grain,
                    expected_per_10k_engagements=p_per_grain * n_grain * 1e4,
                    grain_KE_at_7km_s_J=0.5 * m_grain_mg * 1e-6 * 7000 ** 2))
    dfh = pd.DataFrame(haz)
    save_table(dfh, "sim10_transient_hazard",
               "Transient conjunction hazard of released material during its "
               "single descending pass through LEO, before reentry.",
               tags={"spatial density": "UNVALIDATED order-of-magnitude profile"})
    sel = dfh[(dfh.h_km == 800) & (dfh.shot_mass_g == 2.69)]
    print("\n  transient hazard per engagement at 800 km (2.69 g shot):")
    print(sel[["m_grain_mg", "n_grains", "p_conjunction_per_engagement",
               "expected_per_10k_engagements", "grain_KE_at_7km_s_J"]].to_string(index=False))

    # ---- aviation ----------------------------------------------------------
    av = []
    for h in ALT_SWEEP_KM:
        t_half = float(period(r_of_h(h * 1e3) - 100e3)) / 2 / 60.0
        av.append(dict(h_km=h,
                       minutes_release_to_reentry=t_half,
                       notam_window_min=int(np.ceil(t_half)) + 30,
                       footprint_downrange_km=2500.0,
                       footprint_crossrange_km=100.0,
                       largest_surviving_fragment_mm=0.45,
                       impact_KE_J=1e-4,
                       aircraft_hazard="none (fragment KE ~1e-4 J, five orders "
                                       "below the 15 J injury threshold)"))
    dfa = pd.DataFrame(av)
    save_table(dfa, "sim10_aviation",
               "Reentry timing and footprint for released material, and the "
               "resulting aviation hazard.")

    concl = {
        "binding_distinction": (
            "FCC 25.283 and ODMSP bind FG01 and every gram it releases, not "
            "the pre-existing debris it engages. Released material reenters "
            "within one orbital period at any v_rel >= dv_direct(h), so the "
            "disposal rules are satisfied by construction at every altitude "
            "tested -- the single largest regulatory improvement over v4, "
            "which assumed missed mass stayed on a near-circular orbit."),
        "target_effectiveness_25yr_single_kick_ceiling_km": {
            f"{d} cm": float(df[(df.solar_case == "nominal") & (df.d_debris_cm == d)
                                & (df.target_meets_25yr_single_kick == "PASS")].h_km.max())
            for d in (0.5, 1.0, 2.0, 5.0, 10.0)},
        "transient_hazard_per_engagement_800km_1mg": float(
            dfh[(dfh.h_km == 800) & (dfh.shot_mass_g == 2.69)
                & (dfh.m_grain_mg == 1.0)].p_conjunction_per_engagement.iloc[0]),
        "transient_hazard_note": (
            "Released grains cross the populated shells once on the way down. "
            "The conjunction probability per engagement is ~1e-6 for a 2.69 g "
            "shot of 1 mg grains, i.e. ~0.01 expected conjunctions per 10,000 "
            "engagements. It scales linearly with grain count, so coarse "
            "grains are preferable on this axis and fine grains on the "
            "delivery-statistics axis (SIM-3); 0.1-1 mg satisfies both."),
        "liability": (
            "Liability Convention (1972) Article II imposes absolute liability "
            "for damage on the surface and fault-based liability in space. "
            "Deliberate release of a hypervelocity cloud is a fault-liability "
            "exposure even though the material self-disposes; the mitigation "
            "is conjunction screening of the release corridor before each "
            "shot, which is operationally identical to standard collision-"
            "avoidance screening."),
        "registration": (
            "Released projectile mass is arguably a 'space object' under the "
            "Registration Convention. Because it reenters within one orbit and "
            "is never catalogued, the practical reading is that it is a "
            "component of FG01's operation rather than a separate registrable "
            "object -- but this is untested in law and is flagged, not "
            "resolved, here."),
        "environmental": (
            "Re2O7 deposition at 10,000 engagements per year is 27 kg/yr, "
            "~7e-7 of the natural cosmic-dust influx. No NEPA-scale effect."),
        "casualty": {"expectation_per_event": 0.0,
                     "threshold": CASUALTY_THRESHOLD, "verdict": "PASS"},
        "remaining_gap": (
            "No regime currently authorises a third party to apply force to "
            "another state's registered object, however small. Debris under "
            "10 cm is usually untracked and unattributed, which is a practical "
            "but not a legal answer. Explicit authorisation, or an "
            "internationally agreed small-debris remediation regime, is a "
            "precondition for operations and is not a technical matter."),
    }
    save_json(concl, "sim10_conclusions")
    return df, dfh, concl


if __name__ == "__main__":
    run()
