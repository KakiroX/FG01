"""
ECO-2 (SIM-15) -- Cluster availability in the real debris population.

The make-or-break number: how many co-orbital >=1 cm targets can one platform
actually reach at low Delta-v?

The method is to synthesise each real breakup cloud from its parent orbit using
the NASA SBM ejection-velocity distribution, age it under differential J2 nodal
regression to the present day, and count how many fragments fall inside the
plane-angle cone INT-0 hands over (1.5 deg for a <=10% energy penalty) and
inside an altitude band the platform can reach cheaply.

The result is then cross-checked against the ESA MASTER total, and the
altitude-vs-decay-benefit tension of v5 R2 is applied: a cluster is only useful
if a 200 km kick there actually changes the compliance status of its members.
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from fg01.constants import SEED, R_E, MU, ODMSP_YR, RHO_AL
from fg01.orbital import area_to_mass_sphere, elements_from_apo_peri, v_circ
from fg01.population import (CLUSTERS, synth_cloud, age_cloud, engageable_mask,
                             scale_to_population, MASTER_TOTAL_GE_1CM,
                             MASTER_TOTAL_GE_10CM, MASTER_TRACKED)
from fg01.relmotion import nodal_rate, plane_angle
from fg01.decay import lifetime, YEAR
from fg01.io_utils import save_table, save_fig, save_json, banner, CB

THETA_MAX_DEG = 1.5          # from INT-0 (<=10% energy penalty)
N_SYNTH = 60000
NOW = 2026.6


def natural_lifetime(h_km, d_cm=1.0, inc_deg=98.8):
    aom = float(area_to_mass_sphere(d_cm / 100, RHO_AL))
    t, _ = lifetime(R_E + h_km * 1e3, 0.0, aom, inc_deg=inc_deg,
                    case="nominal", t_max=300 * YEAR)
    return t / YEAR


def post_kick_lifetime(h_km, d_cm=1.0, dh_km=200.0, inc_deg=98.8):
    aom = float(area_to_mass_sphere(d_cm / 100, RHO_AL))
    a, e = elements_from_apo_peri(h_km * 1e3, (h_km - dh_km) * 1e3)
    t, _ = lifetime(a, e, aom, inc_deg=inc_deg, case="nominal", t_max=300 * YEAR)
    return t / YEAR


def run():
    banner("ECO-2  Cluster availability in the real debris population")
    rng = np.random.default_rng(SEED + 15)
    th_max = np.radians(THETA_MAX_DEG)

    rows = []
    spread_rows = []
    for name, c in CLUSTERS.items():
        cloud = synth_cloud(name, N_SYNTH, rng)
        age = NOW - c["year"]
        aged = age_cloud(cloud, age)
        w = scale_to_population(cloud, c["n_ge_1cm"])

        # how far has the cloud spread in RAAN?
        raan_spread = float(np.ptp(np.degrees(aged["raan"])))
        inc_spread = float(np.percentile(np.degrees(aged["inc"]), 95)
                           - np.percentile(np.degrees(aged["inc"]), 5))
        alt_spread = float(np.percentile((aged["a"] - R_E) / 1e3, 95)
                           - np.percentile((aged["a"] - R_E) / 1e3, 5))
        spread_rows.append(dict(
            cluster=name, age_yr=age, raan_spread_deg=raan_spread,
            inc_spread_5_95_deg=inc_spread, alt_spread_5_95_km=alt_spread,
            nodal_rate_deg_day=float(np.degrees(nodal_rate(
                R_E + c["h_km"] * 1e3, np.radians(c["inc_deg"]))) * 86400)))

        # place the platform at the cloud's modal plane and altitude, then
        # count engageable members; scan RAAN to find the best station
        a_plat = float(np.median(aged["a"]))
        inc_plat = float(np.median(aged["inc"]))
        best = (0, 0.0)
        for raan_p in np.linspace(0, 2 * np.pi, 72, endpoint=False):
            m, _ = engageable_mask(aged, a_plat, inc_plat, raan_p, th_max,
                                   dh_band_m=50e3)
            if m.sum() > best[0]:
                best = (int(m.sum()), raan_p)
        n_eng_instant = best[0] * w

        # widen the altitude band: the platform can change altitude cheaply
        for dh_band in (25e3, 50e3, 100e3, 200e3):
            m, _ = engageable_mask(aged, a_plat, inc_plat, best[1], th_max,
                                   dh_band_m=dh_band)
            rows.append(dict(
                cluster=name, h_km=c["h_km"], inc_deg=c["inc_deg"],
                age_yr=age, n_ge_1cm_modelled=c["n_ge_1cm"],
                n_tracked=c["n_tracked"],
                dh_band_km=dh_band / 1e3,
                theta_max_deg=THETA_MAX_DEG,
                n_engageable_instant=float(m.sum() * w),
                frac_of_cloud=float(m.sum() * w / c["n_ge_1cm"]),
                dv_band_m_s=float(abs(v_circ(c["h_km"] * 1e3)
                                      - v_circ(c["h_km"] * 1e3 + dh_band)) * 2)))
    df = pd.DataFrame(rows)
    dfs = pd.DataFrame(spread_rows)
    save_table(dfs, "eco02_cloud_spread",
               "How far each real breakup cloud has spread under differential "
               "J2 nodal regression since its event.",
               tags={"SBM dv": "SOURCED", "ages": "SOURCED"})
    save_table(df, "eco02_engageable",
               "Instantaneously engageable >=1 cm members per cluster, inside "
               "the INT-0 plane cone and an altitude band.",
               tags={"n_ge_1cm": "UNVALIDATED (modelled, not catalogued)"})
    print(dfs.to_string(index=False))
    print("\n  instantaneously engageable members (θ ≤ 1.5°):")
    print(df[df.dh_band_km == 50.0][
        ["cluster", "n_ge_1cm_modelled", "n_engageable_instant",
         "frac_of_cloud"]].to_string(index=False))

    # ---- the flux: members drifting INTO the cone over the mission ----------
    # differential nodal rates sweep fragments through the platform's plane
    flux_rows = []
    for name, c in CLUSTERS.items():
        cloud = synth_cloud(name, N_SYNTH, rng)
        aged = age_cloud(cloud, NOW - c["year"])
        w = scale_to_population(cloud, c["n_ge_1cm"])
        a_plat = float(np.median(aged["a"]))
        inc_plat = float(np.median(aged["inc"]))
        rate_plat = float(nodal_rate(a_plat, inc_plat))
        for dh_band in (50e3, 100e3):
            in_alt = np.abs(aged["a"] - a_plat) <= dh_band
            # relative RAAN drift rate for the altitude-matched subset
            d_rate = nodal_rate(aged["a"][in_alt], aged["inc"][in_alt]) - rate_plat
            # a fragment enters the cone when its relative RAAN sweeps 2*theta
            # cumulative count over the mission = sum over fragments of the
            # number of cone crossings in the mission duration
            for mission_yr in (5.0, 10.0):
                sweep = np.abs(d_rate) * mission_yr * 365.25 * 86400.0
                crossings = np.clip(sweep / (2 * np.pi), 0, None)
                # each full 2pi sweep gives one pass through the cone
                n_pass = float(np.sum(np.minimum(crossings, 50)) * w)
                # plus those already inside
                m0, _ = engageable_mask(aged, a_plat, inc_plat,
                                        float(np.median(aged["raan"])), th_max,
                                        dh_band_m=dh_band)
                flux_rows.append(dict(
                    cluster=name, dh_band_km=dh_band / 1e3,
                    mission_yr=mission_yr,
                    n_in_cone_at_start=float(m0.sum() * w),
                    n_passing_through_cone=n_pass,
                    n_total_opportunities=float(m0.sum() * w) + n_pass,
                    median_rel_nodal_rate_deg_day=float(
                        np.degrees(np.median(np.abs(d_rate))) * 86400)))
    dff = pd.DataFrame(flux_rows)
    save_table(dff, "eco02_flux",
               "Engagement opportunities over a mission: members already in the "
               "plane cone plus those whose differential nodal drift sweeps "
               "them through it.")
    print("\n  engagement opportunities over a 10-yr mission:")
    print(dff[(dff.mission_yr == 10.0) & (dff.dh_band_km == 50.0)][
        ["cluster", "n_in_cone_at_start", "n_passing_through_cone",
         "n_total_opportunities", "median_rel_nodal_rate_deg_day"]].to_string(index=False))

    # ---- decay benefit: is the cluster even worth servicing? ---------------
    ben = []
    for name, c in CLUSTERS.items():
        t_nat = natural_lifetime(c["h_km"], 1.0, c["inc_deg"])
        t_post = post_kick_lifetime(c["h_km"], 1.0, 200.0, c["inc_deg"])
        ben.append(dict(
            cluster=name, h_km=c["h_km"],
            T_natural_yr=t_nat, T_post_kick_yr=t_post,
            speedup=t_nat / t_post if t_post > 0 else np.inf,
            natural_already_compliant=bool(t_nat < ODMSP_YR),
            post_kick_compliant=bool(t_post < ODMSP_YR),
            remediation_value=bool(t_nat >= ODMSP_YR and t_post < ODMSP_YR)))
    dfb = pd.DataFrame(ben)
    save_table(dfb, "eco02_decay_benefit",
               "Does a 200 km kick change the compliance status of a 1 cm "
               "member of each cluster? (v5 R2 tension applied per cluster.)")
    print("\n  decay benefit per cluster (1 cm members):")
    print(dfb.to_string(index=False))

    # ---- the viability-defining number -------------------------------------
    merged = dff[(dff.mission_yr == 10.0) & (dff.dh_band_km == 100.0)].merge(
        dfb[["cluster", "remediation_value", "T_natural_yr", "T_post_kick_yr"]],
        on="cluster")
    merged["useful_opportunities"] = np.where(
        merged.remediation_value, merged.n_total_opportunities, 0.0)
    save_table(merged, "eco02_viability_number",
               "THE viability-defining number: engagement opportunities that "
               "are both reachable and worth taking, per platform per decade.")
    print("\n  useful engagement opportunities per platform per decade:")
    print(merged[["cluster", "n_total_opportunities", "remediation_value",
                  "useful_opportunities"]].to_string(index=False))

    best_cluster = merged.loc[merged.useful_opportunities.idxmax()]
    n_best = float(best_cluster.useful_opportunities)

    concl = {
        "method": (
            "Each real cloud synthesised from its parent orbit with the NASA SBM "
            "ejection-velocity distribution, aged under differential J2 nodal "
            "regression to 2026.6, then filtered by the INT-0 plane cone "
            f"({THETA_MAX_DEG} deg) and an altitude band."),
        "population_context": {
            "MASTER_total_ge_1cm": MASTER_TOTAL_GE_1CM,
            "MASTER_total_ge_10cm": MASTER_TOTAL_GE_10CM,
            "tracked": MASTER_TRACKED,
            "note": ("Only ~3% of the >=1 cm population is tracked, so target "
                     "acquisition cannot rely on the public catalogue -- this is "
                     "v5 R10 restated as an architecture requirement, not a "
                     "footnote.")},
        "cloud_spread_finding": (
            "Every cloud studied has spread to a full 360 deg in RAAN. "
            "Differential nodal regression, driven by the spread in semi-major "
            "axis and inclination the breakup itself imparts, turns a compact "
            "cloud into a complete shell within a few years. The 'cluster' the "
            "v6 plan hoped to service as a co-orbital group no longer exists as "
            "a group: what survives is a shared *inclination band*, not a "
            "shared plane."),
        "engageable_instant": {
            r.cluster: float(r.n_engageable_instant)
            for _, r in df[df.dh_band_km == 50.0].iterrows()},
        "opportunities_per_decade": {
            r.cluster: float(r.n_total_opportunities)
            for _, r in merged.iterrows()},
        "best_cluster": str(best_cluster.cluster),
        "best_useful_opportunities_per_decade": n_best,
        "verdict_vs_thresholds": {
            "niche_threshold_1e4": bool(n_best >= 1e4),
            "laser_competitive_1e5": bool(n_best >= 1e5),
        },
        "R2_tension_per_cluster": (
            "Cosmos-1408 at 480 km is the densest low cluster but its 1 cm "
            "members decay unaided well inside 25 years, so engaging them has "
            "no remediation value. Fengyun-1C at 865 km has real remediation "
            "value and is the largest cloud, which is the combination the "
            "architecture needs -- but it sits at the top of the altitude band "
            "where a single 200 km kick still complies."),
    }
    save_json(concl, "eco02_conclusions")
    print(f"\n  BEST: {best_cluster.cluster} -> {n_best:,.0f} useful "
          f"opportunities per platform-decade")
    print(f"  niche threshold (1e4): {'MET' if n_best >= 1e4 else 'NOT MET'}; "
          f"laser-competitive (1e5): {'MET' if n_best >= 1e5 else 'NOT MET'}")

    # ---- figures -----------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2))
    ax = axes[0]
    cl = list(CLUSTERS)
    x = np.arange(len(cl))
    ax.bar(x - 0.2, [CLUSTERS[c]["n_ge_1cm"] for c in cl], 0.4, color=CB[0],
           label="modelled ≥1 cm members")
    ax.bar(x + 0.2, [float(merged[merged.cluster == c].n_total_opportunities.iloc[0])
                     for c in cl], 0.4, color=CB[1],
           label="reachable opportunities / decade")
    ax.axhline(1e4, color="k", ls="--", lw=1.3)
    ax.text(-0.4, 1.2e4, "niche threshold (10⁴)", fontsize=7.5)
    ax.set_yscale("log"); ax.set_xticks(x)
    ax.set_xticklabels(cl, fontsize=7.5, rotation=20, ha="right")
    ax.set_ylabel("objects")
    ax.set_title("Cloud size vs. what one platform can reach")
    ax.legend(fontsize=7.5)
    ax = axes[1]
    ax.bar(x, dfb.speedup, 0.5, color=[CB[2] if r else CB[3]
                                       for r in dfb.remediation_value])
    ax.set_xticks(x); ax.set_xticklabels(cl, fontsize=7.5, rotation=20, ha="right")
    ax.set_ylabel("lifetime reduction factor from a 200 km kick")
    ax.set_title("Green = the kick changes compliance status")
    fig.suptitle("ECO-2: cluster availability and whether servicing it helps",
                 fontsize=10)
    save_fig(fig, "eco02_clusters",
             "Modelled ≥1 cm membership of each major breakup cloud against the "
             "number one platform can actually reach in a decade, and whether a "
             "200 km kick changes the compliance status of members at that "
             "altitude.")
    return df, dff, dfb, concl


if __name__ == "__main__":
    run()
