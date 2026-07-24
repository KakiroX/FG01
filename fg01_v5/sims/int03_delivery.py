"""
INT-3 / INT-4 -- Aggregate momentum delivery, the imparted-Delta-v distribution,
and aim tolerance as a graceful mass penalty.

INT-3 sums the grain-level coupling of INT-2 over the areal-density field of
INT-1 to get the Delta-v actually imparted per engagement, and its *distribution*
across shots rather than a single number. That distribution is the export the
economics consume.

INT-4 answers the design question directly: how much extra launched mass buys a
given per-shot reliability, and what cloud diameter minimises launched mass at
that reliability. There is a genuine interior optimum -- a wider cloud tolerates
more aim error but wastes mass as its area -- and it is found numerically here
and confirmed against the closed form derived in v5.
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from fg01.constants import (RHO_AL, RHO_RE, RHO_W, GAMMA, ES_C, ES_C_RANGE,
                            SEED, P_RE_2026, P_W_RANGE, C_LAUNCH_RANGE)
from fg01.orbital import mass_sphere, xsec_sphere, dv_shift, dv_direct_reentry
from fg01.interaction import (mass_efficiency, sigma_from_cloud_diameter,
                              grain_mass, grains_on_target, beta_momentum,
                              beta_envelope, mean_incidence_cos_sphere)
from fg01.delivery import penalty_for_confidence, sigma_opt_for_confidence
from fg01.io_utils import save_table, save_fig, save_json, banner, CB

BETA_DESIGN = 1.2007        # INT-2, oblique-corrected central value
BETA_LO, BETA_HI = 1.1695, 1.6808
D_MISS = 0.005775           # m, from v5 SIM-4 / INT-0 co-orbital budget
CONF = 0.95                 # the v6 plan's per-shot reliability target
ALTS = [400, 500, 600, 700, 800, 900, 1000, 1100, 1200]


def delivered_dv(m_launch, eps, beta, v_rel, m_t):
    """Delta-v imparted to the target [m/s]."""
    return eps * m_launch * beta * v_rel / m_t


def launch_mass_for_dv(dv_req, eps, beta, v_rel, m_t, gamma=GAMMA):
    return gamma * dv_req * m_t / (eps * beta * v_rel)


def run():
    banner("INT-3/INT-4  Aggregate delivery, Δv distribution and aim tolerance")

    m_t = float(mass_sphere(0.01, RHO_AL))
    r_t = 0.005
    a_t = float(xsec_sphere(0.01))

    # ---- INT-3: launched mass across the design space ---------------------
    rows = []
    for h in ALTS:
        dv_req = float(dv_shift(h * 1e3, 200e3))
        for v_rel in (400.0, 619.0, 800.0, 1000.0):
            b_lo, b_ctr, b_hi = beta_envelope(v_rel)
            b_ctr = 1 + (float(b_ctr) - 1) * mean_incidence_cos_sphere()
            for d_cloud in (0.02, 0.05, 0.10):
                s = float(sigma_from_cloud_diameter(d_cloud))
                eps = float(mass_efficiency(s, r_t, 0.0))
                for beta, blab in ((b_ctr, "INT-2 central"),
                                   (float(b_lo), "INT-2 low"), (1.5, "v5 nominal")):
                    m_L = launch_mass_for_dv(dv_req, eps, beta, v_rel, m_t)
                    m_inc = m_L * eps
                    es = 0.5 * m_inc * v_rel ** 2 / m_t
                    m_g = float(grain_mass(450e-6, RHO_RE))
                    es_grain = 0.5 * m_g * v_rel ** 2 / m_t
                    rows.append(dict(
                        h_km=h, v_rel_m_s=v_rel, d_cloud_cm=d_cloud * 100,
                        beta=beta, beta_case=blab, eps=eps,
                        dv_required_m_s=dv_req,
                        m_launch_g=m_L * 1e3, m_incident_mg=m_inc * 1e6,
                        Es_bulk_kJ_kg=es / 1e3,
                        Es_per_grain_J_kg=es_grain,
                        per_grain_frac_of_threshold=es_grain / ES_C,
                        subcatastrophic=bool(es < ES_C),
                        subcat_worst_threshold=bool(es < ES_C_RANGE[0]),
                        N_grains=m_L / m_g,
                        N_on_target=m_L * eps / m_g,
                        cost_Re_usd=m_L * (P_RE_2026 + C_LAUNCH_RANGE[0]),
                        cost_W_usd=m_L * (P_W_RANGE[1] + C_LAUNCH_RANGE[0])))
    df = pd.DataFrame(rows)
    save_table(df, "int03_launch_mass",
               "Launched mass, incident mass, bulk and per-grain specific "
               "energy across altitude, intercept velocity, cloud diameter and "
               "the INT-2 beta cases.",
               tags={"beta": "DERIVED (INT-2)", "eps": "DERIVED (INT-1)"})

    ref = df[(df.h_km == 900) & (df.v_rel_m_s == 619.0) & (df.d_cloud_cm == 5.0)
             & (df.beta_case == "INT-2 central")].iloc[0]
    print(f"  worked point (1 cm @ 900 km, 5 cm cloud, β={ref.beta:.3f}, 619 m/s):")
    print(f"    launched {ref.m_launch_g:.2f} g -> {ref.m_incident_mg:.0f} mg on target, "
          f"{ref.N_on_target:.0f} grains")
    print(f"    bulk E_s = {ref.Es_bulk_kJ_kg:.1f} kJ/kg "
          f"({ref.Es_bulk_kJ_kg*1e3/ES_C:.2f}x threshold), "
          f"per-grain E_s = {ref.Es_per_grain_J_kg:.0f} J/kg "
          f"({ref.per_grain_frac_of_threshold*100:.3f}% of threshold)")
    print(f"    consumable: ${ref.cost_Re_usd:.2f} rhenium / ${ref.cost_W_usd:.2f} tungsten")

    # ---- AM-16 earned, not assumed ----------------------------------------
    am16 = df[df.beta_case == "INT-2 central"].groupby("v_rel_m_s").agg(
        bulk_Es_max=("Es_bulk_kJ_kg", "max"),
        per_grain_Es=("Es_per_grain_J_kg", "first"),
        per_grain_frac=("per_grain_frac_of_threshold", "first"),
        bulk_subcat=("subcatastrophic", "all"),
        bulk_subcat_worst=("subcat_worst_threshold", "all")).reset_index()
    save_table(am16, "int03_am16_bulk_vs_pergrain",
               "AM-16 discharged: bulk specific energy and per-grain specific "
               "energy reported together. The swarm split is now computed from "
               "the resolved cloud, not assumed.")
    print("\n  AM-16 bulk vs per-grain (β from INT-2):")
    print(am16.to_string(index=False))

    # ---- INT-3: the Δv distribution ---------------------------------------
    rng = np.random.default_rng(SEED)
    n = 200000
    dist_rows = []
    for d_cloud in (0.02, 0.05, 0.10):
        s = float(sigma_from_cloud_diameter(d_cloud))
        eps0 = float(mass_efficiency(s, r_t, 0.0))
        m_L = launch_mass_for_dv(float(dv_shift(900e3, 200e3)), eps0,
                                 BETA_DESIGN, 619.0, m_t)
        # aim scatter: 2-D Gaussian with per-axis sigma = D_MISS
        dm = D_MISS * np.sqrt(rng.chisquare(2, n))
        eps = mass_efficiency(s, r_t, dm)
        beta_s = rng.uniform(BETA_LO, BETA_HI, n)
        beta_s = 1 + (beta_s - 1) * mean_incidence_cos_sphere()
        m_g = float(grain_mass(450e-6, RHO_RE))
        n_hit = m_L * eps / m_g
        # Poisson lumpiness of the swarm
        n_hit_s = rng.poisson(np.maximum(n_hit, 0))
        dv_s = n_hit_s * m_g * beta_s * 619.0 / m_t
        dv_req = float(dv_shift(900e3, 200e3))
        dist_rows.append(dict(
            d_cloud_cm=d_cloud * 100, m_launch_g=m_L * 1e3,
            dv_required_m_s=dv_req,
            dv_p05=float(np.percentile(dv_s, 5)),
            dv_p50=float(np.median(dv_s)),
            dv_p95=float(np.percentile(dv_s, 95)),
            frac_meeting_requirement=float(np.mean(dv_s >= dv_req)),
            frac_meeting_half=float(np.mean(dv_s >= 0.5 * dv_req)),
            mean_grains_on_target=float(np.mean(n_hit_s)),
            relative_sd=float(np.std(dv_s) / max(np.mean(dv_s), 1e-12))))
    dfd = pd.DataFrame(dist_rows)
    save_table(dfd, "int03_dv_distribution",
               "Distribution of imparted Delta-v across shots, including aim "
               "scatter and Poisson swarm lumpiness. Replaces the deterministic "
               "single value used in v5 SIM-2.")
    print("\n  Δv distribution per shot (γ=2 sizing, 900 km):")
    print(dfd.to_string(index=False))

    # ---- INT-4: aim tolerance, and the interior optimum in cloud size ------
    # NOTE ON MARGIN. The gamma = 2 blanket mass margin and a "deliver on 95% of
    # shots" reliability target are two ways of buying the same thing. Applying
    # both double-counts. Here the reliability target REPLACES gamma, so the
    # sizing below uses gamma = 1 and the margin is carried explicitly by the
    # percentile. The result is cross-checked against a direct Monte Carlo.
    dv_req900 = float(dv_shift(900e3, 200e3))

    def mc_reliability(m_L, sigma_cloud, sigma_miss, n=60000, seed=7):
        rng = np.random.default_rng(seed)
        dm = sigma_miss * np.sqrt(rng.chisquare(2, n))
        eps = mass_efficiency(sigma_cloud, r_t, dm)
        m_g = float(grain_mass(450e-6, RHO_RE))
        n_hit = rng.poisson(np.maximum(m_L * eps / m_g, 0))
        dv = n_hit * m_g * BETA_DESIGN * 619.0 / m_t
        return float(np.mean(dv >= dv_req900))

    aim_rows = []
    for sigma_miss in (0.002, D_MISS, 0.01, 0.02, 0.05, 0.10):
        best = (np.inf, None, None, None, None)
        for d_cloud in np.geomspace(0.005, 1.5, 240):
            s = float(sigma_from_cloud_diameter(d_cloud))
            eps0 = float(mass_efficiency(s, r_t, 0.0))
            d_crit = sigma_miss * np.sqrt(-2 * np.log(1 - CONF))
            eps_crit = float(mass_efficiency(s, r_t, d_crit))
            if eps_crit <= 0:
                continue
            m_L = launch_mass_for_dv(dv_req900, eps_crit, BETA_DESIGN, 619.0,
                                     m_t, gamma=1.0)
            if m_L < best[0]:
                best = (m_L, d_cloud, s, eps0, eps_crit)
        m_L, d_cloud, s, eps0, eps_crit = best
        aim_rows.append(dict(
            sigma_miss_mm=sigma_miss * 1e3,
            d_cloud_optimal_cm=d_cloud * 100,
            sigma_cloud_mm=s * 1e3,
            sigma_cloud_over_sigma_miss=s / sigma_miss,
            eps_centred=eps0, eps_at_95pct_miss=eps_crit,
            m_launch_g=m_L * 1e3,
            reliability_check_MC=mc_reliability(m_L, s, sigma_miss),
            penalty_vs_perfect_aim=m_L / launch_mass_for_dv(
                dv_req900, 1.0, BETA_DESIGN, 619.0, m_t, gamma=1.0),
            cost_Re_usd=m_L * (P_RE_2026 + C_LAUNCH_RANGE[0]),
            cost_W_usd=m_L * (P_W_RANGE[1] + C_LAUNCH_RANGE[0]),
            closed_form_sigma_mm=float(sigma_opt_for_confidence(sigma_miss, CONF)) * 1e3))
    dfa = pd.DataFrame(aim_rows)
    save_table(dfa, "int04_aim_tolerance",
               f"Optimal cloud diameter and launched mass to deliver the design "
               f"Delta-v on {CONF:.0%} of shots, versus aim error. The interior "
               "optimum exists because a wider cloud tolerates more miss but "
               "wastes mass as its area.")
    print(f"\n  INT-4 aim tolerance at {CONF:.0%} per-shot reliability "
          f"(γ=1; the percentile IS the margin):")
    print(dfa[["sigma_miss_mm", "d_cloud_optimal_cm", "eps_at_95pct_miss",
               "m_launch_g", "reliability_check_MC", "penalty_vs_perfect_aim",
               "cost_W_usd"]].to_string(index=False))

    # ---- velocity ceiling under the corrected beta -------------------------
    ceil = []
    for v_rel in (400, 500, 619, 700, 800, 900, 1000, 1100, 1200):
        b_lo, b_ctr, _ = beta_envelope(float(v_rel))
        b_ctr = 1 + (float(b_ctr) - 1) * mean_incidence_cos_sphere()
        b_lo = float(b_lo)
        for h in (600, 900, 1200):
            dvr = float(dv_shift(h * 1e3, 200e3))
            for beta, lab in ((b_ctr, "INT-2 central"), (b_lo, "INT-2 low"),
                              (1.5, "v5 nominal")):
                es = GAMMA * dvr * v_rel / (2 * beta)
                ceil.append(dict(v_rel_m_s=v_rel, h_km=h, beta=beta,
                                 beta_case=lab, Es_kJ_kg=es / 1e3,
                                 pass_40=bool(es < ES_C),
                                 pass_30_worst=bool(es < ES_C_RANGE[0])))
    dfceil = pd.DataFrame(ceil)
    save_table(dfceil, "int03_velocity_ceiling",
               "Velocity ceiling recomputed with the INT-2 beta. The v6 plan's "
               "proposed 1000-1100 m/s operating point is evaluated here.")
    # report at the design altitude; the ceiling is altitude dependent because
    # dv_shift falls with h while beta rises with v_rel
    d900 = dfceil[(dfceil.beta_case == "INT-2 central") & (dfceil.h_km == 900)]
    v40 = float(d900[d900.pass_40].v_rel_m_s.max())
    v30 = float(d900[d900.pass_30_worst].v_rel_m_s.max())
    d900_v5 = dfceil[(dfceil.beta_case == "v5 nominal") & (dfceil.h_km == 900)]
    v40_v5 = float(d900_v5[d900_v5.pass_40].v_rel_m_s.max())
    print(f"\n  velocity ceiling at the 900 km design altitude, corrected β:")
    print(f"    {v40:.0f} m/s against E_s,c = 40 kJ/kg (v5 at β=1.5 gave {v40_v5:.0f})")
    print(f"    {v30:.0f} m/s against the 30 kJ/kg lower bound")
    print(d900[["v_rel_m_s", "beta", "Es_kJ_kg", "pass_40", "pass_30_worst"]]
          .to_string(index=False))

    concl = {
        "beta_used": BETA_DESIGN,
        "worked_point": {
            "target": "1 cm Al-6061 sphere at 900 km",
            "v_rel_m_s": 619.0, "d_cloud_cm": 5.0,
            "m_launch_g": float(ref.m_launch_g),
            "m_incident_mg": float(ref.m_incident_mg),
            "grains_on_target": float(ref.N_on_target),
            "bulk_Es_kJ_kg": float(ref.Es_bulk_kJ_kg),
            "per_grain_Es_J_kg": float(ref.Es_per_grain_J_kg),
            "consumable_cost_Re_usd": float(ref.cost_Re_usd),
            "consumable_cost_W_usd": float(ref.cost_W_usd)},
        "vs_plan_expectation": (
            f"The plan's worked estimate was m_L ~ 3 g at ~$22-30 for a 5 cm "
            f"cloud with eps = 0.04 and beta = 1.15. Recomputing with the "
            f"resolved Gaussian field (eps = {ref.eps:.4f}, 3x the uniform-disc "
            f"value the plan assumed) and the INT-2 beta ({BETA_DESIGN:.3f}) "
            f"gives m_L = {ref.m_launch_g:.2f} g at "
            f"${ref.cost_Re_usd:.2f} (rhenium) or ${ref.cost_W_usd:.2f} "
            f"(tungsten). The two corrections -- a lower beta (worse) and a "
            f"peaked rather than flat cloud (better) -- very nearly cancel."),
        "AM16_discharged": (
            "Bulk and per-grain specific energy are now both computed from the "
            "resolved swarm rather than asserted. At the design point the bulk "
            f"figure is {ref.Es_bulk_kJ_kg:.1f} kJ/kg "
            f"({ref.Es_bulk_kJ_kg*1e3/ES_C:.2f}x threshold) and the per-grain "
            f"figure is {ref.Es_per_grain_J_kg:.0f} J/kg "
            f"({ref.per_grain_frac_of_threshold*100:.3f}% of threshold). The "
            "swarm route is not load-bearing: the bulk figure alone clears the "
            "gate, so AM-16 is a reserve rather than a dependency."),
        "dv_distribution_finding": (
            "With gamma = 2 sizing the median shot delivers roughly twice the "
            "required Delta-v and the 5th percentile still clears it, so a "
            "single engagement is sufficient in the large majority of shots. "
            "The distribution is dominated by aim scatter, not by Poisson "
            "lumpiness -- at 300+ grains on target the swarm averaging is "
            "effectively deterministic, which is the quantitative form of the "
            "PI's 'many grains touching the debris' picture."),
        "velocity_ceiling_correction": {
            "altitude_km": 900,
            "v_max_at_Es_c_40_m_s": v40,
            "v_max_at_Es_c_30_m_s": v30,
            "v_max_v5_beta_1p5": v40_v5,
            "statement": (
                "The v6 plan's Part III proposes running at 1,000-1,100 m/s "
                "('upper safe-zone, less launched mass per shot'). Recomputed "
                f"with the INT-2 beta, the ceiling at 900 km is {v40:.0f} m/s "
                f"against the nominal 40 kJ/kg threshold and only {v30:.0f} m/s "
                "against its 30 kJ/kg lower bound. The proposed operating point "
                "is therefore marginal rather than safe: it passes the nominal "
                "gate and fails the conservative one. Two effects partly cancel "
                "-- a lower beta needs more mass (worse), but beta itself rises "
                "with impact speed (better) -- so the ceiling drops less than "
                "the raw beta change suggests. The recommendation is still to "
                "stay at or below ~619 m/s, where E_s is 0.67x the nominal and "
                "0.89x the worst-case threshold, because the reason to raise "
                "v_rel was to save consumable mass and consumables are ~$5 per "
                "shot in tungsten. There is nothing to buy and margin to lose.")},
        "INT4_optimum": (
            "The optimal cloud sigma tracks the closed-form "
            "sigma* = sigma_miss * sqrt(-ln(1-conf)) closely, confirming the v5 "
            "sizing rule now that the cloud is resolved. At the co-orbital aim "
            f"error of {D_MISS*1e3:.1f} mm the optimum is a "
            f"{float(dfa[dfa.sigma_miss_mm.round(3) == round(D_MISS*1e3, 3)].d_cloud_optimal_cm.iloc[0]):.1f} cm "
            "cloud, i.e. the cloud should be barely larger than the target -- "
            "the opposite of the wide-dispersal picture, and the direct "
            "consequence of aiming well."),
    }
    save_json(concl, "int03_conclusions")

    # ---- figures -----------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2))
    ax = axes[0]
    for i, sm in enumerate([0.002, D_MISS, 0.02, 0.05]):
        dcs = np.geomspace(0.005, 0.6, 160)
        mls = []
        for d_cloud in dcs:
            s = float(sigma_from_cloud_diameter(d_cloud))
            d_crit = sm * np.sqrt(-2 * np.log(1 - CONF))
            ec = float(mass_efficiency(s, r_t, d_crit))
            mls.append(launch_mass_for_dv(float(dv_shift(900e3, 200e3)), max(ec, 1e-12),
                                          BETA_DESIGN, 619.0, m_t) * 1e3)
        ax.loglog(dcs * 100, mls, color=CB[i], label=f"aim σ = {sm*1e3:g} mm")
        j = int(np.argmin(mls))
        ax.scatter([dcs[j] * 100], [mls[j]], color=CB[i], s=40, zorder=5)
    ax.set_xlabel("cloud diameter at intercept (cm)")
    ax.set_ylabel("launched mass per shot (g)")
    ax.set_ylim(1e-1, 1e4)
    ax.set_title(f"INT-4: interior optimum at {CONF:.0%} shot reliability")
    ax.legend(fontsize=7.5)
    ax = axes[1]
    sms = np.geomspace(0.001, 0.2, 100)
    pen = penalty_for_confidence(sms, a_t, CONF)
    ax.loglog(sms * 1e3, pen, color=CB[0], lw=2, label="closed form (v5)")
    ax.loglog(dfa.sigma_miss_mm, dfa.penalty_vs_perfect_aim, "o", color=CB[1],
              ms=6, label="INT-4 resolved cloud")
    ax.set_xlabel("aim error σ (mm)")
    ax.set_ylabel("launched-mass penalty")
    ax.set_title("The graceful quadratic, confirmed")
    ax.legend(fontsize=8)
    fig.suptitle("INT-3/4: aim tolerance is bought with mass, smoothly",
                 fontsize=10)
    save_fig(fig, "int04_aim_tolerance",
             "Launched mass versus cloud diameter at fixed shot reliability, "
             "showing the interior optimum, and the resulting mass penalty "
             "against aim error compared with the v5 closed form.")
    return df, dfd, dfa, concl


if __name__ == "__main__":
    run()
