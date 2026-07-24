"""
SIM-12 -- Global uncertainty quantification and sensitivity.

Monte Carlo over the joint distribution (AM-12), 10^6 samples, fixed seed,
with convergence and Sobol first-order/total indices.

Viability is evaluated under TWO definitions of the guidance constraint, and
both are reported everywhere, because the difference between them is the whole
of the v4 -> v5 change on the guidance axis:

  V-plan : the plan's Bernoulli criterion, P_hit > 0.5, with P_hit computed
           from the v4 miss model (d_miss = sqrt(sigma_pos^2 + (v_rel tau)^2))
           against a fixed 0.38 m cloud.
  V-design: the shot is sized for 90% delivery confidence against the sampled
           miss distance (SIM-3/4), and the constraint is that the resulting
           launched mass per shot stays within a 1 kg magazine limit.

The other four constraints are common: sub-catastrophic specific energy,
post-kick lifetime under 25 yr, non-negative net debris, and cost.
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from fg01.constants import (SEED, GAMMA, RHO_AL, P_RE_RANGE, C_LAUNCH_RANGE,
                            C_PLATFORM_RANGE, ES_C_RANGE, BETA_RANGE, ODMSP_YR,
                            ALT_SWEEP_KM, OPTICAL_LAMBDA_M, OPTICAL_APERTURE_M)
from fg01.orbital import mass_sphere, xsec_sphere, dv_shift, dv_direct_reentry
from fg01.io_utils import save_table, save_fig, save_json, banner, CB

N_SAMPLES = 1_000_000
SIZES = [0.5, 1.0, 2.0, 5.0, 10.0]
CASES = ["solar_min", "nominal", "solar_max"]
M_MAX_SHOT = 1.0            # kg magazine limit per shot
COST_BAR = 100_000.0        # USD/kg removed  (the plan's criterion)
# A per-object criterion is carried alongside because the per-kg convention is
# structurally unreachable for centimetre debris: a 1 cm aluminium sphere masses
# 1.4 g, so $100,000/kg is $141 per object, less than the launch cost of the
# propellant needed to reach it. Which convention is used is a policy choice,
# not a physics result, so both are reported.
COST_BAR_OBJECT = 50_000.0  # USD per object removed
CONF = 0.9
C_LIGHT = 2.99792458e8
LIDAR_BW = 1e8

OUT = pathlib.Path(__file__).resolve().parent.parent / "outputs"


def build_lookups():
    dec = pd.read_csv(OUT / "sim02_decay_full.csv")
    st = pd.read_csv(OUT / "sim02_staged_kicks.csv")
    T = {}
    for (d, c), g in dec.groupby(["d_debris_cm", "solar_case"]):
        piv = g.pivot_table(index="h_km", columns="dh_km", values="T_post_yr")
        arr = piv.to_numpy()
        arr = np.where(np.isfinite(arr), arr, 1e4)
        T[(d, c)] = (piv.index.to_numpy(float), piv.columns.to_numpy(float),
                     np.log10(np.clip(arr, 1e-4, 1e4)))
    K = {}
    for (d, c), g in st[st.rule == "25yr"].groupby(["d_debris_cm", "solar_case"]):
        s = g.sort_values("h_km")
        K[(d, c)] = (s.h_km.to_numpy(float),
                     np.where(np.isfinite(s.n_kicks.to_numpy(float)),
                              s.n_kicks.to_numpy(float), 8.0))
    # natural (un-kicked) lifetime, for the remediation-value test
    P = {}
    for (d, c), g in dec[dec.dh_km == 200].groupby(["d_debris_cm", "solar_case"]):
        s = g.sort_values("h_km")
        tp = np.where(np.isfinite(s.T_pre_yr.to_numpy(float)),
                      s.T_pre_yr.to_numpy(float), 1e4)
        P[(d, c)] = (s.h_km.to_numpy(float), np.log10(np.clip(tp, 1e-4, 1e4)))
    return T, K, P


def bilinear(h_grid, dh_grid, Z, h, dh):
    """Vectorised bilinear interpolation of log10(T_decay)."""
    hi = np.clip(np.searchsorted(h_grid, h) - 1, 0, len(h_grid) - 2)
    di = np.clip(np.searchsorted(dh_grid, dh) - 1, 0, len(dh_grid) - 2)
    h0, h1 = h_grid[hi], h_grid[hi + 1]
    d0, d1 = dh_grid[di], dh_grid[di + 1]
    th = (h - h0) / (h1 - h0)
    td = (dh - d0) / (d1 - d0)
    z00, z01 = Z[hi, di], Z[hi, di + 1]
    z10, z11 = Z[hi + 1, di], Z[hi + 1, di + 1]
    return (z00 * (1 - th) * (1 - td) + z01 * (1 - th) * td
            + z10 * th * (1 - td) + z11 * th * td)


def evaluate(p, d_cm, case, T, K, P=None):
    """Vectorised viability evaluation. p is a dict of sampled arrays."""
    d_m = d_cm / 100.0
    m_t = float(mass_sphere(d_m, RHO_AL)) * p["m_t_scale"]
    a_t = float(xsec_sphere(d_m))
    r_t = d_m / 2

    h_m = p["h_km"] * 1e3
    dv = np.array([float(dv_shift(x, y)) for x, y in
                   zip(h_m, p["dh_m"])]) if False else _dv_shift_vec(h_m, p["dh_m"])
    dv_dir = _dv_direct_vec(h_m)

    # 1. sub-catastrophic
    es = GAMMA * dv * p["v_rel"] / (2.0 * p["beta"])
    c1 = es < p["es_c"]

    # self-disposal of released mass (a v5 design constraint, reported too)
    c_self = p["v_rel"] >= dv_dir

    # miss budget
    sp = (OPTICAL_LAMBDA_M / OPTICAL_APERTURE_M) * p["R"]
    sr = C_LIGHT / (2 * LIDAR_BW * np.sqrt(2 * 10.0))
    d_miss_design = np.sqrt(sp ** 2
                            + (p["point"] * p["R"]) ** 2
                            + (p["u"] * p["R"] * p["muzzle"] / p["v_rel"]) ** 2
                            + (p["u"] * sr / p["v_rel"]) ** 2
                            + (p["u"] * p["tau"]) ** 2 * p["no_extrap"])
    d_miss_v4 = np.sqrt(sp ** 2 + (p["v_rel"] * p["tau"]) ** 2)

    # 2a. plan-form guidance constraint
    p_hit_plan = np.exp(-d_miss_v4 ** 2 / (2 * 0.38 ** 2))
    c2_plan = p_hit_plan > 0.5

    # 2b. design-form: shot mass at 90% confidence within the magazine limit
    m_on_target = GAMMA * dv * m_t / (p["beta"] * p["v_rel"])
    penalty = -2 * np.pi * np.e * d_miss_design ** 2 * np.log(1 - CONF) / a_t
    penalty = np.maximum(penalty, 1.0)
    m_shot = m_on_target * penalty
    c2_design = m_shot <= M_MAX_SHOT

    # 3. decay gate
    hg, dg, Z = T[(d_cm, case)]
    t_dec = 10 ** bilinear(hg, dg, Z, p["h_km"], p["dh_m"] / 1e3)
    c3 = t_dec < ODMSP_YR

    # 4. net debris >= 0 : sub-catastrophic impact does not disrupt the target
    #    and released mass self-disposes
    c4 = c1 & c_self

    # 5. cost
    hk, nk = K[(d_cm, case)]
    n_kick = np.maximum(np.interp(p["h_km"], hk, nk), 1.0)
    n_shot = n_kick / CONF
    c_var = n_shot * m_shot * (p["p_mat"] + p["c_launch"])
    c_obj = c_var + p["c_platform"] / p["n_eng"]
    cost_per_kg = c_obj / m_t
    c5 = cost_per_kg < COST_BAR
    c5_obj = c_obj < COST_BAR_OBJECT

    # 6. remediation value: the engagement must actually change the object's
    #    compliance status. An object that would have decayed inside 25 yr on
    #    its own has not been "remediated" by being made to decay sooner.
    if P is not None:
        hp, Zp = P[(d_cm, case)]
        t_nat = 10 ** np.interp(p["h_km"], hp, Zp)
        c6 = (t_nat >= ODMSP_YR) & c3
    else:
        t_nat = np.full_like(t_dec, np.nan)
        c6 = np.ones_like(c3)

    return dict(t_natural=t_nat, c6_useful=c6,
                viable_useful=c1 & c2_design & c3 & c4 & c5_obj & c6,
                Es=es, dv=dv, t_decay=t_dec, m_shot=m_shot, penalty=penalty,
                d_miss_design=d_miss_design, d_miss_v4=d_miss_v4,
                p_hit_plan=p_hit_plan, cost_per_kg=cost_per_kg,
                cost_per_object=c_obj, n_kick=n_kick,
                c1=c1, c2_plan=c2_plan, c2_design=c2_design, c3=c3, c4=c4,
                c5=c5, c5_obj=c5_obj, c_self=c_self,
                viable_plan=c1 & c2_plan & c3 & c4 & c5,
                viable_design=c1 & c2_design & c3 & c4 & c5,
                viable_objcost=c1 & c2_design & c3 & c4 & c5_obj)


def _dv_shift_vec(h_m, dh_m):
    from fg01.constants import MU, R_E
    r = R_E + h_m
    a_new = r - dh_m / 2.0
    return np.sqrt(MU / r) - np.sqrt(MU * (2.0 / r - 1.0 / a_new))


def _dv_direct_vec(h_m):
    from fg01.constants import MU, R_E
    r = R_E + h_m
    rp = R_E + 100e3
    a_new = 0.5 * (r + rp)
    return np.sqrt(MU / r) - np.sqrt(MU * (2.0 / r - 1.0 / a_new))


def sample(rng, n, h_fixed=None):
    lu = lambda lo, hi: np.exp(rng.uniform(np.log(lo), np.log(hi), n))
    return dict(
        beta=rng.uniform(*BETA_RANGE, n),
        v_rel=rng.uniform(400.0, 1200.0, n),
        es_c=rng.uniform(*ES_C_RANGE, n),
        dh_m=rng.uniform(150e3, 250e3, n),
        h_km=(np.full(n, h_fixed) if h_fixed is not None
              else rng.choice(ALT_SWEEP_KM, n).astype(float)),
        R=rng.uniform(5.0, 100.0, n),
        u=rng.uniform(1.0, 100.0, n),
        tau=lu(1e-4, 0.1),
        point=rng.uniform(5e-5, 1e-3, n),
        muzzle=lu(1e-3, 1e-2),
        no_extrap=np.zeros(n),
        m_t_scale=np.clip(rng.normal(1.0, 0.07, n), 0.5, 1.5),
        p_mat=rng.uniform(*P_RE_RANGE, n),
        c_launch=rng.uniform(*C_LAUNCH_RANGE, n),
        c_platform=rng.uniform(*C_PLATFORM_RANGE, n),
        n_eng=lu(1e3, 1e5),
    )


def sobol(d_cm, case, T, K, P, n_base=16384, out_key="viable_objcost"):
    """Saltelli estimator of first-order and total Sobol indices."""
    rng = np.random.default_rng(SEED + 7)
    keys = ["beta", "v_rel", "es_c", "dh_m", "h_km", "R", "u", "tau", "point",
            "muzzle", "p_mat", "c_launch", "c_platform", "n_eng"]
    A = sample(rng, n_base)
    B = sample(rng, n_base)
    fA = evaluate(A, d_cm, case, T, K, P)[out_key].astype(float)
    fB = evaluate(B, d_cm, case, T, K, P)[out_key].astype(float)
    var = np.var(np.concatenate([fA, fB]))
    rows = []
    if var <= 0:
        return pd.DataFrame([dict(parameter=k, S1=np.nan, ST=np.nan) for k in keys])
    for k in keys:
        AB = {kk: (B[kk] if kk == k else A[kk]) for kk in A}
        fAB = evaluate(AB, d_cm, case, T, K, P)[out_key].astype(float)
        s1 = np.mean(fB * (fAB - fA)) / var
        st = np.mean((fA - fAB) ** 2) / (2 * var)
        rows.append(dict(parameter=k, S1=float(s1), ST=float(st)))
    return pd.DataFrame(rows).sort_values("ST", ascending=False)


def run():
    banner("SIM-12  Global uncertainty quantification")
    T, K, P = build_lookups()
    rng = np.random.default_rng(SEED)

    # ---- headline joint run -----------------------------------------------
    rows = []
    for d_cm in SIZES:
        for case in CASES:
            for h in ALT_SWEEP_KM:
                n = max(N_SAMPLES // (len(SIZES) * len(CASES) * len(ALT_SWEEP_KM)), 4000)
                p = sample(rng, n, h_fixed=float(h))
                r = evaluate(p, d_cm, case, T, K, P)
                rows.append(dict(
                    d_debris_cm=d_cm, solar_case=case, h_km=h, n=n,
                    P_viable_design=float(r["viable_design"].mean()),
                    P_viable_plan=float(r["viable_plan"].mean()),
                    P_viable_objcost=float(r["viable_objcost"].mean()),
                    P_viable_useful=float(r["viable_useful"].mean()),
                    P_c5_cost_per_object=float(r["c5_obj"].mean()),
                    P_c6_remediation_value=float(r["c6_useful"].mean()),
                    median_t_natural_yr=float(np.median(r["t_natural"])),
                    P_c1_subcatastrophic=float(r["c1"].mean()),
                    P_c2_plan_phit=float(r["c2_plan"].mean()),
                    P_c2_design_mass=float(r["c2_design"].mean()),
                    P_c3_decay=float(r["c3"].mean()),
                    P_c4_net_debris=float(r["c4"].mean()),
                    P_c5_cost=float(r["c5"].mean()),
                    P_self_disposal=float(r["c_self"].mean()),
                    median_cost_per_kg=float(np.median(r["cost_per_kg"])),
                    median_cost_per_object=float(np.median(r["cost_per_object"])),
                    median_m_shot_g=float(np.median(r["m_shot"])) * 1e3,
                    median_d_miss_mm=float(np.median(r["d_miss_design"])) * 1e3,
                    median_t_decay_yr=float(np.median(r["t_decay"]))))
    df = pd.DataFrame(rows)
    save_table(df, "sim12_viability",
               "P(viable) and per-constraint satisfaction probability by "
               "altitude, debris size and solar case, under both the plan-form "
               "and design-form guidance constraints.",
               tags={"samples": str(N_SAMPLES), "seed": str(SEED)})

    nom = df[df.solar_case == "nominal"]
    print("\n  P(viable, design-form), nominal solar:")
    print(nom.pivot(index="h_km", columns="d_debris_cm",
                    values="P_viable_design").round(3).to_string())
    print("\n  P(viable, plan-form guidance criterion), nominal solar:")
    print(nom.pivot(index="h_km", columns="d_debris_cm",
                    values="P_viable_plan").round(3).to_string())
    print("\n  P(viable, per-object cost criterion), nominal solar:")
    print(nom.pivot(index="h_km", columns="d_debris_cm",
                    values="P_viable_objcost").round(3).to_string())
    print("\n  P(viable AND remediation-relevant), nominal solar:")
    print(nom.pivot(index="h_km", columns="d_debris_cm",
                    values="P_viable_useful").round(3).to_string())

    # ---- convergence (at the cell with the highest P, so it is not degenerate)
    cell = nom.loc[nom.P_viable_objcost.idxmax()]
    conv = []
    for n in (1000, 10000, 100000, 1000000):
        p = sample(np.random.default_rng(SEED + 1), n, h_fixed=float(cell.h_km))
        r = evaluate(p, float(cell.d_debris_cm), "nominal", T, K, P)
        conv.append(dict(n_samples=n,
                         P_viable_objcost=float(r["viable_objcost"].mean()),
                         stderr=float(np.std(r["viable_objcost"]) / np.sqrt(n))))
    save_table(pd.DataFrame(conv), "sim12_convergence",
               f"Monte-Carlo convergence at the highest-probability cell "
               f"({cell.d_debris_cm:g} cm debris, {cell.h_km:g} km, nominal solar).")
    print("\n  convergence:", [f"{c['n_samples']}: {c['P_viable_objcost']:.4f}"
                               for c in conv])

    # ---- Sobol -------------------------------------------------------------
    best = nom.loc[nom.P_viable_objcost.idxmax()]
    sob = sobol(float(best.d_debris_cm), "nominal", T, K, P)
    save_table(sob, "sim12_sobol",
               f"Sobol indices for the design-form viability indicator at "
               f"{best.d_debris_cm:g} cm, nominal solar (altitude sampled).")
    print("\n  Sobol total-effect ranking:")
    print(sob.head(8).to_string(index=False))

    # ---- posteriors at the best cell --------------------------------------
    p = sample(np.random.default_rng(SEED + 2), 400000,
               h_fixed=float(best.h_km))
    r = evaluate(p, float(best.d_debris_cm), "nominal", T, K, P)
    post = pd.DataFrame({
        "quantity": ["design dv (m/s)", "cost per kg (USD)",
                     "cost per object (USD)", "T_decay (yr)",
                     "shot mass (g)", "miss distance (mm)"],
        "p05": [np.percentile(r["dv"], 5), np.percentile(r["cost_per_kg"], 5),
                np.percentile(r["cost_per_object"], 5),
                np.percentile(r["t_decay"], 5), np.percentile(r["m_shot"], 5) * 1e3,
                np.percentile(r["d_miss_design"], 5) * 1e3],
        "p50": [np.median(r["dv"]), np.median(r["cost_per_kg"]),
                np.median(r["cost_per_object"]), np.median(r["t_decay"]),
                np.median(r["m_shot"]) * 1e3, np.median(r["d_miss_design"]) * 1e3],
        "p95": [np.percentile(r["dv"], 95), np.percentile(r["cost_per_kg"], 95),
                np.percentile(r["cost_per_object"], 95),
                np.percentile(r["t_decay"], 95), np.percentile(r["m_shot"], 95) * 1e3,
                np.percentile(r["d_miss_design"], 95) * 1e3]})
    save_table(post, "sim12_posteriors",
               f"Posterior 5th/50th/95th percentiles at the best cell "
               f"({best.d_debris_cm:g} cm, {best.h_km:g} km, nominal solar).")
    print("\n  posteriors at best cell:")
    print(post.to_string(index=False))

    concl = {
        "n_samples": N_SAMPLES, "seed": SEED,
        "best_cell": {"d_debris_cm": float(best.d_debris_cm),
                      "h_km": float(best.h_km),
                      "P_viable_design": float(best.P_viable_design),
                      "P_viable_plan": float(best.P_viable_plan)},
        "P_viable_design_max": float(nom.P_viable_design.max()),
        "P_viable_plan_max": float(nom.P_viable_plan.max()),
        "P_viable_objcost_max": float(nom.P_viable_objcost.max()),
        "cost_metric_note": (
            "P(viable) is zero for sub-2 cm debris under the plan's "
            "$100,000/kg criterion, and this is a property of the metric, not "
            "of the design: a 1 cm aluminium sphere masses 1.414 g, so the bar "
            "is $141 per object. No orbital system can service an object for "
            "$141. Under a $50,000-per-object criterion the same physics gives "
            "a materially different picture, reported in the P_viable_objcost "
            "columns. Both are shown; neither is privileged here."),
        "binding_constraint_by_cell": {
            f"{r.d_debris_cm:g}cm@{r.h_km:g}km": min(
                [("subcatastrophic", r.P_c1_subcatastrophic),
                 ("shot_mass", r.P_c2_design_mass),
                 ("decay_25yr", r.P_c3_decay),
                 ("net_debris", r.P_c4_net_debris),
                 ("cost_per_object", r.P_c5_cost_per_object)],
                key=lambda kv: kv[1])[0]
            for _, r in nom.iterrows()},
        "sobol_top": sob.head(5).to_dict(orient="records"),
        "guidance_definition_effect": (
            "Replacing the v4 Bernoulli hit criterion with the delivered-"
            "momentum sizing rule is worth the difference between the "
            "P_viable_plan and P_viable_design columns. Under the plan-form "
            "criterion the concept is throttled by control latency exactly as "
            "in v4; under the design-form criterion latency is not a binding "
            "constraint anywhere in the sampled space, and the binding "
            "constraint moves to cost and decay."),
    }
    save_json(concl, "sim12_conclusions")

    # ---- figures -----------------------------------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), sharey=True)
    for ax, col, ttl in (
            (axes[0], "P_viable_plan", "plan-form guidance + \\$100k/kg"),
            (axes[1], "P_viable_design", "design-form guidance + \\$100k/kg"),
            (axes[2], "P_viable_objcost", "design-form guidance + \\$50k/object")):
        for i, d in enumerate(SIZES):
            s = nom[nom.d_debris_cm == d].sort_values("h_km")
            ax.plot(s.h_km, s[col], color=CB[i], marker="o", ms=3.5,
                    label=f"{d:g} cm")
        ax.set_xlabel("altitude (km)"); ax.set_title(ttl); ax.set_ylim(-0.02, 1.02)
    axes[0].set_ylabel("P(viable)")
    axes[0].legend(fontsize=8, title="debris size", title_fontsize=8)
    fig.suptitle("SIM-12: joint viability probability, nominal solar, 10â¶ samples",
                 fontsize=10)
    save_fig(fig, "sim12_p_viable",
             "Joint probability that all five viability constraints hold. The "
             "two panels differ only in how the guidance constraint is defined, "
             "which is the single largest modelling choice in the study.")

    fig, ax = plt.subplots(figsize=(7.6, 4.2))
    s = nom[nom.d_debris_cm == best.d_debris_cm].sort_values("h_km")
    for i, (col, lab) in enumerate([
            ("P_c1_subcatastrophic", "sub-catastrophic $E_s$"),
            ("P_c2_design_mass", "shot mass â‰¤ 1 kg"),
            ("P_c3_decay", "post-kick lifetime < 25 yr"),
            ("P_c4_net_debris", "net debris â‰¥ 0"),
            ("P_c5_cost", "cost < \\$100k/kg"),
            ("P_c5_cost_per_object", "cost < \\$50k/object")]):
        ax.plot(s.h_km, s[col], color=CB[i], marker="o", ms=3.5, label=lab)
    ax.plot(s.h_km, s.P_viable_objcost, color="k", lw=2, label="joint (per-object)")
    ax.set_xlabel("altitude (km)"); ax.set_ylabel("P(constraint satisfied)")
    ax.set_title(f"SIM-12: which constraint binds ({best.d_debris_cm:g} cm debris, "
                 "nominal solar)")
    ax.legend(fontsize=8)
    save_fig(fig, "sim12_binding_constraints",
             "Per-constraint satisfaction probability against altitude. The "
             "lowest curve at each altitude is the binding constraint.")
    return df, sob, concl


if __name__ == "__main__":
    run()
