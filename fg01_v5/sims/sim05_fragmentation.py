"""
SIM-5 -- Fragmentation and secondary-debris hazard.

Confirms (rather than assumes) that the v5 sub-catastrophic design point does
not generate net debris, and quantifies what a sub-catastrophic swarm still
produces.

Three fragment estimates are reported side by side, because they disagree and
the disagreement is itself the finding:

  (a) NASA SBM, non-catastrophic branch, M_eff = m_p * v[km/s]  -- as mandated
      by the plan. At the v5 design point this is an extrapolation ~4x below
      the model's calibration velocity and applied to a target only ~9x the
      projectile mass; it returns a non-zero count of >1 cm fragments from a
      1 cm target, which is not physically possible without disruption.
  (b) NASA SBM, catastrophic branch -- what v4 used, shown for comparison.
  (c) An energy-based cratering bound: crater volume = KE / Y with Y the
      dynamic cratering resistance of Al-6061, giving ejecta mass directly and
      bounding the largest liberated fragment by the crater size.

(c) is the physically applicable one at 0.5-1.1 km/s; (a) is reported because
the plan mandates the SBM and because it is the conservative direction.
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from fg01.constants import RHO_AL, RHO_RE, GAMMA, ES_C, ES_C_RANGE, ALT_SWEEP_KM
from fg01.orbital import mass_sphere, dv_shift, area_to_mass_sphere
from fg01.sbm import n_fragments, specific_energy, crater_ejecta_mass
from fg01.io_utils import save_table, save_fig, save_json, banner, CB

# dynamic cratering resistance of Al-6061 (energy per unit crater volume).
# Laboratory crater-volume scaling for ductile aluminium gives 3e8 - 1e9 J/m^3;
# central value 5e8.  SOURCED range / UNVALIDATED central value.
Y_AL = 5.0e8
Y_AL_RANGE = (3.0e8, 1.0e9)


def run():
    banner("SIM-5  Fragmentation and secondary-debris hazard")
    m_t = float(mass_sphere(0.01, RHO_AL))
    v_grid = [400.0, 619.0, 800.0, 1000.0, 1500.0, 2000.0]
    rows = []
    for h in ALT_SWEEP_KM:
        dv = float(dv_shift(h * 1e3, 200e3))
        for beta in (1.0, 1.5, 2.5):
            for v in v_grid:
                m_p = GAMMA * dv * m_t / (beta * v)
                es = float(specific_energy(m_p, v, m_t))
                for esc in ES_C_RANGE + (ES_C,):
                    n1cm, m_eff, cat, _ = n_fragments(m_p, v, m_t, 0.01, esc)
                    n1mm, _, _, _ = n_fragments(m_p, v, m_t, 0.001, esc)
                    ej, ej_rng = crater_ejecta_mass(m_p, v)
                    ke = 0.5 * m_p * v ** 2
                    v_crater = ke / Y_AL
                    m_crater = v_crater * RHO_AL
                    rows.append(dict(
                        h_km=h, beta=beta, v_rel_m_s=v, Es_c_kJ_kg=esc / 1e3,
                        m_p_g=m_p * 1e3, Es_kJ_kg=es / 1e3,
                        catastrophic=bool(cat), M_eff_kg=float(m_eff),
                        N_gt_1cm_SBM=float(n1cm), N_gt_1mm_SBM=float(n1mm),
                        ejecta_g_gault=float(ej) * 1e3,
                        crater_volume_mm3=v_crater * 1e9,
                        ejecta_g_energy_bound=m_crater * 1e3,
                        ejecta_frac_of_target=m_crater / m_t,
                        KE_J=ke))
    df = pd.DataFrame(rows)
    save_table(df, "sim05_fragmentation_full",
               "Fragment production under the NASA SBM (regime-selected "
               "effective mass, AM-13) and under an energy-based cratering "
               "bound, across altitude, beta, intercept velocity and the "
               "catastrophic-threshold range.",
               tags={"SBM": "SOURCED(johnson2001)", "Y_Al": "UNVALIDATED",
                     "calibration": "SBM extrapolated below its fit velocity"})

    des = df[(df.h_km == 800) & (df.beta == 1.5) & (df.v_rel_m_s == 619.0)
             & (df.Es_c_kJ_kg == 40.0)].iloc[0]
    print(f"  design point (800 km, β=1.5, 619 m/s): m_p={des.m_p_g:.3f} g, "
          f"E_s={des.Es_kJ_kg:.1f} kJ/kg, catastrophic={des.catastrophic}")
    print(f"    SBM non-catastrophic branch: N(>1 cm)={des.N_gt_1cm_SBM:.3f}, "
          f"N(>1 mm)={des.N_gt_1mm_SBM:.1f}")
    print(f"    energy-based crater: {des.crater_volume_mm3:.1f} mm³, ejecta "
          f"{des.ejecta_g_energy_bound:.3f} g "
          f"({des.ejecta_frac_of_target*100:.1f}% of target mass)")

    # ---- per-grain (swarm) branch, AM-16 -----------------------------------
    sw = []
    for m_grain_mg in (0.01, 0.1, 1.0, 10.0, 100.0):
        m_grain = m_grain_mg * 1e-6
        for v in (619.0, 1000.0):
            es_i = float(specific_energy(m_grain, v, m_t))
            n1cm, _, cat, _ = n_fragments(m_grain, v, m_t, 0.01, ES_C)
            sw.append(dict(m_grain_mg=m_grain_mg, v_rel_m_s=v,
                           Es_per_impact_kJ_kg=es_i / 1e3,
                           per_impact_catastrophic=bool(cat),
                           frac_of_threshold=es_i / ES_C,
                           N_gt_1cm_per_impact_SBM=float(n1cm)))
    dfw = pd.DataFrame(sw)
    save_table(dfw, "sim05_swarm_per_impact",
               "AM-16: per-impact specific energy for the swarm. Reported "
               "ALWAYS alongside the bulk figure, never instead of it. The "
               "bulk (single-impactor) energy at the v5 design point is "
               "already sub-catastrophic, so this branch is a robustness "
               "reserve, not a load-bearing assumption.",
               tags={"swarm equivalence": "UNVALIDATED (AM-15/AM-16)"})

    # ---- net debris ledger --------------------------------------------------
    # missed mass persistence comes from SIM-2: at v_rel >= dv_direct the missed
    # grains are on immediate-reentry trajectories, so they contribute zero
    # long-lived objects.
    ledger = []
    for v in v_grid:
        for beta in (1.0, 1.5, 2.5):
            r = df[(df.h_km == 800) & (df.beta == beta) & (df.v_rel_m_s == v)
                   & (df.Es_c_kJ_kg == 40.0)].iloc[0]
            n_sbm = r.N_gt_1cm_SBM
            n_phys = 0.0 if not r.catastrophic else r.N_gt_1cm_SBM
            ledger.append(dict(
                v_rel_m_s=v, beta=beta, catastrophic=bool(r.catastrophic),
                removed=1.0,
                new_gt1cm_SBM_branch=float(n_sbm),
                new_gt1cm_physical=float(n_phys),
                missed_mass_persistent_objects=0.0,
                net_SBM=1.0 - n_sbm,
                net_physical=1.0 - n_phys))
    dfl = pd.DataFrame(ledger)
    save_table(dfl, "sim05_net_debris_ledger",
               "Net debris per engagement at 800 km. 'SBM branch' applies the "
               "mandated model verbatim; 'physical' recognises that a "
               "sub-catastrophic impact cannot liberate a fragment larger than "
               "the target it did not disrupt.")
    print("\n  net debris ledger (800 km):")
    print(dfl.to_string(index=False))

    concl = {
        "design_point_subcatastrophic": bool(not des.catastrophic),
        "design_Es_kJ_kg": float(des.Es_kJ_kg),
        "fraction_of_threshold": float(des.Es_kJ_kg / 40.0),
        "v_rel_at_which_design_turns_catastrophic_m_s": {
            f"beta={b}": float(df[(df.h_km == 800) & (df.beta == b)
                                  & (df.Es_c_kJ_kg == 40.0)
                                  & df.catastrophic].v_rel_m_s.min())
            if len(df[(df.h_km == 800) & (df.beta == b)
                      & (df.Es_c_kJ_kg == 40.0) & df.catastrophic]) else None
            for b in (1.0, 1.5, 2.5)},
        "sbm_non_catastrophic_N_gt_1cm_at_design": float(des.N_gt_1cm_SBM),
        "sbm_caveat": (f"The SBM non-catastrophic branch returns "
                       f"{des.N_gt_1cm_SBM:.3f} fragments >1 cm from a 1 cm, "
                       f"{m_t*1e3:.2f} g target that by construction is not "
                       "disrupted. This is the model being extrapolated outside "
                       "its calibration (hypervelocity impacts, large targets). "
                       "It is reported because it is the mandated model and "
                       "because it errs conservatively."),
        "energy_based_ejecta_g": float(des.ejecta_g_energy_bound),
        "energy_based_ejecta_fraction_of_target": float(des.ejecta_frac_of_target),
        "ejecta_fate": ("Cratering ejecta leave the target at a small fraction "
                        "of the 619 m/s impact speed and have area-to-mass "
                        "ratios 5-50x the parent, so they decay faster than the "
                        "parent object on the same lowered orbit (SIM-2)."),
        "missed_mass_contribution": ("Zero persistent objects: at v_rel >= "
                                     "dv_direct(h) the missed mass is on an "
                                     "immediate-reentry trajectory (SIM-2)."),
        "net_debris_per_engagement": {"SBM_branch": float(dfl[(dfl.v_rel_m_s == 619)
                                      & (dfl.beta == 1.5)].net_SBM.iloc[0]),
                                      "physical": 1.0},
    }
    save_json(concl, "sim05_conclusions")

    fig, ax = plt.subplots()
    for i, beta in enumerate((1.0, 1.5, 2.5)):
        s = df[(df.h_km == 800) & (df.beta == beta)
               & (df.Es_c_kJ_kg == 40.0)].sort_values("v_rel_m_s")
        ax.plot(s.v_rel_m_s, s.Es_kJ_kg, color=CB[i], marker="o", ms=3,
                label=f"β={beta}")
    ax.axhline(40, color="k", lw=1.3)
    ax.axhspan(40, 200, color="#D55E00", alpha=0.10)
    ax.axhspan(30, 50, color="k", alpha=0.06)
    ax.axvline(619, color="0.4", ls=":", lw=1)
    ax.text(625, 5, "design point", fontsize=7.5, color="0.3")
    ax.set_xlabel("intercept velocity $v_{rel}$ (m/s)")
    ax.set_ylabel("bulk specific energy $E_s$ (kJ/kg)")
    ax.set_ylim(0, 130)
    ax.set_title("SIM-5: distance to the catastrophic threshold at 800 km "
                 "(grey band = 30–50 kJ/kg threshold uncertainty)")
    ax.legend(fontsize=8)
    save_fig(fig, "sim05_es_vs_vrel",
             "Bulk specific energy against intercept velocity. The design "
             "point sits at roughly half the catastrophic threshold; the "
             "margin closes at high v_rel and low beta, which is what sets the "
             "velocity ceiling used in SIM-13.")
    return df, dfl, concl


if __name__ == "__main__":
    run()
