"""
CHK-01 -- Is v_rel = 400 m/s consistent with the FG01 simulation gates?

400 m/s came out of the LAUNCHER studies (GUN-01, COIL-01), not out of the
engagement physics. The v5/v6 design point was 619 m/s. This module checks
400 m/s against every gate the earlier work established, so the launcher
choice does not silently break the mission analysis.

Gates checked:
  1. INT-2  beta at 400 m/s (it is velocity dependent -- this is the catch)
  2. SIM-0  disruption gate, E_s = gamma*dv*v_rel/(2 beta)
  3. SIM-2  self-disposal floor: missed mass must reenter, v_rel >= dv_direct(h)
  4. INT-0  firing cone: the plane-match tolerance scales with v_launch
  5. SIM-1  required mass on target, and what a 100 g shot actually delivers
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd

from fg01.constants import (R_E, RHO_AL, ES_C, ES_C_RANGE, GAMMA, ALT_SWEEP_KM,
                            RHO_RE)
from fg01.orbital import (dv_shift, dv_direct_reentry, mass_sphere,
                          xsec_sphere, v_circ)
from fg01.interaction import beta_envelope, mean_incidence_cos_sphere
from fg01.relmotion import v_rel_crossing
from fg01.delivery import penalty_for_confidence
from fg01.io_utils import save_table, save_json, banner

V_A, V_B = 400.0, 619.0
H_DES = 900e3
M_SHOT = 0.100          # total launched mass, COIL-01 / GUN-01
M_RE = 0.081            # rhenium fraction of it


def beta_at(v):
    lo, ctr, hi = beta_envelope(v)
    return (float(lo), float(ctr), float(hi),
            1.0 + (float(ctr) - 1.0) * mean_incidence_cos_sphere())


def run():
    banner("CHK-01  Is 400 m/s consistent with the FG01 gates?")

    # ---- 1. beta is velocity dependent -----------------------------------
    print("\n[1] INT-2 beta (the catch: beta FALLS with impact speed)\n")
    rows = []
    for v in (300.0, 350.0, 400.0, 500.0, 619.0, 700.0):
        lo, ctr, hi, obl = beta_at(v)
        rows.append(dict(v_rel=v, v_over_v_strength=v / 320.0,
                         beta_low=lo, beta_central=ctr, beta_oblique=obl,
                         mass_penalty_vs_619=(1.201 * 619.0) / (obl * v)))
    d1 = pd.DataFrame(rows)
    print(d1.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    b400 = beta_at(V_A)[3]
    b619 = beta_at(V_B)[3]
    print(f"\n    beta(400) = {b400:.3f} vs beta(619) = {b619:.3f}")
    print(f"    400 m/s is only {400/320:.2f}x the Al strength velocity, so the "
          f"impact is essentially pure embedding")

    # ---- 2. disruption gate ----------------------------------------------
    print("\n[2] SIM-0 disruption gate  E_s = gamma*dv*v_rel/(2 beta)\n")
    rows = []
    for h in ALT_SWEEP_KM:
        dv = float(dv_shift(h * 1e3, 200e3))
        for v, b in ((V_A, b400), (V_B, b619)):
            es = GAMMA * dv * v / (2 * b)
            rows.append(dict(h_km=h, v_rel=v, beta=b, dv_shift=dv,
                             Es_kJ_kg=es / 1e3,
                             frac_of_40=es / ES_C,
                             pass_40=bool(es < ES_C),
                             pass_30_worst=bool(es < ES_C_RANGE[0])))
    d2 = pd.DataFrame(rows)
    piv = d2.pivot(index="h_km", columns="v_rel", values="Es_kJ_kg")
    print("    E_s (kJ/kg) by altitude:")
    print(piv.to_string(float_format=lambda x: f"{x:.1f}"))
    print(f"\n    400 m/s passes 40 kJ/kg everywhere: {d2[d2.v_rel==V_A].pass_40.all()}")
    print(f"    400 m/s passes 30 kJ/kg (worst case) everywhere: "
          f"{d2[d2.v_rel==V_A].pass_30_worst.all()}")
    print(f"    619 m/s passes 30 kJ/kg everywhere: "
          f"{d2[d2.v_rel==V_B].pass_30_worst.all()}")

    # ---- 3. self-disposal floor ------------------------------------------
    print("\n[3] SIM-2 self-disposal: missed mass must reenter "
          "(v_rel >= dv_direct)\n")
    rows = []
    for h in ALT_SWEEP_KM:
        dvd = float(dv_direct_reentry(h * 1e3))
        rows.append(dict(h_km=h, dv_direct_m_s=dvd,
                         margin_400=V_A / dvd, margin_619=V_B / dvd,
                         ok_400=bool(V_A >= dvd), ok_619=bool(V_B >= dvd)))
    d3 = pd.DataFrame(rows)
    print(d3.to_string(index=False, float_format=lambda x: f"{x:.2f}"))
    worst = d3.loc[d3.margin_400.idxmin()]
    print(f"\n    tightest margin at {worst.h_km:.0f} km: "
          f"400/{worst.dv_direct_m_s:.0f} = {worst.margin_400:.2f}x "
          f"(619 m/s has {worst.margin_619:.2f}x)")

    # ---- 4. firing cone ---------------------------------------------------
    print("\n[4] INT-0 firing cone: plane-match tolerance scales with "
          "v_launch\n")
    v_orb = float(v_circ(H_DES))
    rows = []
    for v_l in (V_A, V_B):
        for pen in (1.05, 1.10, 1.25):
            # E_s penalty = 1 + (v_geom/v_launch)^2
            v_geom = v_l * np.sqrt(pen - 1.0)
            theta = 2 * np.arcsin(np.clip(v_geom / (2 * v_orb), 0, 1))
            rows.append(dict(v_launch=v_l, energy_penalty=pen,
                             v_geom_max_m_s=v_geom,
                             theta_max_deg=np.degrees(theta)))
    d4 = pd.DataFrame(rows)
    print(d4.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    t400 = float(d4[(d4.v_launch == V_A) & (d4.energy_penalty == 1.10)].theta_max_deg.iloc[0])
    t619 = float(d4[(d4.v_launch == V_B) & (d4.energy_penalty == 1.10)].theta_max_deg.iloc[0])
    print(f"\n    at a 10% energy penalty the cone TIGHTENS from "
          f"{t619:.2f} deg (619 m/s) to {t400:.2f} deg (400 m/s) "
          f"-- {100*(1-t400/t619):.0f}% tighter")

    # ---- 5. what a 100 g shot actually delivers ---------------------------
    print("\n[5] SIM-1 what a 100 g shot delivers, and its matched target\n")
    rows = []
    for v, b in ((V_A, b400), (V_B, b619)):
        p = b * M_SHOT * v                      # delivered momentum, N.s
        dv_req = float(dv_shift(H_DES, 200e3))
        m_matched = p / dv_req                  # target fully deorbited by 1 shot
        d_matched = 2 * (3 * m_matched / (4 * np.pi * RHO_AL)) ** (1 / 3)
        es_matched = 0.5 * M_SHOT * v ** 2 / m_matched
        rows.append(dict(v_rel=v, beta=b, momentum_Ns=p,
                         matched_target_kg=m_matched,
                         matched_diameter_cm=d_matched * 100,
                         Es_on_matched_kJ_kg=es_matched / 1e3,
                         subcatastrophic=bool(es_matched < ES_C)))
    d5 = pd.DataFrame(rows)
    print(d5.to_string(index=False, float_format=lambda x: f"{x:.3f}"))

    # mass needed for a 1 cm target, both velocities
    m_t = float(mass_sphere(0.01, RHO_AL))
    a_t = float(xsec_sphere(0.01))
    dv_req = float(dv_shift(H_DES, 200e3))
    rows = []
    for v, b in ((V_A, b400), (V_B, b619)):
        m_on = GAMMA * dv_req * m_t / (b * v)
        pen = float(penalty_for_confidence(0.005775, a_t, 0.95))
        rows.append(dict(v_rel=v, beta=b, m_on_target_mg=m_on * 1e6,
                         areal_penalty=pen, m_launch_g=m_on * pen * 1e3,
                         shots_from_100g=M_SHOT / (m_on * pen)))
    d6 = pd.DataFrame(rows)
    print("\n    for a 1 cm target at 900 km:")
    print(d6.to_string(index=False, float_format=lambda x: f"{x:.3f}"))

    concl = {
        "beta_400": b400, "beta_619": b619,
        "mass_penalty_400_vs_619": (b619 * V_B) / (b400 * V_A),
        "disruption_gate_400_pass_40": bool(d2[d2.v_rel == V_A].pass_40.all()),
        "disruption_gate_400_pass_30": bool(d2[d2.v_rel == V_A].pass_30_worst.all()),
        "disruption_gate_619_pass_30": bool(d2[d2.v_rel == V_B].pass_30_worst.all()),
        "self_disposal_400_ok_all_altitudes": bool(d3.ok_400.all()),
        "self_disposal_tightest_margin_400": float(d3.margin_400.min()),
        "firing_cone_deg_400": t400, "firing_cone_deg_619": t619,
        "matched_target_kg_400": float(d5[d5.v_rel == V_A].matched_target_kg.iloc[0]),
        "matched_diameter_cm_400": float(d5[d5.v_rel == V_A].matched_diameter_cm.iloc[0]),
    }
    save_json(concl, "chk01_conclusions")
    save_table(d1, "chk01_beta", "beta vs impact velocity")
    save_table(d2, "chk01_disruption", "disruption gate at 400 vs 619 m/s")
    save_table(d3, "chk01_self_disposal", "self-disposal floor by altitude")
    save_table(d4, "chk01_firing_cone", "firing cone vs launch velocity")
    save_table(d5, "chk01_matched_target", "matched target for a 100 g shot")
    return concl


if __name__ == "__main__":
    run()
