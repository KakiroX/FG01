"""
SIM-6 -- Electromagnetic launch technology assessment.

Demonstrated benchmarks only (AM-9). The question is whether the (mass,
velocity) point v5 needs sits inside the envelope that has actually been built.
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from fg01.constants import (EML_BENCHMARKS, COILGUN_ETA_RANGE, P_ELEC, RHO_CU,
                            GAMMA, RHO_AL)
from fg01.orbital import mass_sphere, dv_shift
from fg01.io_utils import save_table, save_fig, save_json, banner, CB

G0 = 9.80665


def coilgun_design(m_p, v, accel_g=10000.0, eta=0.4, l_stage=0.30,
                   bore_d=0.05):
    """
    First-order coilgun sizing.  Barrel length from the acceleration limit,
    stage count from stage length, capacitor mass from pulse-capacitor energy
    density, copper mass from a single-layer solenoid per stage.
    """
    a = accel_g * G0
    l_barrel = v ** 2 / (2 * a)
    n_stage = max(1, int(np.ceil(l_barrel / l_stage)))
    e_kin = 0.5 * m_p * v ** 2
    e_elec = e_kin / eta
    # pulse capacitors: 1.5 J/g is representative of current metallised-film
    # pulse-discharge capacitors
    m_cap = e_elec / 1.5e3
    # copper: single-layer solenoid, bore d, N turns of wire cross-section A_w
    turns_per_m = 200.0
    a_wire = 3.31e-6            # AWG 12, 3.31 mm^2
    l_wire = np.pi * bore_d * turns_per_m * l_barrel
    m_cu = l_wire * a_wire * RHO_CU
    return dict(l_barrel_m=l_barrel, n_stages=n_stage, E_kin_J=e_kin,
                E_elec_J=e_elec, m_cap_kg=m_cap, m_cu_kg=m_cu,
                l_wire_m=l_wire, wire_awg=12,
                energy_cost_usd=e_elec / 3.6e6 * P_ELEC)


def run():
    banner("SIM-6  Electromagnetic launch technology assessment")

    bench = pd.DataFrame([dict(system=b[0], m_p_kg=b[1], v_m_s=b[2],
                               E_muzzle_J=b[3], citation=b[4],
                               E_check_J=0.5 * b[1] * b[2] ** 2)
                          for b in EML_BENCHMARKS])
    save_table(bench, "sim06_demonstrated_benchmarks",
               "Demonstrated electromagnetic-launch benchmarks (AM-9: no "
               "aspirational systems).")

    rows = []
    for m_g in (0.01, 0.1, 1.0, 2.69, 10.0, 100.0, 175.0, 1000.0):
        for v in (400.0, 619.0, 1000.0, 2000.0, 3000.0, 5000.0):
            for eta in (0.30, 0.40, 0.50):
                d = coilgun_design(m_g * 1e-3, v, eta=eta)
                rows.append(dict(m_p_g=m_g, v_m_s=v, eta=eta, **d))
    df = pd.DataFrame(rows)
    save_table(df, "sim06_coilgun_designs",
               "Coilgun sizing across the projectile mass/velocity grid: "
               "barrel length at a 10,000 g acceleration limit, stage count, "
               "stored electrical energy, capacitor and copper mass.",
               tags={"eta": "SOURCED range 0.30-0.50",
                     "cap energy density": "UNVALIDATED (1.5 J/g)"})

    # gap analysis against the demonstrated envelope
    req = []
    for m_g, lab in ((2.69, "SIM-4 design (10 m range, 200 µrad)"),
                     (175.0, "10 cm miss tolerance"),
                     (0.242, "mass on target only (no pattern)")):
        e_kin = 0.5 * (m_g * 1e-3) * 619.0 ** 2
        best = bench.E_muzzle_J.max()
        req.append(dict(case=lab, m_p_g=m_g, v_m_s=619.0, E_muzzle_J=e_kin,
                        demonstrated_max_J=best,
                        margin_x=best / e_kin,
                        inside_envelope=bool(e_kin <= best)))
    dfr = pd.DataFrame(req)
    save_table(dfr, "sim06_gap_analysis",
               "Required muzzle energy versus the largest demonstrated "
               "electromagnetic launch.")
    print(dfr.to_string(index=False))

    d_des = coilgun_design(2.69e-3, 619.0, eta=0.4)
    d_big = coilgun_design(175e-3, 619.0, eta=0.4)
    concl = {
        "design_shot": {"m_p_g": 2.69, "v_m_s": 619.0, **{k: float(v) for k, v in d_des.items()}},
        "large_pattern_shot": {"m_p_g": 175.0, "v_m_s": 619.0, **{k: float(v) for k, v in d_big.items()}},
        "gap_statement": (
            "Required muzzle energy at the v5 design point is 515 J for the "
            "2.69 g pattern and 33.5 kJ for a 175 g pattern. The largest "
            "demonstrated electromagnetic launch is 33 MJ (US Navy EMRG, 2010), "
            "and 2 kJ-class coilguns are routine laboratory hardware. The v5 "
            "launcher requirement is 3-5 orders of magnitude inside the "
            "demonstrated envelope -- the opposite of v4's position, because "
            "the pivot cut both the required velocity and the required mass."),
        "energy_cost_per_shot_usd": float(d_des["energy_cost_usd"]),
        "charge_time_at_5kW_s": float(d_des["E_elec_J"] / 5000.0),
    }
    save_json(concl, "sim06_conclusions")
    print(f"\n  design shot: {d_des['E_kin_J']:.0f} J muzzle, "
          f"{d_des['E_elec_J']:.0f} J stored, barrel {d_des['l_barrel_m']:.2f} m, "
          f"{d_des['n_stages']} stages, Cu {d_des['m_cu_kg']:.2f} kg, "
          f"caps {d_des['m_cap_kg']:.2f} kg")
    print(f"  charge time at 5 kW: {concl['charge_time_at_5kW_s']:.2f} s; "
          f"electricity cost {concl['energy_cost_per_shot_usd']:.2e} USD/shot")

    fig, ax = plt.subplots()
    mm = np.geomspace(1e-3, 10, 100)
    for i, e in enumerate((1e3, 1e4, 1e5, 1e6, 3.3e7)):
        ax.loglog(mm, np.sqrt(2 * e / mm), color=CB[i], lw=1.2,
                  label=f"{e:.0e} J muzzle energy")
    ax.scatter(bench.m_p_kg, bench.v_m_s, c="k", zorder=5, s=28,
               label="demonstrated")
    ax.scatter([2.69e-3, 175e-3], [619, 619], c=CB[1], marker="*", s=160,
               zorder=6, label="FG01 v5 requirement")
    ax.set_xlabel("projectile mass (kg)")
    ax.set_ylabel("muzzle velocity (m/s)")
    ax.set_ylim(100, 1e4)
    ax.set_title("SIM-6: v5 launcher requirement vs. demonstrated envelope")
    ax.legend(fontsize=7)
    save_fig(fig, "sim06_envelope",
             "Required launcher operating point against demonstrated "
             "electromagnetic launches. The v5 requirement sits far inside "
             "hardware that already exists, which is a direct consequence of "
             "the reduced Delta-v and the low intercept velocity the energy "
             "gate prefers.")
    return df, concl


if __name__ == "__main__":
    run()
