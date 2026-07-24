"""
LAS-1 -- Laser-ablation engagement model, and the elevation tension.

The v5/v6 comparison took published laser costs of $100-500/kg at face value
and concluded FG01 loses by five orders of magnitude. That comparison was never
put on the same footing: FG01's engagement was built from physics while the
laser's was a quoted number. This module builds the laser engagement from
physics.

The central structural result is the laser's own version of INT-0. Ablation
recoil pushes the target along the line of sight, away from the illuminator.
For a ground station that direction is mostly radial, and only the component
opposing the target's velocity lowers the orbit. That component is
    retrograde fraction = -R_E sin(phi) / R
which is ZERO at zenith -- where the range is shortest, the airmass lowest and
the beam best -- and largest at low elevation, where the range is longest, the
airmass highest and the spot largest. The two requirements are in direct
opposition, and the optimum is a compromise that gets neither.
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from fg01.constants import R_E, RHO_AL, ES_C, P_ELEC
from fg01.orbital import mass_sphere, xsec_sphere, dv_shift, v_circ
from fg01.laser import (spot_diameter, airmass, atmospheric_transmission,
                        fried_parameter, ao_strehl, fluence,
                        coupling_coefficient, impulse_per_pulse,
                        pass_geometry, angular_rate, F_TH, F_OPT, CM_MAX)
from fg01.io_utils import save_table, save_fig, save_json, banner, CB

H = 900e3
D_TARGET = 0.01
M_T = float(mass_sphere(D_TARGET, RHO_AL))
A_T = float(xsec_sphere(D_TARGET))
DV_REQ = float(dv_shift(H, 200e3))

# Reference ground facility, at the optimistic end of what has been proposed
D_APERTURE = 5.0
LAM = 532e-9
E_PULSE = 20000.0           # J per pulse (ORION-class proposal)
PRF = 10.0                  # Hz
BEAM_Q = 1.3                # near-diffraction-limited high-energy laser


def run():
    banner("LAS-1  Laser-ablation engagement model and the elevation tension")

    # ---- the elevation tension --------------------------------------------
    rows = []
    for phi_deg in np.arange(-40.0, 0.01, 0.5):
        R, elev, retro = pass_geometry(H, np.radians(phi_deg))
        if elev < 5.0:
            continue
        X = float(airmass(elev))
        T = float(atmospheric_transmission(elev))
        r0 = float(fried_parameter(elev, lam=LAM))
        S = float(ao_strehl(D_APERTURE, r0))
        d_spot = float(spot_diameter(R, D_APERTURE, LAM, BEAM_Q, S))
        p, F, cm = impulse_per_pulse(E_PULSE, d_spot, A_T, T)
        rows.append(dict(
            phi_deg=phi_deg, slant_range_km=R / 1e3, elevation_deg=elev,
            retrograde_fraction=retro, airmass=X, transmission=T,
            r0_cm=r0 * 100, ao_strehl=S,
            spot_diameter_m=d_spot,
            fluence_J_cm2=float(F) / 1e4,
            above_threshold=bool(F >= F_TH),
            Cm_N_s_per_J=float(cm),
            impulse_N_s=float(p),
            useful_impulse_N_s=float(p) * max(retro, 0.0),
            dv_per_pulse_m_s=float(p) * max(retro, 0.0) / M_T,
            fraction_on_target=min(A_T / (np.pi * (d_spot / 2) ** 2), 1.0)))
    df = pd.DataFrame(rows)
    save_table(df, "las01_elevation_tension",
               "Ground-station engagement geometry: retrograde efficiency, "
               "range, airmass, adaptive-optics Strehl, spot size, fluence and "
               "useful impulse against position in the pass.",
               tags={"C_m": "SOURCED (Phipps laser-propulsion literature)",
                     "F_th": "SOURCED range 0.5-5 J/cm^2",
                     "r0, tau": "SOURCED (good astronomical site)"})

    show = df[df.phi_deg.isin([-30.0, -20.0, -15.0, -10.0, -7.0, -5.0, -3.0, -1.0])]
    print(show[["phi_deg", "elevation_deg", "slant_range_km", "retrograde_fraction",
                "airmass", "spot_diameter_m", "fluence_J_cm2", "above_threshold",
                "dv_per_pulse_m_s"]].to_string(index=False))

    ok = df[df.above_threshold & (df.retrograde_fraction > 0)]
    if len(ok):
        best = ok.loc[ok.dv_per_pulse_m_s.idxmax()]
        n_pulses = DV_REQ / best.dv_per_pulse_m_s
    else:
        best, n_pulses = None, np.inf
    print(f"\n  best geometry: elevation {best.elevation_deg:.1f} deg, "
          f"range {best.slant_range_km:.0f} km, retrograde fraction "
          f"{best.retrograde_fraction:.3f}")
    print(f"  spot {best.spot_diameter_m:.3f} m on a {D_TARGET*100:.0f} cm target "
          f"-> {best.fraction_on_target:.2e} of the beam intercepted")
    print(f"  dv per pulse {best.dv_per_pulse_m_s*1e3:.4f} mm/s "
          f"-> {n_pulses:,.0f} pulses for the {DV_REQ:.1f} m/s deorbit kick")

    # ---- what aperture is even needed to clear the threshold? -------------
    ap = []
    for D in (1.0, 2.0, 3.0, 5.0, 8.0, 10.0, 15.0):
        for E_p in (2000.0, 5000.0, 20000.0, 100000.0):
            R, elev, retro = pass_geometry(H, np.radians(-7.0))
            T = float(atmospheric_transmission(elev))
            r0 = float(fried_parameter(elev, lam=LAM))
            S = float(ao_strehl(D, r0))
            d_spot = float(spot_diameter(R, D, LAM, BEAM_Q, S))
            p, F, cm = impulse_per_pulse(E_p, d_spot, A_T, T)
            dvp = float(p) * retro / M_T
            ap.append(dict(
                aperture_m=D, pulse_energy_J=E_p, ao_strehl=S,
                spot_diameter_m=d_spot, fluence_J_cm2=float(F) / 1e4,
                above_threshold=bool(F >= F_TH),
                dv_per_pulse_mm_s=dvp * 1e3,
                pulses_needed=DV_REQ / dvp if dvp > 0 else np.inf,
                engagement_time_s=(DV_REQ / dvp / PRF) if dvp > 0 else np.inf,
                laser_energy_MJ=(DV_REQ / dvp * E_p / 1e6) if dvp > 0 else np.inf))
    dfa = pd.DataFrame(ap)
    save_table(dfa, "las01_aperture_scaling",
               "Aperture and pulse energy needed to clear the ablation "
               "threshold at the optimum engagement geometry (7 deg past the "
               "meridian, ~46 deg elevation, 1150 km).")
    print("\n  aperture / pulse-energy scaling at the optimum geometry:")
    print(dfa[dfa.pulse_energy_J == 20000.0][
        ["aperture_m", "ao_strehl", "spot_diameter_m", "fluence_J_cm2",
         "above_threshold", "pulses_needed", "engagement_time_s"]].to_string(index=False))

    # ---- energy and cost per object ---------------------------------------
    e_total = n_pulses * E_PULSE
    wall_plug = 0.2
    e_elec = e_total / wall_plug
    concl = {
        "engagement_geometry": {
            "optimum_elevation_deg": float(best.elevation_deg),
            "slant_range_km": float(best.slant_range_km),
            "retrograde_fraction": float(best.retrograde_fraction),
            "zenith_retrograde_fraction": 0.0,
            "statement": (
                "Ablation recoil is directed along the line of sight, so only "
                "-R_E sin(phi)/R of it opposes the target's velocity. That "
                "fraction is exactly zero at zenith -- the point of best range, "
                "lowest airmass and smallest spot -- and rises toward the "
                "horizon, where the range is 2-3x longer, the airmass 3-6x "
                "higher and the spot correspondingly larger. The laser cannot "
                "have good optics and useful geometry at the same time. This "
                "is the laser's exact analogue of INT-0's finding for FG01, and "
                "unlike FG01 it has no co-orbital escape: a ground station "
                "cannot fly alongside the target.")},
        "beam_waste": {
            "spot_diameter_m": float(best.spot_diameter_m),
            "target_diameter_m": D_TARGET,
            "fraction_intercepted": float(best.fraction_on_target),
            "areal_penalty": float(1.0 / best.fraction_on_target),
            "statement": (
                f"The spot is {best.spot_diameter_m/D_TARGET:.0f}x wider than "
                f"the target, so {(1-best.fraction_on_target)*100:.3f}% of every "
                "pulse misses. The laser pays the same areal penalty FG01 pays "
                f"-- {1.0/best.fraction_on_target:,.0f}x against FG01's 22x -- "
                "but pays it in photons, which are cheap, rather than in "
                "tungsten, which must be launched. That is the laser's real "
                "advantage, and it is an advantage about *consumables*, not "
                "about accuracy.")},
        "pulses_per_object": float(n_pulses),
        "engagement_time_s": float(n_pulses / PRF),
        "laser_energy_MJ": float(e_total / 1e6),
        "electricity_cost_usd": float(e_elec / 3.6e6 * P_ELEC),
        "threshold_is_a_cliff": (
            "Ablation is a threshold process: below F_th the surface is heated "
            "and nothing is ejected, so the impulse is exactly zero. A laser "
            "whose spot grows past the point where fluence falls below "
            "threshold does not deliver a smaller kick -- it delivers none. "
            "Quantified in LAS-2."),
    }
    save_json(concl, "las01_conclusions")
    print(f"\n  energy: {e_total/1e6:.1f} MJ of laser output per object, "
          f"${concl['electricity_cost_usd']:.2f} of electricity")

    # ---- figures -----------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.3))
    ax = axes[0]
    ax.plot(df.elevation_deg, df.retrograde_fraction, color=CB[0], lw=2,
            label="retrograde fraction (useful)")
    ax2 = ax.twinx()
    ax2.plot(df.elevation_deg, df.spot_diameter_m, color=CB[1], lw=2,
             label="spot diameter")
    ax2.plot(df.elevation_deg, df.airmass / 10, color=CB[2], lw=1.4, ls="--",
             label="airmass / 10")
    ax.set_xlabel("elevation angle (deg)")
    ax.set_ylabel("retrograde fraction", color=CB[0])
    ax2.set_ylabel("spot diameter (m) / airmassÃ·10", color=CB[1])
    ax.set_title("The tension: useful geometry vs. good optics")
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=7)
    ax = axes[1]
    ax.plot(df.elevation_deg, df.dv_per_pulse_m_s * 1e3, color=CB[3], lw=2)
    ax.axvline(best.elevation_deg, color="0.4", ls=":", lw=1.2)
    ax.text(best.elevation_deg + 1, df.dv_per_pulse_m_s.max() * 1e3 * 0.5,
            f"optimum\n{best.elevation_deg:.0f}Â°", fontsize=7.5, color="0.3")
    ax.set_xlabel("elevation angle (deg)")
    ax.set_ylabel("useful Î”v per pulse (mm/s)")
    ax.set_title(f"Net result: {n_pulses:,.0f} pulses per object")
    fig.suptitle("LAS-1: a ground laser cannot have useful geometry and good "
                 "optics at once", fontsize=10)
    save_fig(fig, "las01_elevation_tension",
             "Ablation recoil points away from the illuminator, so the useful "
             "retrograde component vanishes at zenith where the optics are "
             "best, and peaks near the horizon where range, airmass and spot "
             "size are worst. The product has a shallow optimum that achieves "
             "neither.")
    return df, concl


if __name__ == "__main__":
    run()
