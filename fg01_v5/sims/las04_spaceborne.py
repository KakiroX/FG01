"""
LAS-4 -- The space-based laser: closing the strongest counterargument.

LAS-1..3 compared FG01 against a GROUND laser and found the laser 46x beyond
demonstrated hardware, fighting an elevation tension it cannot escape, and
pointing at 1.6x its own diffraction limit. The obvious rebuttal is: put the
laser in orbit. No atmosphere, no airmass, no turbulence, no terminator window,
no elevation tension. LIMITATIONS L26 flagged this as the largest gap in the
comparison and the first place a critic would attack.

This module closes it, and the answer is not the one that flatters FG01 most.

Three questions, in order:

  1. Does an orbiting laser escape the OPTICS problem?  YES, completely.
     Without atmosphere the spot is diffraction-limited, and at short range it
     shrinks below the target size, so essentially all the energy lands. Pulse
     energy collapses from ~9 kJ to a few joules and the 46x technology gap
     disappears. On the physics of beam delivery a space laser is far better
     than a ground laser AND better than FG01's projectile.

  2. Does it escape the DELIVERY-GEOMETRY problem?  NO.
     Ablation recoil is directed along the line of sight. To push a target
     retrograde the illuminator must therefore lie ahead of it along its own
     velocity vector -- which is exactly FG01's co-orbital leading geometry.
     A laser in a different plane or a different altitude sees a line of sight
     that is mostly radial or mostly cross-track, and its useful retrograde
     fraction collapses in precisely the way INT-0 established for a kinetic
     impactor. The ground laser escapes the transportation cost only because
     Earth's rotation carries targets past it; an orbiting laser has no such
     free ride.

  3. Does it escape the ECONOMICS?  NO -- for the same reason FG01 does not.
     A co-orbital laser must visit its targets, so it inherits ECO-1's
     Delta-v-per-engagement invariant unchanged: pi*v/(3.5*Omega_dot*T).

The honest conclusion is therefore uncomfortable for FG01 in one respect and
strongly favourable in another, and both are reported.
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from fg01.constants import (R_E, MU, RHO_AL, RHO_W, P_W_RANGE, C_LAUNCH_RANGE,
                            P_ELEC, G0)
from fg01.orbital import mass_sphere, xsec_sphere, dv_shift, v_circ
from fg01.laser import (spot_diameter, impulse_per_pulse, coupling_coefficient,
                        F_TH, F_OPT, CM_MAX)
from fg01.relmotion import nodal_rate
from fg01.io_utils import save_table, save_fig, save_json, banner, CB

H = 900e3
D_TARGET = 0.01
M_T = float(mass_sphere(D_TARGET, RHO_AL))
A_T = float(xsec_sphere(D_TARGET))
DV_REQ = float(dv_shift(H, 200e3))
DAY = 86400.0

BETA_FG = 1.2007
V_REL_FG = 619.0
FG_SHOT_G = 2.14           # g of tungsten per engagement (INT-4, 95% conf)
FG_POINTING_URAD = 876.0   # FG01 tolerance


def run():
    banner("LAS-4  The space-based laser")

    # ---- 1. optics in vacuum ----------------------------------------------
    rows = []
    for D_ap in (0.25, 0.5, 1.0, 2.0):
        for lam in (1.06e-6, 532e-9):
            for R in (100.0, 1e3, 1e4, 1e5, 1e6):
                # vacuum: Strehl 1, only laser beam quality degrades it
                d_spot = float(spot_diameter(R, D_ap, lam, beam_quality=1.3,
                                             strehl=1.0))
                a_spot = np.pi * (d_spot / 2) ** 2
                frac = min(A_T / a_spot, 1.0)
                # Pulse energy to put the TARGET at the optimum coupling
                # fluence.  Fluence is energy per unit SPOT area, so the pulse
                # must be sized to the spot, not to the target: sizing to the
                # target drives fluence below threshold whenever the spot is
                # larger, which is the regime that matters.
                e_pulse = F_OPT * a_spot
                p, F, cm = impulse_per_pulse(e_pulse, d_spot, A_T, 1.0)
                dvp = float(p) / M_T
                n_pulses = DV_REQ / dvp if dvp > 0 else np.inf
                rows.append(dict(
                    aperture_m=D_ap, wavelength_nm=lam * 1e9, range_m=R,
                    spot_diameter_m=d_spot,
                    spot_vs_target=d_spot / D_TARGET,
                    fraction_on_target=frac,
                    pulse_energy_J=e_pulse,
                    dv_per_pulse_m_s=dvp,
                    pulses_per_object=n_pulses,
                    total_energy_J=n_pulses * e_pulse if np.isfinite(n_pulses) else np.inf,
                    avg_power_at_10Hz_W=e_pulse * 10.0,
                    pointing_tolerance_urad=(d_spot / 2) / R * 1e6))
    df = pd.DataFrame(rows)
    save_table(df, "las04_vacuum_optics",
               "Space-based laser: spot size, pulse energy and pointing "
               "tolerance in vacuum. Without atmosphere the beam is "
               "diffraction-limited and the spot falls below the target size at "
               "short range, so essentially all energy lands.",
               tags={"vacuum optics": "DERIVED", "C_m": "SOURCED"})
    sel = df[(df.aperture_m == 1.0) & (df.wavelength_nm == 1060.0)]
    print("  1 m aperture, 1.06 um, in vacuum:")
    print(sel[["range_m", "spot_diameter_m", "spot_vs_target",
               "fraction_on_target", "pulse_energy_J", "pulses_per_object",
               "avg_power_at_10Hz_W", "pointing_tolerance_urad"]].to_string(index=False))

    # the good news for the laser, stated plainly
    close = df[(df.aperture_m == 1.0) & (df.wavelength_nm == 1060.0)
               & (df.range_m == 1e4)].iloc[0]
    print(f"\n  at 10 km a 1 m aperture puts a {close.spot_diameter_m*100:.1f} cm spot "
          f"on a 1 cm target: {close.pulses_per_object:.0f} pulses, "
          f"{close.pulse_energy_J:.1f} J each, {close.avg_power_at_10Hz_W:.0f} W average")

    # ---- 2. delivery geometry: does recoil direction save it? -------------
    # Impulse is along the line of sight. Useful retrograde fraction is
    # -LOS_hat . v_hat_target. Parameterise by where the illuminator sits
    # relative to the target in the target's local frame.
    geo = []
    for lab, along, radial, cross in (
            ("directly ahead (co-orbital, leading)", 1.0, 0.0, 0.0),
            ("ahead and 10 deg above", np.cos(np.radians(10)), np.sin(np.radians(10)), 0.0),
            ("ahead and 45 deg above", np.cos(np.radians(45)), np.sin(np.radians(45)), 0.0),
            ("directly above (higher orbit, same plane)", 0.0, 1.0, 0.0),
            ("directly below", 0.0, -1.0, 0.0),
            ("abeam (cross-track, different plane)", 0.0, 0.0, 1.0),
            ("directly behind (trailing)", -1.0, 0.0, 0.0)):
        # (along, radial, cross) is the illuminator's position RELATIVE TO the
        # target, with +along = the target's velocity direction.
        n = np.linalg.norm([along, radial, cross])
        r_laser_to_target = -np.array([along, radial, cross]) / n
        # the target is pushed along the beam, i.e. directly AWAY from the
        # illuminator, so the push is parallel to r_laser_to_target
        push = r_laser_to_target
        useful = -push[0]      # component opposing +along = retrograde
        geo.append(dict(
            illuminator_position=lab,
            useful_retrograde_fraction=float(useful),
            energy_penalty_for_same_kick=(1.0 / useful if useful > 1e-9 else np.inf),
            direction=("retrograde (lowers orbit)" if useful > 1e-9 else
                       "prograde (RAISES orbit)" if useful < -1e-9 else
                       "no along-track component"),
            usable=bool(useful > 0.1)))
    dfg = pd.DataFrame(geo)
    save_table(dfg, "las04_recoil_geometry",
               "Useful retrograde fraction of ablation recoil as a function of "
               "where the illuminator sits relative to the target. Recoil is "
               "along the line of sight, so only an illuminator AHEAD of the "
               "target along its velocity vector produces a retrograde push.")
    print("\n  where the laser must be to push the target retrograde:")
    print(dfg[["illuminator_position","useful_retrograde_fraction","direction","usable"]].to_string(index=False))

    # ---- 3. does it escape the Delta-v invariant? -------------------------
    def dv_per_engagement(h_km, inc_deg, mission_yr):
        a = R_E + h_km * 1e3
        v = float(v_circ(h_km * 1e3))
        om = abs(float(nodal_rate(a, np.radians(inc_deg))))
        return np.pi * v / (3.5 * om * mission_yr * 365.25 * DAY)

    inv = []
    for lab, coorbital in (("FG01 kinetic kick", True),
                           ("Space laser, co-orbital", True),
                           ("Ground laser", False)):
        dvpe = dv_per_engagement(865, 98.8, 10.0) if coorbital else 0.0
        inv.append(dict(
            architecture=lab, must_visit_targets=coorbital,
            dv_per_engagement_m_s=dvpe,
            engagements_per_2000_m_s=(2000.0 / dvpe) if dvpe > 0 else np.inf,
            dv_for_10000_engagements_km_s=(1e4 * dvpe / 1e3) if dvpe > 0 else 0.0,
            note=("inherits ECO-1 invariant unchanged" if coorbital
                  else "Earth's rotation brings targets to it -- no transport cost")))
    dfi = pd.DataFrame(inv)
    save_table(dfi, "las04_dv_invariant",
               "Whether each architecture must transport itself to its targets. "
               "This, not the beam physics, is what decides small-debris "
               "economics.")
    print("\n  the Delta-v invariant applies to any co-orbital architecture:")
    print(dfi[["architecture", "must_visit_targets", "dv_per_engagement_m_s",
               "dv_for_10000_engagements_km_s"]].to_string(index=False))

    # ---- 4. consumables: the axis where the laser genuinely wins ----------
    cons = []
    for n_eng in (24, 100, 1000, 10000):
        # FG01: finite magazine
        fg_mass = n_eng * FG_SHOT_G * 1e-3
        fg_cost = fg_mass * (P_W_RANGE[1] + C_LAUNCH_RANGE[0])
        # space laser: photons from sunlight; only the power system is a cost
        e_per_obj = float(close.total_energy_J)
        las_energy = n_eng * e_per_obj / 0.2       # wall-plug
        cons.append(dict(
            n_engagements=n_eng,
            fg01_tungsten_kg=fg_mass,
            fg01_consumable_usd=fg_cost,
            fg01_magazine_limited=bool(fg_mass > 30.0),
            laser_energy_per_object_J=e_per_obj,
            laser_total_electrical_MJ=las_energy / 1e6,
            laser_solar_charge_time_per_object_s=e_per_obj / 0.2 / 5000.0,
            laser_consumable_usd=0.0))
    dfc = pd.DataFrame(cons)
    save_table(dfc, "las04_consumables",
               "Consumable comparison. A laser's magazine is effectively "
               "infinite: it converts sunlight. FG01 carries finite mass and "
               "must be resupplied.")
    print("\n  consumables -- where the space laser genuinely beats FG01:")
    print(dfc[["n_engagements", "fg01_tungsten_kg", "fg01_consumable_usd",
               "laser_solar_charge_time_per_object_s",
               "laser_consumable_usd"]].to_string(index=False))

    # ---- 5. what FG01 still holds: pointing -------------------------------
    pt = []
    for R in (100.0, 1e3, 1e4, 1e5):
        s = df[(df.aperture_m == 1.0) & (df.wavelength_nm == 1060.0)
               & (df.range_m == R)].iloc[0]
        pt.append(dict(
            engagement_range_m=R,
            space_laser_pointing_urad=float(s.pointing_tolerance_urad),
            fg01_pointing_urad=FG_POINTING_URAD,
            fg01_advantage=FG_POINTING_URAD / float(s.pointing_tolerance_urad),
            laser_spot_cm=float(s.spot_diameter_m) * 100))
    dfp = pd.DataFrame(pt)
    save_table(dfp, "las04_pointing",
               "Pointing tolerance of a space-based laser against FG01's, at "
               "matched engagement range. FG01 keeps a large advantage because "
               "its 'spot' is a physical cloud sized to the miss distance, "
               "while the laser's is set by diffraction.")
    print("\n  pointing, space laser vs FG01:")
    print(dfp.to_string(index=False))

    concl = {
        "question_1_optics": {
            "verdict": "the space laser WINS this outright",
            "detail": (
                f"In vacuum a 1 m aperture at 10 km puts a "
                f"{close.spot_diameter_m*100:.1f} cm spot on a 1 cm target, "
                f"needing {close.pulse_energy_J:.1f} J per pulse and "
                f"{close.avg_power_at_10Hz_W:.0f} W of average power -- against "
                "9.3 kJ and 93 kW for the ground version. The 46x technology "
                "gap of LAS-3 disappears entirely: this is a space-qualifiable "
                "laser, not a national-facility one. The elevation tension, "
                "airmass, turbulence, adaptive optics and terminator window all "
                "vanish with it. Any argument for FG01 that rests on the "
                "ground laser's optics problem does NOT survive this."),
        },
        "question_2_geometry": {
            "verdict": "the space laser does NOT escape",
            "detail": (
                "Ablation recoil is directed along the line of sight, so a "
                "retrograde push requires the illuminator to lie ahead of the "
                "target along its own velocity vector. An illuminator directly "
                "above, below or abeam produces a purely radial or cross-track "
                "impulse with ZERO useful retrograde component. The efficient "
                "geometry for a laser is therefore the same co-orbital leading "
                "geometry INT-0 derived for a kinetic impactor -- the two "
                "methods are constrained identically, because the constraint "
                "comes from the direction the momentum must point, not from how "
                "it is delivered."),
        },
        "question_3_economics": {
            "verdict": "the space laser does NOT escape",
            "dv_per_engagement_m_s": float(dv_per_engagement(865, 98.8, 10.0)),
            "detail": (
                "Because it must be co-orbital to push retrograde, it must "
                "visit its targets, so ECO-1's invariant applies unchanged: "
                f"{dv_per_engagement(865, 98.8, 10.0):.0f} m/s per engagement, "
                "1,070 km/s for 10,000 engagements. A space laser services the "
                "same ~24 objects per platform-decade as FG01. The ground "
                "laser is the only architecture that escapes the invariant, and "
                "it does so only because Earth's rotation carries targets past "
                "it for free -- which is also why it is stuck with the "
                "elevation tension and the 46x technology gap. There is no "
                "configuration that escapes both."),
        },
        "where_the_space_laser_beats_FG01": (
            "Consumables, decisively. A laser converts sunlight and has an "
            "effectively infinite magazine; FG01 carries tungsten and runs out. "
            "At 10,000 engagements FG01 needs 21.4 kg of tungsten launched to "
            "the operating orbit, while the laser needs "
            f"{close.total_energy_J/0.2/5000.0:.1f} s of solar charging per "
            "object and nothing else. If the Delta-v problem were ever solved, "
            "the space laser would be the better architecture and FG01's "
            "magazine would become its binding limit."),
        "where_FG01_still_wins": (
            "Pointing, by roughly two orders of magnitude at matched range, and "
            "technology readiness. FG01's effective 'spot' is a physical cloud "
            "that can be sized to the expected miss distance; the laser's is "
            "set by diffraction and cannot be widened without losing fluence. "
            "FG01 also keeps the graceful-degradation property and needs one "
            "aiming solution rather than tens of consecutive ones."),
        "the_unifying_result": (
            "Momentum must point retrograde. Anything that delivers retrograde "
            "momentum must approach from ahead along the velocity vector, which "
            "forces co-orbital geometry, which forces the Delta-v invariant. "
            "Ablation, kinetic impact and a tug are all subject to it. The only "
            "escape is to stay on the ground and let Earth's rotation do the "
            "transporting -- and that escape costs a 46x technology gap and an "
            "unavoidable geometry penalty. This is a statement about "
            "small-debris remediation in general, not about FG01."),
        "effect_on_the_FG01_accuracy_argument": (
            "It narrows it. Against a GROUND laser the accuracy case is strong "
            "on every axis. Against a SPACE laser only the pointing and "
            "readiness axes survive, and the consumables axis reverses. The "
            "argument should therefore be made against the ground laser, which "
            "is the system the published $100-500/kg figure actually describes "
            "-- and it should concede the space-laser case explicitly rather "
            "than be caught by it."),
    }
    save_json(concl, "las04_conclusions")
    print(f"\n  VERDICT: space laser wins the optics outright; escapes neither "
          f"the recoil geometry nor the Î”v invariant.")

    # ---- figures -----------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.3))
    ax = axes[0]
    for i, D_ap in enumerate((0.25, 0.5, 1.0, 2.0)):
        s = df[(df.aperture_m == D_ap) & (df.wavelength_nm == 1060.0)].sort_values("range_m")
        ax.loglog(s.range_m / 1e3, s.spot_diameter_m * 100, color=CB[i],
                  marker="o", ms=4, label=f"{D_ap:g} m aperture")
    ax.axhline(D_TARGET * 100, color="k", ls="--", lw=1.4)
    ax.text(0.12, 1.25, "1 cm target", fontsize=7.5)
    ax.axhline(50.5, color=CB[5], ls=":", lw=1.4)
    ax.text(0.12, 57, "ground laser spot (50 cm)", fontsize=7.5, color=CB[5])
    ax.set_xlabel("engagement range (km)")
    ax.set_ylabel("spot diameter (cm)")
    ax.set_title("In vacuum the spot can be smaller than the target")
    ax.legend(fontsize=7.5)
    ax = axes[1]
    labels = ["FG01\nkinetic", "Space laser\n(co-orbital)", "Ground laser"]
    vals = [dv_per_engagement(865, 98.8, 10.0), dv_per_engagement(865, 98.8, 10.0), 0.0]
    ax.bar(labels, vals, color=[CB[0], CB[1], CB[2]], width=0.55)
    ax.set_ylabel("Î”v per engagement (m/s)")
    ax.set_title("Who has to visit their targets")
    ax.text(2, 4, "0 â€” Earth's rotation\ndelivers the targets", ha="center",
            fontsize=8, color=CB[2])
    ax.text(0.5, 112, "both inherit ECO-1's invariant â†’ ~24 objects/decade",
            ha="center", fontsize=8)
    fig.suptitle("LAS-4: an orbiting laser fixes the optics and inherits the "
                 "economics", fontsize=10)
    save_fig(fig, "las04_spaceborne",
             "Left: in vacuum a modest aperture drives the spot below the "
             "target size at short range, eliminating the ground laser's "
             "energy and technology problem. Right: but any architecture that "
             "must push retrograde must be co-orbital, and therefore pays the "
             "same Î”v per engagement â€” only the ground station escapes, by "
             "letting Earth's rotation do the transporting.")
    return df, dfg, concl


if __name__ == "__main__":
    run()
