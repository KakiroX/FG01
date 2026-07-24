"""
LAS-3 -- Technology readiness, corrected economics, and the honest verdict.

v6 concluded that FG01 "loses to a ground-based laser on price" by five orders
of magnitude, using a published $100-500/kg that was never re-derived. LAS-1 and
LAS-2 built the laser engagement from physics. This module closes the loop by
asking the two questions the published figure hides:

  1. Does the required laser exist? (SIM-6 asked exactly this of FG01's
     launcher and found it 64,000x inside demonstrated hardware.)
  2. What does the facility cost, and how many objects can it actually service
     given twilight, weather, and the catalogue it depends on?

The answer to (1) is the strongest single result in the laser comparison, and it
is the mirror image of SIM-6.
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from fg01.constants import (R_E, RHO_AL, P_ELEC, C_LAUNCH_RANGE, P_W_RANGE,
                            EML_BENCHMARKS)
from fg01.orbital import mass_sphere, xsec_sphere, dv_shift
from fg01.laser import (spot_diameter, atmospheric_transmission, fried_parameter,
                        ao_strehl, impulse_per_pulse, pass_geometry, F_TH)
from fg01.population import MASTER_TOTAL_GE_1CM, MASTER_TRACKED
from fg01.io_utils import save_table, save_fig, save_json, banner, CB

H = 900e3
D_TARGET = 0.01
M_T = float(mass_sphere(D_TARGET, RHO_AL))
A_T = float(xsec_sphere(D_TARGET))
DV_REQ = float(dv_shift(H, 200e3))
LAM = 532e-9
BEAM_Q = 1.3
PRF = 10.0

# High-energy repetitive nanosecond lasers.  Average power is the binding
# figure of merit for a debris-removal duty cycle.  AM-9: only DEMONSTRATED
# systems set the benchmark; development-stage concepts are carried separately
# as a labelled sensitivity.
LASER_BENCHMARKS = [
    ("LLNL Mercury", 100.0, 10.0, True, "diode-pumped ns, demonstrated"),
    ("DiPOLE-100 (STFC)", 100.0, 10.0, True, "diode-pumped ns, demonstrated"),
    ("HAPLS / ELI-Beamlines", 200.0, 10.0, True, "demonstrated"),
    ("kJ-class repetitive concept", 1000.0, 10.0, False, "development only"),
]

from fg01.laser import F_OPT


def required_pulse_energy(D, target_fluence=F_OPT):
    """
    Pulse energy to put the target at `target_fluence`.

    Sizing to the bare ablation threshold is wrong: the coupling coefficient is
    identically zero there.  A real design targets the optimum-coupling fluence
    F_opt, where momentum per joule is maximised.
    """
    R, elev, retro = pass_geometry(H, np.radians(-7.0))
    T = float(atmospheric_transmission(elev))
    r0 = float(fried_parameter(elev, lam=LAM))
    S = float(ao_strehl(D, r0))
    d_spot = float(spot_diameter(R, D, LAM, BEAM_Q, S))
    area = np.pi * (d_spot / 2) ** 2
    return target_fluence * area / T, d_spot, S, R, retro, T


def run():
    banner("LAS-3  Technology readiness, corrected economics, honest verdict")

    # ---- 1. does the required laser exist? --------------------------------
    demo_max_avg = max(e * f for _, e, f, dem, _ in LASER_BENCHMARKS if dem)
    dev_max_avg = max(e * f for _, e, f, _, _ in LASER_BENCHMARKS)
    save_table(pd.DataFrame([dict(system=n, pulse_J=e, prf_hz=f,
                                  average_power_kW=e * f / 1e3,
                                  demonstrated=dem, note=note)
                             for n, e, f, dem, note in LASER_BENCHMARKS]),
               "las03_laser_benchmarks",
               "High-energy repetitive nanosecond lasers. Only demonstrated "
               "systems set the benchmark (AM-9).")
    rows = []
    for D in (3.0, 5.0, 8.0, 10.0, 15.0, 20.0):
        E_req, d_spot, S, R, retro, T = required_pulse_energy(D)
        # delivered dv per pulse at that energy
        p, F, cm = impulse_per_pulse(E_req, d_spot, A_T, T)
        dvp = float(p) * retro / M_T
        n_pulses = DV_REQ / dvp if dvp > 0 else np.inf
        avg_power = E_req * PRF
        rows.append(dict(
            aperture_m=D, ao_strehl=S, spot_diameter_m=d_spot,
            pulse_energy_required_J=E_req,
            average_power_required_kW=avg_power / 1e3,
            demonstrated_avg_power_kW=demo_max_avg / 1e3,
            gap_beyond_demonstrated=avg_power / demo_max_avg,
            gap_beyond_development_concepts=avg_power / dev_max_avg,
            fluence_J_cm2=F_OPT / 1e4,
            pulses_per_object=n_pulses,
            engagement_time_s=n_pulses / PRF,
            laser_energy_per_object_MJ=n_pulses * E_req / 1e6))
    df = pd.DataFrame(rows)
    save_table(df, "las03_technology_readiness",
               "Laser aperture and average power required to clear the "
               "ablation threshold on a 1 cm target at 900 km, against the "
               "highest demonstrated average power for repetitive nanosecond "
               "lasers.",
               tags={"benchmarks": "SOURCED (demonstrated systems)",
                     "requirement": "DERIVED (LAS-1)"})
    print(df[["aperture_m", "spot_diameter_m", "pulse_energy_required_J",
              "average_power_required_kW", "gap_beyond_demonstrated",
              "pulses_per_object"]].to_string(index=False))

    best = df.loc[df.gap_beyond_demonstrated.idxmin()]
    # FG01's own technology gap, from SIM-6
    fg_muzzle_J = 0.5 * 2.14e-3 * 619.0 ** 2
    fg_demo_J = max(b[3] for b in EML_BENCHMARKS)
    print(f"\n  best case for the laser: {best.aperture_m:.0f} m aperture -> "
          f"{best.average_power_required_kW:.0f} kW average, "
          f"{best.gap_beyond_demonstrated:.0f}x beyond demonstrated")
    print(f"  FG01 launcher: {fg_muzzle_J:.0f} J muzzle vs {fg_demo_J/1e6:.0f} MJ "
          f"demonstrated -> {fg_demo_J/fg_muzzle_J:,.0f}x INSIDE demonstrated")

    # ---- 2. facility cost and throughput ----------------------------------
    thr = []
    for D, c_tel in ((5.0, 150e6), (10.0, 600e6), (15.0, 1500e6)):
        E_req, d_spot, S, R, retro, T = required_pulse_energy(D)
        p, F, cm = impulse_per_pulse(E_req, d_spot, A_T, T)
        dvp = float(p) * retro / M_T
        n_pulses = DV_REQ / dvp if dvp > 0 else np.inf
        t_engage = n_pulses / PRF
        # laser cost scales with average power; ~$50/W for high-energy pulsed
        c_laser = E_req * PRF * 50.0
        c_fac = c_tel + c_laser
        for twilight_hr, weather in ((3.0, 0.6), (4.0, 0.7)):
            # acquisition + track setup per object
            t_overhead = 60.0
            per_night_s = twilight_hr * 3600.0 * weather
            n_per_night = per_night_s / (t_engage + t_overhead)
            n_decade = n_per_night * 365.25 * 10
            thr.append(dict(
                aperture_m=D, facility_cost_musd=c_fac / 1e6,
                telescope_musd=c_tel / 1e6, laser_musd=c_laser / 1e6,
                engagement_time_s=t_engage,
                twilight_hr=twilight_hr, weather_frac=weather,
                objects_per_night=n_per_night,
                objects_per_decade=n_decade,
                cost_per_object_usd=c_fac / max(n_decade, 1e-9),
                electricity_per_object_usd=(n_pulses * E_req / 0.2) / 3.6e6 * P_ELEC))
    dft = pd.DataFrame(thr)
    save_table(dft, "las03_facility_throughput",
               "Ground-laser facility cost and throughput. Throughput is "
               "limited by the terminator window (the target must be sunlit "
               "while the site is dark for the tracking camera to work), by "
               "weather, and by acquisition overhead.",
               tags={"facility cost": "UNVALIDATED (parametric)",
                     "twilight/weather": "SOURCED (site statistics)"})
    print("\n  facility cost and throughput:")
    print(dft[["aperture_m", "facility_cost_musd", "engagement_time_s",
               "objects_per_night", "objects_per_decade",
               "cost_per_object_usd"]].to_string(index=False))

    # ---- 3. the catalogue constraint --------------------------------------
    cat = pd.DataFrame([
        dict(system="Ground laser", needs_position_to_m=0.2525,
             at_range_km=1186.0,
             population_tracked_at_all=MASTER_TRACKED,
             population_ge_1cm=MASTER_TOTAL_GE_1CM,
             frac_tracked=MASTER_TRACKED / MASTER_TOTAL_GE_1CM,
             can_self_acquire=False,
             note="Must be handed a 0.25 m solution; cannot improve on it. "
                  "No catalogue provides this for 1 cm objects."),
        dict(system="FG01", needs_position_to_m=1000.0, at_range_km=50.0,
             population_tracked_at_all=MASTER_TRACKED,
             population_ge_1cm=MASTER_TOTAL_GE_1CM,
             frac_tracked=MASTER_TRACKED / MASTER_TOTAL_GE_1CM,
             can_self_acquire=True,
             note="Needs km-level cueing only; its own sensor closes the loop "
                  "from 100 m inward at 10^11 photoelectrons/s."),
    ])
    save_table(cat, "las03_catalogue_dependency",
               "Both systems depend on an external catalogue that does not yet "
               "exist for the 1 cm population. They depend on it to very "
               "different degrees.")
    print("\n  catalogue dependency:")
    print(cat[["system", "needs_position_to_m", "can_self_acquire"]].to_string(index=False))

    # ---- 4. the corrected comparison --------------------------------------
    las_best = dft.loc[dft.cost_per_object_usd.idxmin()]
    fg01_cost_lo, fg01_cost_hi = 9.09e5, 4.17e6
    comp = pd.DataFrame([
        dict(metric="Cost per 1 cm object removed",
             fg01="$0.9M – 4.2M",
             laser=f"${dft.cost_per_object_usd.min():,.0f} – "
                   f"{dft.cost_per_object_usd.max():,.0f}",
             favours=f"laser ({fg01_cost_lo/dft.cost_per_object_usd.max():,.0f}–"
                     f"{fg01_cost_hi/dft.cost_per_object_usd.min():,.0f}×)"),
        dict(metric="Published figure this study replaces",
             fg01="—", laser="$0.14–0.71 per object ($100–500/kg)",
             favours="—"),
        dict(metric="Correction to the published laser figure",
             fg01="—",
             laser=f"{dft.cost_per_object_usd.min()/0.71:,.0f}–"
                   f"{dft.cost_per_object_usd.max()/0.71:,.0f}× more expensive",
             favours="—"),
        dict(metric="Pointing requirement",
             fg01="876 µrad", laser="0.21 µrad", favours="FG01 (4,115×)"),
        dict(metric="Margin over own diffraction limit",
             fg01="131×", laser="1.64×", favours="FG01 (80×)"),
        dict(metric="Target at the sensor",
             fg01="resolved 149×, 2×10¹¹ e⁻/s",
             laser="unresolved 16×, 4×10⁴ e⁻/s", favours="FG01"),
        dict(metric="Aiming solutions per object",
             fg01="1", laser=f"{best.pulses_per_object:,.0f} consecutive",
             favours=f"FG01 ({best.pulses_per_object:,.0f}×)"),
        dict(metric="External cueing accuracy needed",
             fg01="~1 km", laser="0.25 m", favours="FG01 (3,960×)"),
        dict(metric="Failure mode on aiming error",
             fg01="graceful (45% delivered at 2× error)",
             laser="threshold cliff (zero beyond 1.01×)", favours="FG01"),
        dict(metric="Self-perturbation of the target",
             fg01="none (single shot)",
             laser="loses track in 10 pulses without feedback modelling",
             favours="FG01"),
        dict(metric="Prime hardware vs. demonstrated",
             fg01=f"{fg_demo_J/fg_muzzle_J:,.0f}× inside",
             laser=f"{best.gap_beyond_demonstrated:.0f}× beyond",
             favours="FG01"),
        dict(metric="Atmosphere in the path",
             fg01="none", laser="turbulence, extinction, thermal blooming",
             favours="FG01"),
        dict(metric="Duty cycle",
             fg01="continuous", laser="terminator window only, weather-limited",
             favours="FG01"),
    ])
    save_table(comp, "las03_head_to_head",
               "Corrected head-to-head. The laser retains a large cost "
               "advantage; every accuracy and readiness metric favours FG01.")
    print("\n  head-to-head:")
    print(comp.to_string(index=False))

    concl = {
        "technology_readiness": {
            "laser_best_case": {
                "aperture_m": float(best.aperture_m),
                "average_power_kW": float(best.average_power_required_kW),
                "demonstrated_kW": float(demo_max_avg / 1e3),
                "gap": float(best.gap_beyond_demonstrated)},
            "fg01_launcher_gap": float(fg_demo_J / fg_muzzle_J),
            "statement": (
                f"The laser needs {best.average_power_required_kW:.0f} kW of "
                "average power in nanosecond pulses on a "
                f"{best.aperture_m:.0f} m adaptive-optics aperture. The highest "
                f"demonstrated repetitive nanosecond laser is "
                f"{demo_max_avg/1e3:.0f} kW, so the requirement is "
                f"{best.gap_beyond_demonstrated:.0f}x beyond the state of the "
                "art -- and that is the BEST case across apertures; a 5 m "
                f"aperture needs {float(df[df.aperture_m==5].gap_beyond_demonstrated.iloc[0]):.0f}x. "
                "FG01's launcher, by the identical test in SIM-6, sits "
                f"{fg_demo_J/fg_muzzle_J:,.0f}x INSIDE demonstrated hardware. "
                "This is the single sharpest contrast between the two concepts, "
                "and it runs the opposite way to the cost comparison."),
        },
        "corrected_laser_cost": {
            "cost_per_object_usd": float(las_best.cost_per_object_usd),
            "published_per_object_usd": [0.14, 0.71],
            "correction_factor": float(las_best.cost_per_object_usd / 0.71),
            "statement": (
                f"Built from physics, a ground laser costs about "
                f"${las_best.cost_per_object_usd:,.0f} per 1 cm object -- "
                f"{las_best.cost_per_object_usd/0.71:,.0f}x the $100-500/kg "
                "figure v6 compared against. That published figure is a "
                "marginal-energy cost: the electricity really is only "
                f"${float(las_best.electricity_per_object_usd):.2f} per object. "
                "It omits the facility, whose cost is set by the aperture and "
                "average power the accuracy requirement forces. Correcting it "
                "closes most of the gap v6 reported."),
        },
        "cost_gap_after_correction": {
            "laser_usd_per_object_range": [float(dft.cost_per_object_usd.min()),
                                           float(dft.cost_per_object_usd.max())],
            "fg01_usd_per_object_range": [fg01_cost_lo, fg01_cost_hi],
            "gap_low": float(fg01_cost_lo / dft.cost_per_object_usd.max()),
            "gap_high": float(fg01_cost_hi / dft.cost_per_object_usd.min()),
            "v6_reported_gap": 1.3e6},
        "verdict": (
            f"The laser still wins on cost, by "
            f"{fg01_cost_lo/dft.cost_per_object_usd.max():,.0f}x to "
            f"{fg01_cost_hi/dft.cost_per_object_usd.min():,.0f}x "
            f"(${dft.cost_per_object_usd.min():,.0f}-"
            f"{dft.cost_per_object_usd.max():,.0f} per object against FG01's "
            "$0.9M-4.2M) -- a large gap, but two to four orders of magnitude "
            "smaller than the six-order-of-magnitude gap v6 reported against "
            "the published figure. And it wins that "
            "comparison only by assuming a laser that does not exist: "
            f"{best.gap_beyond_demonstrated:.0f}x beyond demonstrated average "
            "power, pointing at 1.6x its own diffraction limit, tracking an "
            "unresolved point source through turbulence while modelling its own "
            "perturbation of the target. FG01 needs none of that -- its "
            "launcher is 64,000x inside demonstrated hardware and its pointing "
            "requirement has 131x margin. The honest statement is not 'FG01 is "
            "cheaper' but 'FG01 is buildable now and the laser is not, and the "
            "laser's cost advantage is contingent on solving an accuracy "
            "problem that has no demonstrated solution.'"),
        "what_this_does_not_fix": (
            "None of this repairs FG01's own binding constraint. ECO-1's "
            "Delta-v-per-engagement invariant is untouched by the laser "
            "comparison: FG01 still services only ~24 objects per platform "
            "decade, and that -- not the laser benchmark -- is what makes "
            "small-debris removal uneconomic for it. Being more buildable than "
            "a competitor that also does not work is a weak claim, and it "
            "should be made as a claim about RISK and TIMELINE, not about "
            "cost."),
        "where_the_argument_is_strongest": (
            "The accuracy case is strongest not for 1 cm debris but for "
            "intermediate objects, where FG01's areal penalty has already "
            "collapsed to 1 while the laser's problem gets strictly worse: "
            "fluence is energy per unit AREA, so a bigger target needs "
            "proportionally more total impulse with no relief on the fluence "
            "requirement, whereas FG01's mass penalty falls as 1/A_target and "
            "bottoms out at 3 cm. Combined with ECO-5's finding that the "
            "momentum kick is only cost-competitive against a tug between "
            "~100 kg and ~3.7 t, the accuracy argument and the economic "
            "argument point at the SAME target class: tens to hundreds of "
            "kilograms, roughly 30-70 cm. That is the only regime where FG01 "
            "is simultaneously more buildable than a laser and cheaper than a "
            "tug."),
    }
    save_json(concl, "las03_conclusions")

    # ---- 5. how the comparison scales with target size --------------------
    sz = []
    for d_t in (0.01, 0.03, 0.10, 0.30, 1.0):
        a_t = float(xsec_sphere(d_t))
        m_t = float(mass_sphere(d_t, RHO_AL))
        E_req, d_spot, S, R, retro, T = required_pulse_energy(10.0)
        p, F, cm = impulse_per_pulse(E_req, d_spot, a_t, T)
        dvp = float(p) * retro / m_t
        n_pulses = DV_REQ / dvp if dvp > 0 else np.inf
        # FG01 areal penalty at the design miss
        pen = max(np.pi * np.e * 0.005775 ** 2 / a_t, 1.0)
        sz.append(dict(
            target_diameter_cm=d_t * 100, target_mass_g=m_t * 1e3,
            laser_pulses_per_object=n_pulses,
            laser_engagement_time_s=n_pulses / PRF,
            laser_beam_fraction_on_target=min(a_t / (np.pi * (d_spot / 2) ** 2), 1.0),
            fg01_areal_penalty=pen,
            fg01_launched_mass_g=DV_REQ * m_t / (1.2007 * 619.0) * pen * 1e3))
    dfs = pd.DataFrame(sz)
    save_table(dfs, "las03_target_size_scaling",
               "How the comparison scales with target size. FG01's mass penalty "
               "falls as 1/A_target and reaches 1; the laser's engagement time "
               "grows because a larger target needs proportionally more total "
               "impulse while the fluence per unit area is unchanged.")
    print("\n  scaling with target size:")
    print(dfs.to_string(index=False))

    # ---- figures -----------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.3))
    ax = axes[0]
    ax.semilogy(df.aperture_m, df.average_power_required_kW, color=CB[1],
                marker="o", ms=5, lw=2, label="laser average power required")
    ax.axhline(demo_max_avg / 1e3, color="k", ls="--", lw=1.5,
               label=f"highest demonstrated ({demo_max_avg/1e3:.0f} kW)")
    ax.fill_between(df.aperture_m, demo_max_avg / 1e3,
                    df.average_power_required_kW, color=CB[1], alpha=0.12)
    ax.set_xlabel("telescope aperture (m)")
    ax.set_ylabel("average laser power required (kW)")
    ax.set_title("The laser is beyond the state of the art")
    ax.legend(fontsize=7.5)
    ax = axes[1]
    labels = ["FG01\nlauncher", "Ground laser\n(10 m aperture)"]
    vals = [fg_demo_J / fg_muzzle_J, 1.0 / best.gap_beyond_demonstrated]
    colors = [CB[2], CB[1]]
    ax.bar(labels, vals, color=colors, width=0.55)
    ax.axhline(1.0, color="k", lw=1.4)
    ax.set_yscale("log")
    ax.set_ylabel("demonstrated capability ÷ requirement")
    ax.text(0, vals[0] * 2, f"{vals[0]:,.0f}× inside", ha="center", fontsize=9)
    ax.text(1, vals[1] * 2, f"{1/vals[1]:.0f}× beyond", ha="center", fontsize=9)
    ax.text(1.45, 1.3, "buildable", fontsize=7.5, color="0.3")
    ax.set_title("Prime hardware vs. what exists")
    fig.suptitle("LAS-3: the laser's cost advantage is contingent on hardware "
                 "that does not exist", fontsize=10)
    save_fig(fig, "las03_readiness",
             "Left: average laser power required to clear the ablation "
             "threshold, against the highest demonstrated repetitive "
             "nanosecond laser. Right: each concept's prime hardware measured "
             "against demonstrated capability — FG01's coilgun sits far inside "
             "the envelope, the laser far outside it.")
    return df, comp, concl


if __name__ == "__main__":
    run()
