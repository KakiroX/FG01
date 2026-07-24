"""
SIM-14 support -- space-based laser-ablation error budget, for direct
comparison against the FG01 terminal-GNC budget in gnc.py.

Source for every laser-side number: `laser accuracy comm.md` (an informal
brief summarising secondary/non-peer-reviewed material, not independently
verified in this study -- see docs/CITATION_LOG.md, key `laser_brief`). Three
things are taken directly from that brief and used as stated, not tuned:

  1. Orbit-prediction requirement "better than 1 m", with an optimistic claim
     of "<1 cm within 10 s of detection".
  2. Pointing requirement described only as "sub-arcsecond" (no exact figure
     given in the brief; a 0.1-2.0 arcsec envelope is used here and tagged
     UNVALIDATED/engineering estimate).
  3. Post-Newtonian (p-N) correction: the brief states that omitting it
     produces an aiming bias "on the order of the size of the debris objects
     themselves". This is implemented literally, as an additive residual
     equal to the target radius when p-N correction is switched off.

Engagement range is NOT given in the brief (it discusses angular/temporal
accuracy only). A representative range grid for space-based laser ablation
concepts is used here (100-1000 km) and tagged as this study's own
assumption, not sourced from the brief.

The laser's failure mode is modelled as the brief's own comparison table
states it: a miss larger than the target radius delivers no ablation
impulse at all (binary), in contrast to FG01's graceful areal-capture model
(delivery.py). This is the central structural asymmetry the brief itself
identifies as the "Achilles' heel" of the laser approach.
"""
import numpy as np

ARCSEC_TO_RAD = np.pi / (180.0 * 3600.0)


def pointing_linear_error(pointing_arcsec, R_m):
    """Linear pointing error at range R from an angular pointing budget."""
    return np.asarray(pointing_arcsec, float) * ARCSEC_TO_RAD * np.asarray(R_m, float)


def laser_miss_budget(sigma_orbit_pred_m, pointing_arcsec, R_m, pN_residual_m=0.0):
    """
    Total 1-sigma miss budget for a space-based laser-ablation shot.

        sigma_total = sqrt(sigma_orbit_pred^2 + (pointing_arcsec * R)^2 + pN_residual^2)

    All three terms are taken from / motivated directly by `laser accuracy
    comm.md`: orbit-prediction error, angular pointing error scaled by range,
    and the p-N correction residual (zero if corrected, ~target radius if
    not -- the brief's own stated consequence of skipping it).
    """
    point = pointing_linear_error(pointing_arcsec, R_m)
    terms = dict(
        orbit_pred=np.asarray(sigma_orbit_pred_m, float),
        pointing=point,
        pN_residual=np.asarray(pN_residual_m, float),
    )
    total = np.sqrt(sum(v ** 2 for v in terms.values()))
    return total, terms


def p_effective_shot(sigma_total, r_target):
    """
    Probability that a single shot delivers any ablation effect at all, i.e.
    that a Rayleigh-distributed 2-D miss (per-axis std sigma_total) falls
    within the target radius. This is the brief's own binary failure mode:
    "Miss entirely; no effect."
    """
    sigma_total = np.asarray(sigma_total, float)
    r_target = np.asarray(r_target, float)
    return 1.0 - np.exp(-r_target ** 2 / (2.0 * sigma_total ** 2))
