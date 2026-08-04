"""
COIL-03 -- Full radial stack-up of the coilgun tube, bore to outside.

STANDALONE AND ADDITIVE. Imports fg01/ read-only, writes only coil03_* outputs,
not registered in run_all.py.

COIL-02 sized the structural walls but treated the winding as a thin shell. The
winding is in fact the THICKEST layer in the stack, because a 40 mm long coil
cannot fit 22 turns of AWG 10 in a single layer -- it needs two, and each layer
is a full wire diameter of radial build.

Stack, inside to outside:

    canister (30 mm OD)
    | running clearance
    | BORE TUBE          insulating, structural-ish
    | former / ground insulation
    | WINDING            n_layers x wire diameter   <- dominant
    | ground insulation
    | BANDING            reacts B^2/2mu0
    | thermal gap
    | OUTER CASE
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd

from fg01.io_utils import save_table, save_json, banner

RHO_CU = 8960.0
BORE_ID = 0.030
N_STAGES = 31
L_COIL = 0.040

# AWG: (diameter m, area m2)
WIRE = {
    "AWG 8  round":  (3.264e-3, 8.37e-6),
    "AWG 10 round":  (2.588e-3, 5.26e-6),
    "AWG 12 round":  (2.053e-3, 3.31e-6),
    "AWG 14 round":  (1.628e-3, 2.08e-6),
    "3x2 mm rect":   (2.000e-3, 6.00e-6),   # radial build 2 mm, 3 mm axial
    "4x1.5 mm rect": (1.500e-3, 6.00e-6),   # radial build 1.5 mm, 4 mm axial
}
# for rectangular, the axial pitch differs from the radial build
AXIAL = {"3x2 mm rect": 3.0e-3, "4x1.5 mm rect": 4.0e-3}

LAYERS = [
    ("running clearance", 0.5e-3, None),
    ("bore tube (G10)", 1.5e-3, 1850.0),
    ("former / ground insulation", 0.5e-3, 1400.0),
    ("WINDING", None, RHO_CU),            # computed
    ("ground insulation", 0.5e-3, 1400.0),
    ("banding (S-glass)", 0.93e-3, 2000.0),
    ("thermal gap", 1.0e-3, None),
    ("outer case (Ti)", 0.8e-3, 4430.0),
]


def winding_build(n_turns, wire, l_coil=L_COIL):
    """Radial build of the winding, and the layer count."""
    d_rad, area = WIRE[wire]
    pitch = AXIAL.get(wire, d_rad)
    per_layer = max(int(np.floor(l_coil / pitch)), 1)
    n_layers = int(np.ceil(n_turns / per_layer))
    return n_layers * d_rad, n_layers, per_layer


def stackup(n_turns, wire):
    r = BORE_ID / 2
    rows = []
    w_build, n_lay, per_lay = winding_build(n_turns, wire)
    for name, t, rho in LAYERS:
        if name == "WINDING":
            t = w_build
            note = f"{n_lay} layers x {per_lay} turns, {wire}"
        else:
            note = ""
        r_in, r_out = r, r + t
        vol = np.pi * (r_out ** 2 - r_in ** 2) * L_COIL if rho else 0.0
        rows.append(dict(layer=name, thickness_mm=t * 1e3,
                         r_inner_mm=r_in * 1e3, r_outer_mm=r_out * 1e3,
                         mass_per_stage_g=(vol * rho * 1e3) if rho else 0.0,
                         note=note))
        r = r_out
    return pd.DataFrame(rows), r


def run():
    banner("COIL-03  Full radial stack-up of the coilgun tube")

    # ---- headline stack --------------------------------------------------
    N, WIRE_SEL = 22, "AWG 10 round"
    df, r_out = stackup(N, WIRE_SEL)
    total_wall = r_out - BORE_ID / 2
    print(f"\n[1] STACK-UP  ({N} turns, {WIRE_SEL}, {L_COIL*1e3:.0f} mm coil)\n")
    print(df.to_string(index=False, float_format=lambda x: f"{x:.2f}"))
    print(f"\n    bore inner radius      {BORE_ID/2*1e3:6.2f} mm")
    print(f"    outer radius           {r_out*1e3:6.2f} mm")
    print(f"    TOTAL RADIAL WALL      {total_wall*1e3:6.2f} mm")
    print(f"    OUTER DIAMETER         {2*r_out*1e3:6.2f} mm  "
          f"(for a {BORE_ID*1e3:.0f} mm bore)")
    w_build, n_lay, _ = winding_build(N, WIRE_SEL)
    print(f"    winding is {100*w_build/total_wall:.0f}% of the wall "
          f"({n_lay} layers)")
    save_table(df, "coil03_stackup",
               f"Radial stack-up, bore to outside, {N} turns of {WIRE_SEL}.")

    # ---- consistency with COIL-01 ---------------------------------------
    r_w_in = BORE_ID / 2 + 0.5e-3 + 1.5e-3 + 0.5e-3
    r_mean = r_w_in + w_build / 2
    print(f"\n[2] CONSISTENCY CHECK vs COIL-01")
    print(f"    winding mean radius from this stack {r_mean*1e3:.2f} mm")
    print(f"    R_COIL assumed in COIL-01           20.00 mm")
    print(f"    agreement                           "
          f"{100*abs(r_mean-0.020)/0.020:.1f}%  -> COIL-01 electromagnetics stand")

    # ---- wire / turns trade ----------------------------------------------
    print(f"\n[3] WIRE AND TURN-COUNT TRADE\n")
    rows = []
    for wire in WIRE:
        for n in (15, 22, 30, 45):
            wb, nl, pl = winding_build(n, wire)
            _, ro = stackup(n, wire)
            cu_vol = n * 2 * np.pi * (BORE_ID / 2 + 2.5e-3 + wb / 2) * WIRE[wire][1]
            rows.append(dict(wire=wire, turns=n, layers=nl, turns_per_layer=pl,
                             winding_mm=wb * 1e3,
                             total_wall_mm=(ro - BORE_ID / 2) * 1e3,
                             OD_mm=2 * ro * 1e3,
                             cu_kg_all_stages=cu_vol * RHO_CU * N_STAGES))
    d2 = pd.DataFrame(rows)
    save_table(d2, "coil03_wire_trade",
               "Winding build, total wall and outer diameter across wire gauge "
               "and turn count.")
    print(d2[d2.turns.isin([22, 30])].to_string(
        index=False, float_format=lambda x: f"{x:.2f}"))

    # ---- envelope --------------------------------------------------------
    print(f"\n[4] LAUNCHER ENVELOPE\n")
    total_mass = float(df.mass_per_stage_g.sum()) * N_STAGES / 1e3
    L_total = 1.364
    print(f"    outer diameter         {2*r_out*1e3:6.1f} mm")
    print(f"    length                 {L_total*1e3:6.0f} mm")
    print(f"    tube mass (all layers) {total_mass:6.2f} kg over {N_STAGES} stages")
    print(f"    swept volume           "
          f"{np.pi*r_out**2*L_total*1e3:6.2f} litres")

    concl = {
        "total_radial_wall_mm": float(total_wall * 1e3),
        "outer_diameter_mm": float(2 * r_out * 1e3),
        "bore_mm": BORE_ID * 1e3,
        "winding_build_mm": float(w_build * 1e3),
        "winding_share_of_wall_pct": float(100 * w_build / total_wall),
        "n_layers": int(n_lay),
        "headline": (
            f"For a {BORE_ID*1e3:.0f} mm bore the full radial build is "
            f"{total_wall*1e3:.1f} mm, giving a {2*r_out*1e3:.0f} mm outer "
            f"diameter. The WINDING dominates at {w_build*1e3:.1f} mm "
            f"({100*w_build/total_wall:.0f}% of the wall) because 22 turns of "
            f"AWG 10 need {n_lay} layers in a 40 mm coil. The structural walls "
            f"COIL-02 sized are only {total_wall*1e3 - w_build*1e3:.1f} mm of it."),
        "for_cost_model": (
            "Envelope: 52 mm OD x 1364 mm. Use this for stowage and for the "
            "outer-case area; the mass model in COIL-02 is unchanged, since "
            "copper mass depends on turns and mean radius, not on layer count."),
        "lever": (
            "To shrink the tube, use finer wire -- but check the cost. AWG 12 "
            "gives 49.7 mm OD (from 51.8) and AWG 14 at 22 turns fits a SINGLE "
            "layer for 44.7 mm OD. However conductor area falls from 5.26 mm^2 "
            "(AWG 10) to 3.31 (AWG 12) to 2.08 (AWG 14), so coil resistance "
            "rises 1.6x and 2.5x respectively and ohmic loss rises with it. "
            "COIL-01's 19.6% efficiency was computed at AWG 10; going finer "
            "requires re-running it. Rectangular conductor is the better lever: "
            "3x2 mm at 22 turns gives 49.5 mm OD at 6.0 mm^2 area, i.e. a "
            "smaller tube AND lower resistance than AWG 10."),
    }
    save_json(concl, "coil03_conclusions")
    return df, d2, concl


if __name__ == "__main__":
    run()
