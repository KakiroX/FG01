"""
COIL-01 -- Multistage induction coilgun for the FG01 3 cm rhenium canister.

Standalone: python coil01_design.py   (numpy / pandas / matplotlib only)

--------------------------------------------------------------------------
WHY AN *INDUCTION* COILGUN AND NOT A RELUCTANCE ONE
--------------------------------------------------------------------------
A reluctance coilgun pulls a FERROMAGNETIC slug into a coil. Rhenium is
paramagnetic (chi ~ +7e-5) -- for practical purposes non-magnetic -- so a
reluctance gun cannot grip the payload at all. It would need an iron armature,
which saturates at ~2.1 T and therefore caps the achievable magnetic pressure
at B^2/2mu0 = 1.75 MPa. That is 6x too low for this design (see below), before
counting the dead mass of the iron.

An INDUCTION launcher drives eddy currents in a CONDUCTING armature and is not
saturation-limited -- its ceiling is set by coil structural strength, not by a
material property. Rhenium's own conductivity is poor (5.42e6 S/m, 11x worse
than copper), so the payload cannot be its own armature either. The design is
therefore: a thin ALUMINIUM canister that both holds the rhenium grains and
acts as the armature, and which flies with the payload (no discarding sabot --
that would create exactly the debris FG01 exists to remove).

MASS BUDGET, given the 100 g total accelerated mass:
    aluminium canister/armature (3 cm dia, 4 cm long, 2 mm wall)  ~19 g
    rhenium payload                                              ~81 g
The armature is 19% of the launch mass and is NOT optional -- it is what the
field pushes on.

--------------------------------------------------------------------------
MODEL
--------------------------------------------------------------------------
Per stage, a capacitor discharges into a drive coil that is magnetically
coupled to the armature. Two coupled circuits plus motion:

    lambda_c = L_c I_c + M(x) I_a          (coil flux linkage)
    lambda_a = M(x) I_c + L_a I_a          (armature flux linkage)

    V_c = R_c I_c + dlambda_c/dt
    0   = R_a I_a + dlambda_a/dt           (armature is a shorted turn)
    m dv/dt = I_c I_a dM/dx                (Lorentz force from co-energy)
    dV_c/dt = -I_c/C

Expanding the flux derivatives gives the 2x2 system solved each step:

    [L_c  M ] [dI_c/dt]   [V_c - R_c I_c - v (dM/dx) I_a]
    [M   L_a] [dI_a/dt] = [    -R_a I_a - v (dM/dx) I_c ]

Mutual inductance uses the coaxial-dipole profile, which captures the essential
behaviour (peaks at alignment, falls off over a length scale ~R):

    M(x) = M_peak / (1 + (x/R)^2)^{3/2},    M_peak = k_max sqrt(L_c L_a)

The armature sits AHEAD of the coil, so dM/dx < 0; Lenz makes I_a opposite in
sign to I_c, and the product I_c I_a dM/dx is positive -- the armature is
repelled forward. This is the standard coaxial induction launcher geometry.

--------------------------------------------------------------------------
ASSUMPTIONS
--------------------------------------------------------------------------
A1  Lumped circuit: coil and armature each a single lumped L-R, no skin-depth
    or proximity-effect current redistribution within a pulse.
A2  Dipole M(x) profile with a fitted peak coupling k_max = 0.6 (a good but
    not heroic value for a coaxial launcher).
A3  Rigid armature, no melting, no deformation under magnetic pressure.
A4  Ideal switching: each stage fires instantly at its trigger position, and
    the capacitor is not crowbarred (current is allowed to ring).
A5  Room-temperature resistivities, constant. Real coil resistance rises with
    temperature during the pulse -- this makes the model slightly optimistic.
A6  No stage-to-stage magnetic interaction.
A7  Vacuum, no air drag.

NON-LINEAR EFFECTS NOT MODELLED: skin effect in the armature (raises effective
R_a and reduces efficiency), armature heating over a burst, and flux diffusion
through the canister wall.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pathlib
import json

OUT = pathlib.Path(__file__).resolve().parent.parent / "outputs"
FIG = pathlib.Path(__file__).resolve().parent.parent / "figures"
OUT.mkdir(exist_ok=True)
FIG.mkdir(exist_ok=True)

CB = ["#0072B2", "#D55E00", "#009E73", "#CC79A7",
      "#E69F00", "#56B4E9", "#F0E442", "#000000"]
plt.rcParams.update({"figure.dpi": 130, "font.size": 9.5, "axes.grid": True,
                     "grid.alpha": 0.25, "axes.spines.top": False,
                     "axes.spines.right": False, "legend.frameon": False,
                     "savefig.bbox": "tight"})

MU0 = 4e-7 * np.pi
G0 = 9.80665
RHO_CU = 1.68e-8            # ohm.m
RHO_AL = 2.65e-8            # ohm.m
RHO_AL_D = 2700.0           # kg/m3
RHO_RE_D = 21010.0          # kg/m3
RHO_CU_D = 8960.0           # kg/m3

# ==========================================================================
# GEOMETRY -- fixed by the 3 cm canister
# ==========================================================================
BORE = 0.030                # m, canister outer diameter
M_PROJ = 0.100              # kg, TOTAL accelerated mass (payload + canister)
T_WALL = 0.002              # m, aluminium canister wall
L_ARM = 0.040               # m, canister length
R_ARM = (BORE - T_WALL) / 2         # armature mean radius, m
R_COIL = BORE / 2 + 0.005           # coil mean radius (5 mm of tube + former)

K_MAX = 0.60                # peak coupling coefficient
CAP_J_PER_KG = 1500.0       # pulse film capacitors, J/kg
ETA_CHARGER = 0.90          # PSU efficiency, capacitor charging


def armature_params():
    """Aluminium canister acting as the shorted-turn armature."""
    a_cross = T_WALL * L_ARM                     # current-carrying cross-section
    R_a = RHO_AL * 2 * np.pi * R_ARM / a_cross
    a_eq = np.sqrt(a_cross / np.pi)              # equivalent conductor radius
    L_a = MU0 * R_ARM * (np.log(8 * R_ARM / a_eq) - 2.0)
    m_arm = 2 * np.pi * R_ARM * T_WALL * L_ARM * RHO_AL_D
    return L_a, R_a, m_arm


def coil_params(N, l_coil, wire_area):
    """Solenoid inductance (Wheeler) and DC resistance."""
    L_c = MU0 * N ** 2 * np.pi * R_COIL ** 2 / (l_coil + 0.9 * R_COIL)
    len_wire = 2 * np.pi * R_COIL * N
    R_c = RHO_CU * len_wire / wire_area
    m_cu = len_wire * wire_area * RHO_CU_D
    return L_c, R_c, m_cu


def M_of_x(x, M_peak, R):
    return M_peak / (1.0 + (x / R) ** 2) ** 1.5


def dM_of_x(x, M_peak, R):
    return -3.0 * M_peak * x / R ** 2 / (1.0 + (x / R) ** 2) ** 2.5


def stage_sim_vec(v_in, C, V0, N, l_coil, wire_area, m, x_trig,
                  dt=1e-7, x_end_factor=5.0, t_max=6e-3):
    """
    Vectorised stage simulation: x_trig (and optionally C, V0) may be arrays,
    and all configurations are integrated simultaneously. This is ~50x faster
    than looping, which matters because the design search is nested three deep.
    """
    L_a, R_a, _ = armature_params()
    L_c, R_c, _ = coil_params(N, l_coil, wire_area)
    M_peak = K_MAX * np.sqrt(L_c * L_a)

    x_trig = np.atleast_1d(np.asarray(x_trig, float))
    n = x_trig.size
    C = np.broadcast_to(np.asarray(C, float), (n,)).astype(float)
    V0 = np.broadcast_to(np.asarray(V0, float), (n,)).astype(float)

    x = x_trig.copy()
    v = np.full(n, float(v_in))
    Ic = np.zeros(n)
    Ia = np.zeros(n)
    Vc = V0.copy()
    F_max = np.zeros(n)
    q_coil = np.zeros(n)
    q_arm = np.zeros(n)
    x_end = x_end_factor * R_COIL
    alive = np.ones(n, dtype=bool)

    def deriv(x, v, Ic, Ia, Vc):
        u = 1.0 + (x / R_COIL) ** 2
        M = M_peak / u ** 1.5
        dM = -3.0 * M_peak * x / R_COIL ** 2 / u ** 2.5
        d = L_c * L_a - M ** 2
        ec = Vc - R_c * Ic - v * dM * Ia
        ea = -R_a * Ia - v * dM * Ic
        dIc = (L_a * ec - M * ea) / d
        dIa = (L_c * ea - M * ec) / d
        F = Ic * Ia * dM
        return v, F / m, dIc, dIa, -Ic / C, F

    steps = int(t_max / dt)
    for _ in range(steps):
        k1 = deriv(x, v, Ic, Ia, Vc)
        F = k1[5]
        F_max = np.where(alive, np.maximum(F_max, np.abs(F)), F_max)
        q_coil += np.where(alive, R_c * Ic ** 2 * dt, 0.0)
        q_arm += np.where(alive, R_a * Ia ** 2 * dt, 0.0)
        alive &= (x <= x_end)
        if not alive.any():
            break
        k2 = deriv(x + .5 * dt * k1[0], v + .5 * dt * k1[1], Ic + .5 * dt * k1[2],
                   Ia + .5 * dt * k1[3], Vc + .5 * dt * k1[4])
        k3 = deriv(x + .5 * dt * k2[0], v + .5 * dt * k2[1], Ic + .5 * dt * k2[2],
                   Ia + .5 * dt * k2[3], Vc + .5 * dt * k2[4])
        k4 = deriv(x + dt * k3[0], v + dt * k3[1], Ic + dt * k3[2],
                   Ia + dt * k3[3], Vc + dt * k3[4])
        a = alive.astype(float)
        x = x + a * (dt / 6) * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
        v = v + a * (dt / 6) * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        Ic = Ic + a * (dt / 6) * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2])
        Ia = Ia + a * (dt / 6) * (k1[3] + 2 * k2[3] + 2 * k3[3] + k4[3])
        Vc = Vc + a * (dt / 6) * (k1[4] + 2 * k2[4] + 2 * k3[4] + k4[4])

    E_cap0 = 0.5 * C * V0 ** 2
    dKE = 0.5 * m * (v ** 2 - float(v_in) ** 2)
    B_peak = np.sqrt(2 * MU0 * np.maximum(F_max, 0) / (np.pi * R_COIL ** 2))
    return dict(v_out=v, dKE=dKE, E_cap0=E_cap0,
                eta_stage=dKE / np.maximum(E_cap0, 1e-12),
                a_max_g=F_max / m / G0, B_peak=B_peak,
                q_coil=q_coil, q_arm=q_arm, F_max=F_max)


def stage_sim(v_in, C, V0, N, l_coil, wire_area, m, x_trig,
              dt=2e-8, x_end_factor=6.0):
    """
    Fire one stage. Armature enters at v_in, positioned x_trig AHEAD of the
    coil centre when the switch closes. Integrates until the armature is far
    enough downstream that coupling is negligible.

    Returns dict with exit velocity, energy drawn, losses, peak force, etc.
    """
    L_a, R_a, _ = armature_params()
    L_c, R_c, _ = coil_params(N, l_coil, wire_area)
    M_peak = K_MAX * np.sqrt(L_c * L_a)
    det = L_c * L_a - M_peak ** 2
    if det <= 0:
        return None

    x, v, Ic, Ia, Vc = x_trig, v_in, 0.0, 0.0, V0
    x_end = x_end_factor * R_COIL
    E_cap0 = 0.5 * C * V0 ** 2
    q_coil = q_arm = 0.0
    F_max = 0.0
    t = 0.0
    hist = []

    def deriv(s):
        x, v, Ic, Ia, Vc = s
        M = M_of_x(x, M_peak, R_COIL)
        dM = dM_of_x(x, M_peak, R_COIL)
        d = L_c * L_a - M ** 2
        ec = Vc - R_c * Ic - v * dM * Ia
        ea = -R_a * Ia - v * dM * Ic
        dIc = (L_a * ec - M * ea) / d
        dIa = (L_c * ea - M * ec) / d
        F = Ic * Ia * dM
        return np.array([v, F / m, dIc, dIa, -Ic / C]), F

    s = np.array([x, v, Ic, Ia, Vc])
    for _ in range(int(4e-3 / dt)):
        k1, F = deriv(s)
        F_max = max(F_max, abs(F))
        q_coil += R_c * s[2] ** 2 * dt
        q_arm += R_a * s[3] ** 2 * dt
        hist.append((t, s[0], s[1], s[2], s[3], s[4], F))
        if s[0] > x_end:
            break
        k2, _ = deriv(s + 0.5 * dt * k1)
        k3, _ = deriv(s + 0.5 * dt * k2)
        k4, _ = deriv(s + dt * k3)
        s = s + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        t += dt

    v_out = s[1]
    E_cap_res = 0.5 * C * s[4] ** 2
    dKE = 0.5 * m * (v_out ** 2 - v_in ** 2)
    E_used = E_cap0 - E_cap_res
    return dict(v_out=v_out, dKE=dKE, E_cap0=E_cap0, E_used=max(E_used, 1e-12),
                eta_stage=dKE / E_cap0, eta_of_used=dKE / max(E_used, 1e-12),
                q_coil=q_coil, q_arm=q_arm, F_max=F_max,
                a_max_g=F_max / m / G0,
                B_peak=np.sqrt(2 * MU0 * F_max / (np.pi * R_COIL ** 2))
                if F_max > 0 else 0.0,
                t_stage=t, hist=np.array(hist))


def best_trigger(v_in, C, V0, N, l_coil, wire_area, m):
    """Scan the trigger offset for maximum energy transfer."""
    best = None
    for xt in np.linspace(0.15, 1.4, 18) * R_COIL:
        r = stage_sim(v_in, C, V0, N, l_coil, wire_area, m, xt,
                      dt=5e-8)
        if r and (best is None or r["dKE"] > best["dKE"]):
            best = r
            best["x_trig"] = xt
    return best


# ==========================================================================
# DESIGN CONSTRAINTS -- the coilgun's analogues of chamber pressure
# ==========================================================================
# 1. Canister structure. The base of the aluminium canister must carry the
#    inertial load of the rhenium column above it:
#       sigma = rho_Re * a * L_column
#    81 g of rhenium in a 3 cm bore is a 7.3 mm column, so at 50,000 g the
#    base sees 75 MPa against Al-6061's 276 MPa yield -- SF 3.7. This is the
#    binding limit on the payload, not on the gun.
# 2. Coil structure. Magnetic pressure B^2/(2 mu0) tries to burst the drive
#    coil exactly as chamber pressure bursts a barrel. At 15 T that is 90 MPa,
#    which bare copper (yield 70-200 MPa) cannot hold alone -- the coil needs
#    steel or composite banding. Above ~20 T this becomes pulsed-magnet
#    engineering and is ruled out here.
A_MAX_G = 50_000.0
B_MAX_T = 15.0

L_RE_COLUMN = None          # computed in main()


def canister_base_stress(a_g, m_payload):
    """Compressive stress at the canister base from the rhenium column [Pa]."""
    a_cross = np.pi * (BORE / 2 - T_WALL) ** 2
    l_col = m_payload / (RHO_RE_D * a_cross)
    return RHO_RE_D * (a_g * G0) * l_col, l_col


def stage_best_config(v_in, N, l_coil, wire_area, m=M_PROJ,
                      a_max_g=A_MAX_G, b_max=B_MAX_T):
    """
    Optimise one stage at entry velocity v_in.

    The pulse must be matched to the transit time through the coupling region
    (~2R/v), so the optimum capacitance falls as 1/v^2 -- later stages want
    shorter, sharper pulses. For each candidate C the voltage is pushed up to
    whichever constraint binds first (acceleration or field), then the C giving
    the most kinetic energy is chosen.
    """
    # OBJECTIVE: maximise stage EFFICIENCY (dKE / E_stored), not dKE. Maximising
    # dKE alone drives the optimiser to enormous low-voltage capacitors that
    # store hundreds of kJ to deliver a few hundred J -- the right answer to the
    # wrong question. Total capacitor mass is what the spacecraft actually pays.
    L_c, _, _ = coil_params(N, l_coil, wire_area)
    v_ref = max(v_in, 15.0)
    C_match = (2.0 * R_COIL / v_ref) ** 2 / L_c

    Cs = C_match * np.array([0.005, 0.012, 0.03, 0.07, 0.15, 0.35, 0.8, 1.6])
    Cs = np.clip(Cs, 2e-6, 0.05)
    Vs = np.array([1e3, 2e3, 4e3, 7e3, 12e3, 18e3, 26e3, 36e3])
    Ts = np.linspace(0.08, 1.2, 8) * R_COIL

    CC, VV, TT = np.meshgrid(Cs, Vs, Ts, indexing="ij")
    CC, VV, TT = CC.ravel(), VV.ravel(), TT.ravel()
    r = stage_sim_vec(v_in, CC, VV, N, l_coil, wire_area, m, TT)

    ok = (r["a_max_g"] <= a_max_g) & (r["B_peak"] <= b_max) & (r["dKE"] > 0)
    if not ok.any():
        return None
    score = np.where(ok, r["eta_stage"], -np.inf)
    j = int(np.argmax(score))
    return dict(v_out=float(r["v_out"][j]), dKE=float(r["dKE"][j]),
                E_cap0=float(r["E_cap0"][j]),
                eta_stage=float(r["eta_stage"][j]),
                a_max_g=float(r["a_max_g"][j]), B_peak=float(r["B_peak"][j]),
                q_coil=float(r["q_coil"][j]), q_arm=float(r["q_arm"][j]),
                C=float(CC[j]), V=float(VV[j]), x_trig=float(TT[j]))


def multistage(v_target, N, l_coil, wire_area, m=M_PROJ,
               stage_pitch=None, max_stages=200,
               a_max_g=A_MAX_G, b_max=B_MAX_T):
    """Chain optimally-graded stages until v_target is reached."""
    if stage_pitch is None:
        stage_pitch = 2.2 * R_COIL
    v = 0.0
    stages = []
    E_tot = 0.0
    for i in range(max_stages):
        r = stage_best_config(v, N, l_coil, wire_area, m, a_max_g, b_max)
        if r is None or r["dKE"] <= 0:
            break
        v_prev = v
        v = r["v_out"]
        E_tot += r["E_cap0"]
        stages.append(dict(stage=i + 1, v_in=v_prev, v_out=v,
                           C_uF=r["C"] * 1e6, V_kV=r["V"] / 1e3,
                           E_J=r["E_cap0"], dKE=r["dKE"],
                           eta_pct=100 * r["eta_stage"],
                           a_max_g=r["a_max_g"], B_peak=r["B_peak"],
                           q_coil=r["q_coil"], q_arm=r["q_arm"]))
        if v >= v_target:
            break
    if not stages:
        return None
    n = len(stages)
    _, _, m_cu_stage = coil_params(N, l_coil, wire_area)
    return dict(n_stages=n, v_final=v, E_stored_total=E_tot,
                length=n * stage_pitch,
                eta_total=0.5 * m * v ** 2 / E_tot,
                a_max_g=max(s["a_max_g"] for s in stages),
                B_peak=max(s["B_peak"] for s in stages),
                q_total=sum(s["q_coil"] + s["q_arm"] for s in stages),
                m_cu=n * m_cu_stage,
                m_cap=E_tot / CAP_J_PER_KG,
                E_peak_stage=max(s["E_J"] for s in stages),
                stages=stages)


# ==========================================================================
# MAIN
# ==========================================================================
def main():
    print("=" * 76)
    print("  COIL-01  Multistage induction coilgun, 3 cm bore, 100 g total mass")
    print("=" * 76)

    L_a, R_a, m_arm = armature_params()
    s50, l_col = canister_base_stress(A_MAX_G, M_PROJ - m_arm)
    print(f"\n[0] CONFIGURATION")
    print(f"    bore                       {BORE*1e3:.0f} mm")
    print(f"    total accelerated mass     {M_PROJ*1e3:.0f} g")
    print(f"    aluminium canister/armature{m_arm*1e3:7.1f} g  "
          f"({100*m_arm/M_PROJ:.0f}% of launch mass -- NOT optional)")
    print(f"    rhenium payload            {(M_PROJ-m_arm)*1e3:7.1f} g")
    print(f"    rhenium column length      {l_col*1e3:.1f} mm")
    print(f"    armature L_a={L_a*1e9:.1f} nH  R_a={R_a*1e6:.1f} uOhm  "
          f"L/R={L_a/R_a*1e3:.2f} ms")
    print(f"    base stress at {A_MAX_G:,.0f} g   {s50/1e6:.0f} MPa "
          f"(Al-6061 yield 276, SF {276e6/s50:.1f})")

    # ---- headline designs ------------------------------------------------
    print(f"\n[1] DESIGN POINTS (a_max <= {A_MAX_G:,.0f} g, B <= {B_MAX_T:.0f} T)\n")
    rows = []
    designs = {}
    for vt in (300.0, 400.0, 500.0, 619.0, 700.0):
        r = multistage(vt, 30, 0.04, 5.26e-6)
        if not r:
            continue
        designs[vt] = r
        rows.append(dict(v_target=vt, v_actual=r["v_final"],
                         n_stages=r["n_stages"], length_m=r["length"],
                         E_stored_kJ=r["E_stored_total"] / 1e3,
                         eta_pct=100 * r["eta_total"],
                         KE_J=0.5 * M_PROJ * r["v_final"] ** 2,
                         m_cap_kg=r["m_cap"], m_cu_kg=r["m_cu"],
                         a_max_g=r["a_max_g"], B_peak_T=r["B_peak"],
                         recoil_Ns=M_PROJ * r["v_final"],
                         charge_s_at_5kW=r["E_stored_total"] / 5000.0
                         / ETA_CHARGER))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "coil01_design_points.csv", index=False)
    print(df.to_string(index=False, float_format=lambda x: f"{x:.2f}"))

    # ---- acceleration-limit trade ----------------------------------------
    print("\n[2] ACCELERATION LIMIT vs LENGTH AND ENERGY (target 400 m/s)\n")
    rows = []
    for ag in (20_000.0, 50_000.0, 100_000.0, 200_000.0):
        r = multistage(400.0, 30, 0.04, 5.26e-6, a_max_g=ag)
        if not r:
            continue
        st, _ = canister_base_stress(ag, M_PROJ - m_arm)
        rows.append(dict(a_limit_g=ag, n_stages=r["n_stages"],
                         length_m=r["length"],
                         E_stored_kJ=r["E_stored_total"] / 1e3,
                         eta_pct=100 * r["eta_total"],
                         m_cap_kg=r["m_cap"], B_peak_T=r["B_peak"],
                         base_stress_MPa=st / 1e6,
                         canister_SF=276e6 / st))
    da = pd.DataFrame(rows)
    da.to_csv(OUT / "coil01_accel_trade.csv", index=False)
    print(da.to_string(index=False, float_format=lambda x: f"{x:.2f}"))

    # ---- coil turns ------------------------------------------------------
    print("\n[3] COIL TURNS (target 400 m/s)\n")
    rows = []
    for N in (15, 22, 30, 45, 60):
        r = multistage(400.0, N, 0.04, 5.26e-6)
        if not r:
            continue
        Lc, Rc, mcu = coil_params(N, 0.04, 5.26e-6)
        rows.append(dict(N_turns=N, L_coil_uH=Lc * 1e6, R_coil_mOhm=Rc * 1e3,
                         n_stages=r["n_stages"], length_m=r["length"],
                         E_stored_kJ=r["E_stored_total"] / 1e3,
                         eta_pct=100 * r["eta_total"],
                         V_typ_kV=np.median([s["V_kV"] for s in r["stages"]]),
                         m_cu_kg=r["m_cu"], m_cap_kg=r["m_cap"]))
    dn = pd.DataFrame(rows)
    dn.to_csv(OUT / "coil01_turns.csv", index=False)
    print(dn.to_string(index=False, float_format=lambda x: f"{x:.2f}"))

    # ---- stage-by-stage profile of the recommended design ----------------
    best = designs.get(400.0)
    ds = pd.DataFrame(best["stages"])
    ds.to_csv(OUT / "coil01_stage_profile.csv", index=False)
    print(f"\n[4] STAGE PROFILE, 400 m/s design ({best['n_stages']} stages)\n")
    print(ds.iloc[::4].to_string(index=False, float_format=lambda x: f"{x:.2f}"))

    # ---- coilgun vs chemical gun -----------------------------------------
    print("\n[5] COILGUN vs CHEMICAL GUN (both 100 g to 400 m/s)\n")
    # chemical figures come from GUN-01 (25 mm bore, 0.6 m, 12.0 g Retumbo)
    chem = dict(length_m=0.60, prime_mass_kg=0.4 + 5.0,   # barrel + breech
                consumable_g_per_shot=12.0, recoil_Ns=43.4,
                heat_per_shot_kJ=10.7, cadence_s=150.0, erosion="yes")
    coil = dict(length_m=best["length"],
                prime_mass_kg=best["m_cap"] + best["m_cu"] + 0.30 *
                (best["m_cap"] + best["m_cu"]),
                consumable_g_per_shot=0.0,
                recoil_Ns=M_PROJ * best["v_final"],
                heat_per_shot_kJ=best["q_total"] / 1e3,
                cadence_s=best["E_stored_total"] / 5000.0 / ETA_CHARGER,
                erosion="no")
    cmp_rows = []
    for k in ("length_m", "prime_mass_kg", "consumable_g_per_shot",
              "recoil_Ns", "heat_per_shot_kJ", "cadence_s", "erosion"):
        cmp_rows.append(dict(metric=k, coilgun=coil[k], chemical=chem[k]))
    dc = pd.DataFrame(cmp_rows)
    dc.to_csv(OUT / "coil01_vs_chemical.csv", index=False)
    print(dc.to_string(index=False))

    # mass crossover
    n_cross = (coil["prime_mass_kg"] - chem["prime_mass_kg"]) / \
              (chem["consumable_g_per_shot"] * 2.0e-3)   # x2 for cases
    print(f"\n    mass crossover: the coilgun becomes lighter after "
          f"{n_cross:,.0f} shots")
    print(f"    (chemical consumable charged at 2x propellant mass to include "
          f"cases/primers)")
    print(f"    FG01 realised mission is ~24 engagements per platform-decade "
          f"(ECO-1),")
    print(f"    so the crossover is NOT reached by a wide margin.")

    summary = dict(
        configuration=dict(bore_mm=BORE * 1e3, m_total_g=M_PROJ * 1e3,
                           m_armature_g=m_arm * 1e3,
                           m_rhenium_g=(M_PROJ - m_arm) * 1e3),
        design_400=dict(n_stages=best["n_stages"], length_m=best["length"],
                        E_stored_kJ=best["E_stored_total"] / 1e3,
                        eta_pct=100 * best["eta_total"],
                        m_cap_kg=best["m_cap"], m_cu_kg=best["m_cu"],
                        a_max_g=best["a_max_g"], B_peak_T=best["B_peak"]),
        mass_crossover_shots=float(n_cross))
    if 619.0 in designs:
        d6 = designs[619.0]
        summary["design_619"] = dict(
            n_stages=d6["n_stages"], length_m=d6["length"],
            E_stored_kJ=d6["E_stored_total"] / 1e3,
            eta_pct=100 * d6["eta_total"], m_cap_kg=d6["m_cap"])
    with open(OUT / "coil01_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2, default=float)

    # ---- figures ---------------------------------------------------------
    fig, ax = plt.subplots(1, 3, figsize=(13.5, 3.9))
    ax[0].plot(df.v_actual, df.E_stored_kJ, color=CB[0], marker="o", ms=4)
    ax[0].set_xlabel("muzzle velocity (m/s)")
    ax[0].set_ylabel("capacitor energy stored (kJ)")
    ax[0].set_title("Stored energy vs velocity")
    ax[1].plot(df.v_actual, df.length_m, color=CB[1], marker="s", ms=4)
    ax[1].set_xlabel("muzzle velocity (m/s)")
    ax[1].set_ylabel("coilgun length (m)")
    ax[1].set_title("Length vs velocity")
    ax[2].plot(da.a_limit_g / 1e3, da.length_m, color=CB[2], marker="o", ms=4,
               label="length (m)")
    a2 = ax[2].twinx()
    a2.plot(da.a_limit_g / 1e3, da.E_stored_kJ, color=CB[3], marker="s", ms=4,
            label="energy (kJ)")
    ax[2].set_xlabel("acceleration limit (thousand g)")
    ax[2].set_ylabel("length (m)", color=CB[2])
    a2.set_ylabel("stored energy (kJ)", color=CB[3])
    ax[2].set_title("Acceleration limit trade")
    fig.suptitle("COIL-01: 100 g in a 3 cm bore", fontsize=10)
    for e in ("svg", "png"):
        fig.savefig(FIG / f"coil01_design.{e}")
    plt.close(fig)

    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    ax[0].plot(ds.stage, ds.v_out, color=CB[0], marker="o", ms=3)
    ax[0].set_xlabel("stage number"); ax[0].set_ylabel("velocity (m/s)")
    ax[0].set_title("Velocity build-up")
    ax[1].plot(ds.stage, ds.E_J, color=CB[1], marker="o", ms=3, label="energy/stage")
    a2 = ax[1].twinx()
    a2.plot(ds.stage, ds.eta_pct, color=CB[2], marker="s", ms=3, label="efficiency")
    ax[1].set_xlabel("stage number"); ax[1].set_ylabel("stage energy (J)", color=CB[1])
    a2.set_ylabel("stage efficiency (%)", color=CB[2])
    ax[1].set_title("Stage grading")
    fig.suptitle(f"COIL-01 stage profile, {best['n_stages']} stages to 400 m/s",
                 fontsize=10)
    for e in ("svg", "png"):
        fig.savefig(FIG / f"coil01_stages.{e}")
    plt.close(fig)

    print("\nWrote outputs/coil01_*.csv and figures/coil01_*.svg")
    return df, da, best


run = main          # run_all.py calls .run()


if __name__ == "__main__":
    main()
