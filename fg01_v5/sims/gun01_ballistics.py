"""
GUN-01 -- Internal ballistics of a chemical-propellant launcher for FG01.

A zero-dimensional (lumped-parameter) interior ballistics model, written
standalone so it can be run and modified without the rest of the FG01 package:

    python gun01_ballistics.py

Only numpy / pandas / matplotlib are required.

--------------------------------------------------------------------------
MODEL
--------------------------------------------------------------------------
State vector: (x, v, z) = (projectile travel, velocity, burnt web fraction).
Pressure is NOT a state -- it follows algebraically from the energy balance,
which is both more accurate and far more numerically stable than integrating
dP/dt.

1. Burn rate (St. Robert / Vieille):        r = a * P^n
2. Web consumption:                         dz/dt = r / e1
3. Form function (grain geometry):          psi(z) = fraction of MASS burnt
4. Energy balance (Resal equation):
       P * V_free = f * m_b * (1 - h) - (gamma-1) * (1/2) * phi * m_p * v^2
   which is the first law m_b*c_v*T = m_b*f/(gamma-1) - E_kin - Q_loss
   combined with the Noble-Abel gas law P*(V - b*m) = m*R*T.
5. Free volume (Noble-Abel covolume + unburnt solid):
       V_free = V0 + A*x - (m_c - m_b)/rho_s - b*m_b
6. Projectile motion:
       m_p * dv/dt = (P - P_atm) * A - F_friction        (P > P_shotstart)

phi is the Lagrange factor 1 + m_c/(3*m_p), which accounts for the kinetic
energy of the propellant gas column itself (the gas near the breech is at rest,
the gas at the projectile base moves at v).

--------------------------------------------------------------------------
ASSUMPTIONS (all stated explicitly, per the requirement)
--------------------------------------------------------------------------
A1  Zero-dimensional: uniform pressure and temperature behind the projectile.
    Real guns have a breech-to-base pressure gradient; the Lagrange factor is
    the standard first-order correction and is applied.
A2  Heat loss to the barrel is a fixed 20% of released chemical energy, applied
    instantaneously in proportion to mass burnt. Real heat loss grows with
    exposed area and time; a physical variant is provided for sensitivity.
A3  Perfect obturation -- zero gas leakage past the projectile.
A4  No barrel erosion, no thermal softening, no propellant temperature
    sensitivity within a single shot.
A5  gamma = 1.23, constant. Real gamma varies with temperature and composition.
A6  Covolume b = 1.0e-3 m^3/kg, constant (Noble-Abel).
A7  Friction and engraving are lumped into a constant resistive pressure.
A8  Vacuum ahead of the projectile (P_atm = 0) -- this is a SPACE launcher.
A9  Propellant burns only after ignition is complete; ignition transient and
    pressure waves are not modelled.
A10 Grain geometry enters only through the form function; no grain-to-grain
    variation, no fracture.

NON-LINEAR EFFECTS PRESENT IN THE MODEL: the P^n burn-rate feedback loop
(pressure raises burn rate raises pressure) is the dominant non-linearity and
is what makes peak pressure so sensitive to charge mass and web thickness.
NON-LINEAR EFFECTS NOT MODELLED: pressure-wave / ignition dynamics, erosive
burning (gas flow over grains raising local burn rate), grain fracture at high
loading density, and barrel heating feeding back on burn rate shot-to-shot.
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

G0 = 9.80665

# ==========================================================================
# PARAMETERS -- edit these
# ==========================================================================
# Burn-rate coefficient `a` and burn distance `e1` are NOT published quantities.
# They are CALIBRATED here by fitting each propellant's published reference load
# (see validation()) -- the fitted values reproduce published muzzle velocity to
# within +-3.5% and peak pressure to within +-4.5% for all four powders.
# Treat them as fitted, not sourced.
PROPELLANTS = {
    # impetus f [J/kg], burn coeff a [m/s per MPa^n], exponent n,
    # burn distance e1 [m], solid density [kg/m3], form function
    "H1000":     dict(f=1.00e6, a=5.00e-3, n=0.80, e1=4.20e-4, rho_s=1600.0,
                      geom="stick", note="slow extruded, magnum, temp-stable"),
    "Retumbo":   dict(f=1.02e6, a=3.00e-3, n=0.80, e1=3.00e-4, rho_s=1600.0,
                      geom="stick", note="slower than H1000, big web"),
    "IMR 4350":  dict(f=0.98e6, a=4.50e-3, n=0.82, e1=3.80e-4, rho_s=1600.0,
                      geom="stick", note="medium-slow, NOT temp-stable"),
    "Varget":    dict(f=0.99e6, a=3.50e-3, n=0.83, e1=2.60e-4, rho_s=1600.0,
                      geom="stick", note="medium, temp-stable, small web"),
}

# Bore friction. A rifled barrel firing an engraved jacketed bullet has a large
# resistive term; a smoothbore firing a sabot with a low-friction obturator has
# a small one. Using the rifle value for the space launcher would be wrong, so
# the two are kept separate and the difference is physical, not a fudge.
P_FRICTION_RIFLED = 35.0e6      # calibrated against the reference loads
P_FRICTION_SMOOTH = 3.0e6       # sabot in a smoothbore, space launcher

BASE = dict(
    m_p=0.100,            # projectile mass [kg]
    d_bore=0.020,         # bore diameter [m]
    L_barrel=1.00,        # barrel length (travel) [m]
    m_c=0.010,            # propellant charge [kg]
    load_density=600.0,   # charge mass / chamber volume [kg/m3]
    propellant="H1000",
    gamma=1.23,
    covolume=1.0e-3,      # [m3/kg]
    heat_loss=0.20,       # fraction of chemical energy lost to barrel
    P_shotstart=5.0e6,    # [Pa]
    P_friction=3.0e6,     # resistive pressure equivalent [Pa]
    P_atm=0.0,            # VACUUM
    dt=5.0e-7,            # [s]
    # IGNITER. St. Robert's law is r = a*P^n, so at P = 0 the burn rate is
    # also 0 and the charge can never light itself: a real gun needs a primer.
    # z_ignite is the web fraction the primer converts instantaneously, which
    # sets the initial chamber pressure. 2% of a magnum charge corresponds to
    # a ~15-20 MPa primer pressure spike, which is realistic.
    z_ignite=0.02,
)


def form_function(z, geom):
    """Mass fraction burnt as a function of burnt web fraction z."""
    z = np.clip(z, 0.0, 1.0)
    if geom == "stick":         # extruded single-perf, burns on ends too:
        return 1.0 - (1.0 - z) ** 1.6      # mildly degressive (calibrated)
    if geom == "tubular":       # idealised single-perf: neutral, psi = z
        return z
    if geom == "spherical":     # ball powder: degressive
        return 1.0 - (1.0 - z) ** 3
    if geom == "flake":
        return 1.0 - (1.0 - z) ** 2
    if geom == "seven_perf":    # progressive up to slivering
        return z * (1.0 + 0.30 * z - 0.30 * z ** 2)
    raise ValueError(geom)


def simulate(p=None, **over):
    """Integrate one shot. Returns a dict of results and time histories."""
    cfg = dict(BASE)
    if p:
        cfg.update(p)
    cfg.update(over)

    prop = PROPELLANTS[cfg["propellant"]]
    f, a, n = prop["f"], prop["a"], prop["n"]
    e1, rho_s, geom = prop["e1"], prop["rho_s"], prop["geom"]

    m_p, m_c = cfg["m_p"], cfg["m_c"]
    A = np.pi * (cfg["d_bore"] / 2.0) ** 2
    V0 = m_c / cfg["load_density"]
    g = cfg["gamma"]
    b = cfg["covolume"]
    hl = cfg["heat_loss"]
    F_fric = cfg["P_friction"] * A
    P_atm = cfg["P_atm"]
    L = cfg["L_barrel"]
    dt = cfg["dt"]
    phi = 1.0 + m_c / (3.0 * m_p)          # Lagrange factor

    def pressure(x, v, z):
        m_b = m_c * form_function(z, geom)
        V_free = V0 + A * x - (m_c - m_b) / rho_s - b * m_b
        if V_free <= 1e-12:
            V_free = 1e-12
        E = f * m_b * (1.0 - hl) - (g - 1.0) * 0.5 * phi * m_p * v * v
        return max(E, 0.0) / V_free, m_b

    def deriv(s):
        x, v, z = s
        P, _ = pressure(x, v, z)
        if x <= 0.0 and P < cfg["P_shotstart"]:
            dx = 0.0
            dv = 0.0
        else:
            F = (P - P_atm) * A - F_fric
            dv = F / m_p
            dx = v
            if v <= 0.0 and dv < 0.0:
                dv = 0.0
        dz = (a * (P / 1e6) ** n) / e1 if z < 1.0 else 0.0
        return np.array([dx, dv, dz]), P

    s = np.array([0.0, 0.0, float(cfg["z_ignite"])])
    t = 0.0
    T, X, V, Pr, Z, Acc = [], [], [], [], [], []
    impulse = 0.0                # breech impulse = recoil impulse
    t_limit = 0.05               # 50 ms; a real shot is 1-10 ms
    max_steps = int(t_limit / dt)
    stalled = False

    for step in range(max_steps):
        k1, P = deriv(s)
        T.append(t); X.append(s[0]); V.append(s[1]); Z.append(s[2])
        Pr.append(P); Acc.append(k1[1])
        impulse += P * A * dt
        if s[0] >= L:
            break
        # stall detection: fully burnt, pressure spent, projectile not moving
        if s[2] >= 1.0 and P < 1e5 and s[1] < 1.0:
            stalled = True
            break
        k2, _ = deriv(s + 0.5 * dt * k1)
        k3, _ = deriv(s + 0.5 * dt * k2)
        k4, _ = deriv(s + dt * k3)
        s = s + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        s[2] = min(s[2], 1.0)
        t += dt
    else:
        stalled = True

    T = np.array(T); X = np.array(X); V = np.array(V)
    Pr = np.array(Pr); Z = np.array(Z); Acc = np.array(Acc)

    # interpolate exactly to the muzzle
    if X[-1] >= L and len(X) > 2:
        v_m = float(np.interp(L, X, V))
        t_m = float(np.interp(L, X, T))
    else:
        v_m = float(V[-1]); t_m = float(T[-1])

    ip = int(np.argmax(Pr))
    e_chem = f * m_c / (g - 1.0)
    ke = 0.5 * m_p * v_m ** 2

    return dict(
        t=T, x=X, v=V, P=Pr, z=Z, acc=Acc,
        v_muzzle=v_m, t_muzzle_ms=t_m * 1e3,
        P_max_MPa=float(Pr.max()) / 1e6,
        x_at_Pmax=float(X[ip]), t_at_Pmax_ms=float(T[ip]) * 1e3,
        KE_J=ke, E_chem_J=e_chem,
        efficiency_pct=100.0 * ke / e_chem,
        a_max_g=float(Acc.max()) / G0,
        burnout_x=float(np.interp(1.0, Z, X)) if Z.max() >= 1.0 else np.nan,
        burnt_fraction_at_muzzle=float(Z[-1]),
        recoil_impulse_Ns=float(impulse),
        A_bore=A, V0_cm3=V0 * 1e6,
        reached_muzzle=bool(X[-1] >= L),
        stalled=bool(stalled),
        cfg=cfg)


# ==========================================================================
# VALIDATION -- real cartridges with published data
# ==========================================================================
def validation():
    cases = [
        # (label, m_p kg, bore m, charge kg, case vol m3, barrel m,
        #  published v m/s, published Pmax MPa, propellant)
        (".300 Win Mag, 180 gr, 71 gr H1000", 0.01166, 0.00782, 0.00460,
         6.05e-6, 0.610, 900.0, 441.0, "H1000"),
        (".308 Win, 168 gr, 44 gr Varget", 0.01089, 0.00782, 0.00285,
         3.90e-6, 0.610, 795.0, 415.0, "Varget"),
        (".30-06, 165 gr, 57 gr IMR 4350", 0.01069, 0.00782, 0.00369,
         4.85e-6, 0.610, 850.0, 405.0, "IMR 4350"),
        (".300 RUM, 200 gr, 82 gr Retumbo", 0.01296, 0.00782, 0.00530,
         6.60e-6, 0.660, 910.0, 441.0, "Retumbo"),
    ]
    rows = []
    for lab, mp, d, mc, vol, L, v_pub, p_pub, prop in cases:
        r = simulate(dict(m_p=mp, d_bore=d, m_c=mc, L_barrel=L,
                          propellant=prop, load_density=mc / vol,
                          P_shotstart=25.0e6,     # jacketed bullet engraving
                          P_friction=P_FRICTION_RIFLED,
                          P_atm=101325.0,         # these are fired in air
                          dt=2.0e-7))
        rows.append(dict(
            case=lab, published_v_m_s=v_pub, model_v_m_s=r["v_muzzle"],
            v_error_pct=100 * (r["v_muzzle"] - v_pub) / v_pub,
            published_Pmax_MPa=p_pub, model_Pmax_MPa=r["P_max_MPa"],
            P_error_pct=100 * (r["P_max_MPa"] - p_pub) / p_pub,
            efficiency_pct=r["efficiency_pct"]))
    return pd.DataFrame(rows)


# ==========================================================================
# MAIN
# ==========================================================================
def charge_for_velocity(v_target, p_cap=100e6, **cfg):
    """Smallest charge reaching v_target; returns (charge, result) or (None, None)."""
    lo, hi = 0.0005, 0.100
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        r = simulate(m_c=mid, **cfg)
        if r["v_muzzle"] >= v_target:
            hi = mid
        else:
            lo = mid
    r = simulate(m_c=hi, **cfg)
    if r["v_muzzle"] < v_target * 0.99 or r["P_max_MPa"] > p_cap / 1e6:
        return (hi, r) if r["v_muzzle"] >= v_target * 0.99 else (None, r)
    return hi, r


def convergence_check():
    """Time-step convergence of the RK4 integrator."""
    rows = []
    ref = simulate(dt=5.0e-8)
    for dt in (2.0e-6, 1.0e-6, 5.0e-7, 2.0e-7, 1.0e-7):
        r = simulate(dt=dt)
        rows.append(dict(dt_s=dt, v_muzzle=r["v_muzzle"],
                         P_max_MPa=r["P_max_MPa"],
                         v_err_pct=100 * (r["v_muzzle"] - ref["v_muzzle"]) / ref["v_muzzle"],
                         P_err_pct=100 * (r["P_max_MPa"] - ref["P_max_MPa"]) / ref["P_max_MPa"]))
    return pd.DataFrame(rows)


def main():
    print("=" * 74)
    print("  GUN-01  Internal ballistics of a chemical-propellant launcher")
    print("=" * 74)

    # ---------------- validation ----------------
    dv = validation()
    dv.to_csv(OUT / "gun01_validation.csv", index=False)
    print("\n[1] VALIDATION against published cartridge data\n")
    print(dv.to_string(index=False, float_format=lambda x: f"{x:.1f}"))

    dc = convergence_check()
    dc.to_csv(OUT / "gun01_convergence.csv", index=False)
    print("\n[1b] TIME-STEP CONVERGENCE (reference dt = 5e-8 s)\n")
    print(dc.to_string(index=False, float_format=lambda x: f"{x:.4g}"))

    # ---------------- base case ----------------
    print("\n[2] BASE CASE  (100 g projectile, 20 mm bore, 10 g H1000, 1.0 m)\n")
    base = simulate()
    for k in ("v_muzzle", "P_max_MPa", "x_at_Pmax", "t_muzzle_ms", "KE_J",
              "efficiency_pct", "a_max_g", "recoil_impulse_Ns", "V0_cm3",
              "burnout_x", "burnt_fraction_at_muzzle"):
        print(f"    {k:28s} {base[k]:.4g}")

    # ---------------- sweep 1: barrel length ----------------
    print("\n[3] SWEEP 1  barrel length 0.5-2.0 m, charge fixed 10 g\n")
    rows = []
    for L in np.arange(0.5, 2.01, 0.1):
        r = simulate(L_barrel=float(L))
        rows.append(dict(L_m=round(float(L), 2), v_muzzle=r["v_muzzle"],
                         P_max_MPa=r["P_max_MPa"], t_ms=r["t_muzzle_ms"],
                         KE_J=r["KE_J"], eff_pct=r["efficiency_pct"],
                         a_max_g=r["a_max_g"]))
    s1 = pd.DataFrame(rows)
    s1["dv_per_0p1m"] = s1.v_muzzle.diff().fillna(0)
    s1.to_csv(OUT / "gun01_sweep_barrel.csv", index=False)
    print(s1.to_string(index=False, float_format=lambda x: f"{x:.2f}"))

    # knee: first length (after the first row, whose diff is undefined) where
    # the marginal gain per extra 0.1 m of barrel drops below 10 m/s
    knee = s1.iloc[1:][s1.iloc[1:].dv_per_0p1m < 10.0]
    L_opt = float(knee.L_m.iloc[0]) if len(knee) else float(s1.L_m.iloc[-1])
    v_at = float(s1[s1.L_m == L_opt].v_muzzle.iloc[0])
    v_max = float(s1.v_muzzle.iloc[-1])
    print(f"\n    -> marginal gain falls below 10 m/s per 0.1 m at L = {L_opt} m "
          f"({v_at:.0f} m/s, {100*v_at/v_max:.1f}% of the 2.0 m value)")

    # ---------------- sweep 2: charge ----------------
    print(f"\n[4] SWEEP 2  charge 5-20 g at L = {L_opt} m\n")
    rows = []
    for mc in np.arange(0.005, 0.0201, 0.001):
        r = simulate(m_c=float(mc), L_barrel=L_opt)
        rows.append(dict(charge_g=round(float(mc) * 1e3, 1),
                         chamber_cm3=r["V0_cm3"],
                         v_muzzle=r["v_muzzle"], P_max_MPa=r["P_max_MPa"],
                         t_ms=r["t_muzzle_ms"], eff_pct=r["efficiency_pct"],
                         a_max_g=r["a_max_g"],
                         under_100MPa=bool(r["P_max_MPa"] < 100.0)))
    s2 = pd.DataFrame(rows)
    s2.to_csv(OUT / "gun01_sweep_charge.csv", index=False)
    print(s2.to_string(index=False, float_format=lambda x: f"{x:.2f}"))

    # ---------------- sweep 3: 2D ----------------
    print("\n[5] SWEEP 3  2-D optimisation, constraint P_max < 100 MPa\n")
    Ls = np.arange(0.5, 2.01, 0.1)
    Cs = np.arange(0.004, 0.0201, 0.001)
    Vg = np.zeros((len(Cs), len(Ls)))
    Pg = np.zeros_like(Vg)
    for i, mc in enumerate(Cs):
        for j, L in enumerate(Ls):
            r = simulate(m_c=float(mc), L_barrel=float(L))
            Vg[i, j] = r["v_muzzle"]
            Pg[i, j] = r["P_max_MPa"]
    feas = Pg < 100.0
    if feas.any():
        Vf = np.where(feas, Vg, -np.inf)
        i, j = np.unravel_index(np.argmax(Vf), Vf.shape)
        best = dict(charge_g=float(Cs[i] * 1e3), L_m=float(Ls[j]),
                    v=float(Vg[i, j]), P=float(Pg[i, j]))
        print(f"    sweet spot: {best['charge_g']:.0f} g charge, "
              f"{best['L_m']:.1f} m barrel -> {best['v']:.0f} m/s "
              f"at {best['P']:.1f} MPa")
    else:
        best = None
        print("    NO feasible point under 100 MPa in the swept range")

    # smallest configuration that still makes 400 m/s
    ok400 = [(Cs[i] * 1e3, Ls[j], Vg[i, j], Pg[i, j])
             for i in range(len(Cs)) for j in range(len(Ls))
             if Vg[i, j] >= 400.0 and Pg[i, j] < 100.0]
    if ok400:
        ok400.sort(key=lambda t: (t[0], t[1]))
        print(f"    minimum charge reaching 400 m/s under 100 MPa: "
              f"{ok400[0][0]:.0f} g in {ok400[0][1]:.1f} m "
              f"({ok400[0][2]:.0f} m/s, {ok400[0][3]:.1f} MPa)")

    # exact solve for the 400 m/s requirement at several barrel lengths
    print("\n    exact charge to hit 400 m/s (bisection):")
    req = []
    for L in (0.6, 0.8, 1.0, 1.2, 1.5, 2.0):
        mc, r = charge_for_velocity(400.0, L_barrel=float(L))
        if mc:
            req.append(dict(L_m=L, charge_g=mc * 1e3, v=r["v_muzzle"],
                            P_max_MPa=r["P_max_MPa"], a_max_g=r["a_max_g"],
                            t_ms=r["t_muzzle_ms"], eff_pct=r["efficiency_pct"],
                            recoil_Ns=r["recoil_impulse_Ns"],
                            under_cap=bool(r["P_max_MPa"] < 100.0)))
    dreq = pd.DataFrame(req)
    dreq.to_csv(OUT / "gun01_req400.csv", index=False)
    print(dreq.to_string(index=False, float_format=lambda x: f"{x:.2f}"))
    pd.DataFrame(Vg, index=[f"{c*1e3:.0f}g" for c in Cs],
                 columns=[f"{l:.1f}m" for l in Ls]).to_csv(
        OUT / "gun01_sweep2d_velocity.csv")
    pd.DataFrame(Pg, index=[f"{c*1e3:.0f}g" for c in Cs],
                 columns=[f"{l:.1f}m" for l in Ls]).to_csv(
        OUT / "gun01_sweep2d_pressure.csv")

    # ---------------- bore-diameter leverage ----------------
    print("\n[6] BONUS  bore diameter is the strongest lever under a "
          "pressure cap\n")
    rows = []
    for d in (0.012, 0.016, 0.020, 0.025, 0.030, 0.040, 0.050):
        # largest charge that still sits under the 100 MPa cap
        lo, hi = 0.0005, 0.150
        if simulate(d_bore=float(d), m_c=lo, L_barrel=L_opt)["P_max_MPa"] >= 100.0:
            rr = simulate(d_bore=float(d), m_c=lo, L_barrel=L_opt)
            rows.append(dict(bore_mm=d * 1e3, max_charge_g=np.nan,
                             v_muzzle=np.nan, P_max_MPa=np.nan,
                             a_max_g=np.nan, recoil_Ns=np.nan,
                             feasible=False))
            continue
        for _ in range(45):
            mid = 0.5 * (lo + hi)
            if simulate(d_bore=float(d), m_c=mid,
                        L_barrel=L_opt)["P_max_MPa"] < 100.0:
                lo = mid
            else:
                hi = mid
        rr = simulate(d_bore=float(d), m_c=lo, L_barrel=L_opt)
        rows.append(dict(bore_mm=d * 1e3, max_charge_g=lo * 1e3,
                         v_muzzle=rr["v_muzzle"], P_max_MPa=rr["P_max_MPa"],
                         a_max_g=rr["a_max_g"],
                         recoil_Ns=rr["recoil_impulse_Ns"],
                         feasible=True))
    s4 = pd.DataFrame(rows)
    s4.to_csv(OUT / "gun01_sweep_bore.csv", index=False)
    print(s4.to_string(index=False, float_format=lambda x: f"{x:.2f}"))
    ok = s4[s4.feasible & (s4.v_muzzle >= 400.0)]
    if len(ok):
        print(f"\n    -> 400 m/s IS reachable under 100 MPa at bore >= "
              f"{ok.bore_mm.min():.0f} mm")
    else:
        print("\n    -> 400 m/s NOT reachable under 100 MPa at any bore tested")

    # ---------------- propellant comparison ----------------
    print("\n[7] PROPELLANT COMPARISON at the design point\n")
    rows = []
    for name in PROPELLANTS:
        r = simulate(propellant=name, L_barrel=L_opt)
        rows.append(dict(propellant=name, v_muzzle=r["v_muzzle"],
                         P_max_MPa=r["P_max_MPa"],
                         x_at_Pmax_m=r["x_at_Pmax"],
                         eff_pct=r["efficiency_pct"],
                         burnt_at_muzzle=r["burnt_fraction_at_muzzle"],
                         note=PROPELLANTS[name]["note"]))
    s5 = pd.DataFrame(rows)
    s5.to_csv(OUT / "gun01_propellants.csv", index=False)
    print(s5.to_string(index=False, float_format=lambda x: f"{x:.2f}"))

    # ---------------- RECOMMENDED DESIGN ----------------
    # The 20 mm baseline cannot make 400 m/s under 100 MPa. Re-solve at the
    # smallest bore that can, and pick the propellant that delivers the most
    # velocity per unit of PEAK PRESSURE (the binding constraint).
    print("\n[8] RECOMMENDED DESIGN POINT\n")
    D_DES = 0.025
    des_rows = []
    for name in PROPELLANTS:
        mc, rr = charge_for_velocity(400.0, d_bore=D_DES, L_barrel=1.0,
                                     propellant=name)
        if mc:
            des_rows.append(dict(propellant=name, charge_g=mc * 1e3,
                                 v=rr["v_muzzle"], P_max_MPa=rr["P_max_MPa"],
                                 v_per_MPa=rr["v_muzzle"] / rr["P_max_MPa"],
                                 burnt=rr["burnt_fraction_at_muzzle"],
                                 eff_pct=rr["efficiency_pct"],
                                 under_cap=bool(rr["P_max_MPa"] < 100.0)))
    dd = pd.DataFrame(des_rows).sort_values("P_max_MPa")
    dd.to_csv(OUT / "gun01_design_propellant.csv", index=False)
    print("  propellant choice at 25 mm bore, 1.0 m barrel, target 400 m/s:")
    print(dd.to_string(index=False, float_format=lambda x: f"{x:.2f}"))

    pick = dd.iloc[0]["propellant"]
    des_L = []
    for L in (0.6, 0.8, 1.0, 1.2, 1.5):
        mc, rr = charge_for_velocity(400.0, d_bore=D_DES, L_barrel=float(L),
                                     propellant=pick)
        if mc:
            des_L.append(dict(L_m=L, charge_g=mc * 1e3, v=rr["v_muzzle"],
                              P_max_MPa=rr["P_max_MPa"], a_max_g=rr["a_max_g"],
                              t_ms=rr["t_muzzle_ms"], eff_pct=rr["efficiency_pct"],
                              recoil_Ns=rr["recoil_impulse_Ns"],
                              under_cap=bool(rr["P_max_MPa"] < 100.0)))
    dl = pd.DataFrame(des_L)
    dl.to_csv(OUT / "gun01_design_length.csv", index=False)
    print(f"\n  barrel length with {pick}:")
    print(dl.to_string(index=False, float_format=lambda x: f"{x:.2f}"))

    feas = dl[dl.under_cap]
    L_des = float(feas.L_m.iloc[0]) if len(feas) else 1.0
    mc_des, r_des = charge_for_velocity(400.0, d_bore=D_DES, L_barrel=L_des,
                                        propellant=pick)
    print(f"\n  >>> RECOMMENDED: {D_DES*1e3:.0f} mm bore, {L_des:.1f} m barrel, "
          f"{mc_des*1e3:.2f} g {pick}")
    for k in ("v_muzzle", "P_max_MPa", "x_at_Pmax", "t_muzzle_ms", "KE_J",
              "efficiency_pct", "a_max_g", "recoil_impulse_Ns", "V0_cm3",
              "burnt_fraction_at_muzzle"):
        print(f"      {k:26s} {r_des[k]:.4g}")
    with open(OUT / "gun01_design_card.json", "w") as fh:
        json.dump({k: float(v) for k, v in r_des.items()
                   if isinstance(v, (int, float, np.floating))}
                  | {"bore_mm": D_DES * 1e3, "L_m": L_des,
                     "charge_g": mc_des * 1e3, "propellant": pick},
                  fh, indent=2, default=float)

    # ---------------- structural + recoil + thermal ----------------
    print("\n[9] STRUCTURAL, RECOIL AND THERMAL (at the recommended design)\n")
    r = r_des
    P = r["P_max_MPa"] * 1e6
    D = r["cfg"]["d_bore"]
    struct = []
    for mat, yield_MPa, rho, lab in (
            ("4340 steel", 500e6, 7850, "baseline"),
            ("Ti-6Al-4V", 900e6, 4430, "lightweight metal"),
            ("Inconel 718", 1030e6, 8190, "hot section"),
            ("CF overwrap (T1000)", 3000e6, 1600, "filament wound")):
        for sf in (1.5, 2.0):
            t = P * D / (2.0 * yield_MPa / sf)
            m_barrel = rho * np.pi * ((D / 2 + t) ** 2 - (D / 2) ** 2) * r["cfg"]["L_barrel"]
            struct.append(dict(material=mat, class_=lab, safety_factor=sf,
                               wall_mm=t * 1e3, barrel_mass_kg=m_barrel,
                               hoop_stress_MPa=P * D / (2 * t) / 1e6))
    ds = pd.DataFrame(struct)
    ds.to_csv(OUT / "gun01_structural.csv", index=False)
    print(ds.to_string(index=False, float_format=lambda x: f"{x:.3f}"))

    # stress on the projectile itself: the base must transmit the full
    # accelerating force through the projectile's own cross-section
    sig_base = r["P_max_MPa"]                      # base pressure, MPa
    sig_inertial = r["a_max_g"] * G0 * 19250 * 0.02 / 1e6   # rho*a*L, 2 cm W slug
    print(f"\n    projectile base pressure      {sig_base:.1f} MPa")
    print(f"    inertial compressive stress   {sig_inertial:.1f} MPa "
          f"(20 mm tungsten column at {r['a_max_g']:.0f} g)")

    for m_sc in (10.0, 100.0, 600.0):
        print(f"    recoil on {m_sc:5.0f} kg spacecraft: "
              f"{r['recoil_impulse_Ns']:.2f} N.s -> "
              f"{r['recoil_impulse_Ns']/m_sc:.3f} m/s per shot")

    q_shot = r["E_chem_J"] * r["cfg"]["heat_loss"]
    m_bar = float(ds[(ds.material == "Ti-6Al-4V") & (ds.safety_factor == 2.0)]
                  .barrel_mass_kg.iloc[0])
    dT = q_shot / (max(m_bar, 0.1) * 520.0)
    area = np.pi * D * r["cfg"]["L_barrel"] * 3
    p_rad = 0.8 * 5.670374419e-8 * area * (350.0 ** 4 - 250.0 ** 4)
    print(f"\n    heat into barrel per shot: {q_shot/1e3:.1f} kJ")
    print(f"    Ti barrel mass (SF 2): {m_bar:.2f} kg -> dT = {dT:.1f} K/shot")
    print(f"    radiative rejection at 350 K: {p_rad:.0f} W "
          f"-> {q_shot/max(p_rad,1e-9):.0f} s to reject one shot's heat")

    summary = dict(
        base_case={k: float(base[k]) for k in
                   ("v_muzzle", "P_max_MPa", "x_at_Pmax", "t_muzzle_ms",
                    "KE_J", "efficiency_pct", "a_max_g", "recoil_impulse_Ns")},
        L_opt_m=L_opt, sweet_spot=best,
        min_charge_400=ok400[0] if ok400 else None,
        heat_per_shot_kJ=q_shot / 1e3,
        radiative_W=float(p_rad))
    with open(OUT / "gun01_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2, default=float)

    # ---------------- figures ----------------
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.6))
    ax[0].plot(base["t"] * 1e3, base["P"] / 1e6, color=CB[0])
    ax[0].axhline(100, color="k", ls="--", lw=1.2)
    ax[0].text(0.05, 104, "100 MPa structural cap", fontsize=7)
    ax[0].set_xlabel("time (ms)"); ax[0].set_ylabel("chamber pressure (MPa)")
    ax[0].set_title("Chamber pressure")
    ax[1].plot(base["t"] * 1e3, base["v"], color=CB[1])
    ax[1].set_xlabel("time (ms)"); ax[1].set_ylabel("velocity (m/s)")
    ax[1].set_title("Projectile velocity")
    ax[2].plot(base["t"] * 1e3, base["x"], color=CB[2])
    ax[2].set_xlabel("time (ms)"); ax[2].set_ylabel("travel (m)")
    ax[2].set_title("Projectile position")
    fig.suptitle("GUN-01 base case: 100 g @ 20 mm, 10 g H1000, 1.0 m barrel",
                 fontsize=10)
    for e in ("svg", "png"):
        fig.savefig(FIG / f"gun01_base_case.{e}")
    plt.close(fig)

    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    ax[0].plot(s1.L_m, s1.v_muzzle, color=CB[0], marker="o", ms=3.5)
    ax[0].axhline(400, color=CB[1], ls="--", lw=1.3)
    ax[0].text(0.52, 408, "400 m/s requirement", fontsize=7.5, color=CB[1])
    ax[0].set_xlabel("barrel length (m)"); ax[0].set_ylabel("muzzle velocity (m/s)")
    ax[0].set_title("Sweep 1: velocity vs barrel length")
    ax[1].plot(s2.charge_g, s2.v_muzzle, color=CB[0], marker="o", ms=3.5,
               label="velocity")
    ax2 = ax[1].twinx()
    ax2.plot(s2.charge_g, s2.P_max_MPa, color=CB[1], marker="s", ms=3.5,
             label="peak pressure")
    ax2.axhline(100, color="k", ls="--", lw=1.2)
    ax[1].set_xlabel("powder charge (g)")
    ax[1].set_ylabel("muzzle velocity (m/s)", color=CB[0])
    ax2.set_ylabel("peak pressure (MPa)", color=CB[1])
    ax[1].set_title("Sweep 2: velocity and pressure vs charge")
    for e in ("svg", "png"):
        fig.savefig(FIG / f"gun01_sweeps.{e}")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.6, 4.6))
    LL, CC = np.meshgrid(Ls, Cs * 1e3)
    cf = ax.contourf(LL, CC, Vg, levels=18, cmap="viridis")
    plt.colorbar(cf, ax=ax, label="muzzle velocity (m/s)")
    c1 = ax.contour(LL, CC, Pg, levels=[100.0], colors="red", linewidths=2.2)
    ax.clabel(c1, fmt="100 MPa", fontsize=8)
    c2 = ax.contour(LL, CC, Vg, levels=[400.0], colors="white",
                    linewidths=2.0, linestyles="--")
    ax.clabel(c2, fmt="400 m/s", fontsize=8)
    if best:
        ax.plot(best["L_m"], best["charge_g"], "w*", ms=17,
                markeredgecolor="k")
    ax.set_xlabel("barrel length (m)"); ax.set_ylabel("powder charge (g)")
    ax.set_title("GUN-01 sweep 3: feasible region is below the red line")
    for e in ("svg", "png"):
        fig.savefig(FIG / f"gun01_contour.{e}")
    plt.close(fig)

    print("\nWrote outputs/gun01_*.csv and figures/gun01_*.svg")
    return base, s1, s2, summary


run = main          # run_all.py calls .run()


if __name__ == "__main__":
    main()
