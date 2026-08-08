# Launcher Trade Study — Methodology

**FG01 / KIDR · momentum-transfer debris remediation**
Chemical-propellant and electromagnetic launcher options for the 100 g rhenium
canister.

Document 1 of 3 · [Results](LAUNCHER_RESULTS.md) · [Conclusions](LAUNCHER_CONCLUSIONS.md)

---

## 1. Scope

### 1.1 Question

FG01 must impart a retrograde Δv of 48.7–57.7 m/s (200 km orbit lowering,
400–1200 km altitude band) to debris by launching a rhenium grain cloud at it.
This study determines **which launcher technology to use**, and sizes it.

Requirement, inherited from the completed v5/v6 engagement analysis:

| Quantity | Value | Origin |
|---|---|---|
| Launched mass | 100 g total | design constraint |
| Bore | 30 mm | canister geometry |
| Muzzle velocity | 619 m/s design, 400 m/s alternate | v6 SIM-13 / CHK-01 |
| Muzzle energy | 19.2 kJ (619 m/s) · 8.0 kJ (400 m/s) | `DERIVED` |
| Chamber-pressure limit | 100 MPa (initial assumption, tested in §Results) | design brief |

### 1.2 What was modelled

| Module | Device | Status |
|---|---|---|
| **GUN-01** | Chemical propellant, 0-D internal ballistics | modelled, validated |
| **GUN-02** | Barrel length, wall taper, chamber, cartridge | modelled |
| **COIL-01** | Multistage **induction coilgun**, coupled-circuit | modelled |
| **COIL-02** | Coilgun structural walls | modelled |
| **COIL-03** | Coilgun radial stack-up | modelled |
| **CHK-01** | Velocity consistency against v5/v6 engagement gates | modelled |
| **CART-01** | Cartridge charge law (rhenium mass vs. target) | modelled |

### 1.3 ⚠ What was NOT modelled — the railgun

**A railgun is a different device from a coilgun and was not simulated.** The
distinction is load-bearing and is stated explicitly to avoid mislabelling:

- **Coilgun (modelled).** Drive coils induce eddy currents in a conducting
  armature. Force is `I_c·I_a·dM/dx`. **No sliding electrical contact.**
- **Railgun (not modelled).** Current passes through the projectile via two
  sliding rail contacts. Force is `½·L′·I²` with L′ the inductance gradient.

A scoping calculation (not a simulation) indicates a railgun is the wrong device
at this operating point:

| L′ | Peak current for 60.9 kN |
|---|---|
| 0.30 µH/m | 637 kA |
| 0.45 µH/m | 520 kA |
| 0.60 µH/m | 451 kA |

Delivering **19.2 kJ** of muzzle energy would require **~0.5 MA through sliding
contacts**. For comparison the US Navy EMRG used 3–5 MA to deliver 32 MJ — three
orders of magnitude more energy for only ~8× the current, because railguns are
efficient at km/s velocities and inefficient at hundreds of m/s. Rail erosion at
these currents limits published rail life to tens–hundreds of shots.

**Conclusion of the scoping check:** a railgun buys nothing a coilgun does not,
and adds a sliding-contact erosion mechanism. It is excluded on these grounds and
is **not** claimed to have been quantitatively evaluated. `SCOPING ONLY`

### 1.4 Traceability convention

Inherited from the v5 plan. Every quantity carries:

- `DERIVED` — follows from stated physics
- `SOURCED` — specific external reference
- `CALIBRATED` — fitted in this study against published data (§4)
- `UNVALIDATED` — assumption with no source at the relevant scale

---

## 2. Chemical launcher model (GUN-01)

### 2.1 Formulation

Zero-dimensional (lumped-parameter) interior ballistics. State vector
**(x, v, z)** = travel, velocity, burnt web fraction. **Pressure is not
integrated** — it follows algebraically from the energy balance at each step,
which is more accurate and far more numerically stable than integrating dP/dt.

| # | Physics | Equation | Tag |
|---|---|---|---|
| 1 | Burn rate (St. Robert / Vieille) | `r = a·Pⁿ` | `SOURCED` |
| 2 | Web consumption | `dz/dt = r/e₁` | `DERIVED` |
| 3 | Form function | `ψ(z) = 1 − (1−z)^1.6` | `CALIBRATED` |
| 4 | Energy balance (Resal) | `P·V_free = f·m_b·(1−h) − (γ−1)·½·φ·m_p·v²` | `SOURCED` |
| 5 | Free volume (Noble–Abel) | `V_free = V₀ + A·x − (m_c−m_b)/ρ_s − b·m_b` | `SOURCED` |
| 6 | Projectile motion | `m_p·dv/dt = (P − P_atm)·A − F_friction` | `DERIVED` |

Equation 4 combines the first law `m_b·c_v·T = m_b·f/(γ−1) − E_kin − Q_loss` with
the Noble–Abel gas law `P·(V − b·m) = m·R·T`. **φ = 1 + m_c/(3·m_p)** is the
Lagrange factor, accounting for the kinetic energy of the propellant gas column
(gas at the breech is at rest, gas at the projectile base moves at v).

### 2.2 Ignition

St. Robert's law gives `r = 0` at `P = 0`, so a charge **cannot ignite itself** —
the model requires a primer. `z_ignite = 0.02` (2% of web converted
instantaneously) reproduces a realistic ~15–20 MPa magnum-primer spike.

> *Recorded because it is a genuine modelling trap:* without this term the
> integration stalls at zero pressure indefinitely rather than failing visibly.

### 2.3 Numerical method

- **4th-order Runge–Kutta**, fixed step
- Adaptive termination: muzzle exit, demise, or stall detection
- Muzzle state obtained by interpolation onto `x = L`, not by the last step

---

## 3. Electromagnetic launcher model (COIL-01)

### 3.1 Device selection

Two facts constrain the choice before any modelling:

1. **Rhenium is paramagnetic** (χ ≈ +7×10⁻⁵). A *reluctance* coilgun has nothing
   to grip. An iron armature caps magnetic pressure at saturation:
   `B ≤ 2.1 T → B²/2μ₀ = 1.75 MPa`, ~6× too low. `SOURCED`
2. **Rhenium conducts poorly** (5.42×10⁶ S/m, 11× worse than copper), so it
   cannot serve as its own induction armature. `SOURCED`

**Therefore: induction launcher with an aluminium canister that is simultaneously
armature and payload container**, flying with the payload. A discarding sabot is
excluded — it would create exactly the debris FG01 removes.

### 3.2 Formulation

Two coupled circuits plus motion, per stage:

```
λ_c = L_c·I_c + M(x)·I_a          λ_a = M(x)·I_c + L_a·I_a
V_c = R_c·I_c + dλ_c/dt           0 = R_a·I_a + dλ_a/dt
m·dv/dt = I_c·I_a·dM/dx           dV_c/dt = −I_c/C
```

solved as the 2×2 system

```
[L_c  M ] [dI_c/dt]   [V_c − R_c·I_c − v·(dM/dx)·I_a]
[M   L_a] [dI_a/dt] = [     −R_a·I_a − v·(dM/dx)·I_c]
```

Mutual inductance uses the coaxial-dipole profile
`M(x) = M_peak/(1+(x/R)²)^{3/2}`, `M_peak = k·√(L_c·L_a)`, k = 0.60
`UNVALIDATED`.

**Geometry.** The armature sits *ahead* of the coil, so `dM/dx < 0`; Lenz makes
`I_a` oppose `I_c`, and `I_c·I_a·dM/dx > 0` — repelled forward. Standard coaxial
induction launcher.

### 3.3 Optimisation objective

**Stage efficiency `ΔKE/E_stored` is maximised, not ΔKE.**

> *Recorded because the alternative was tried and was wrong:* maximising ΔKE
> drives the optimiser to 326,000 µF capacitors storing 168 kJ to deliver 376 J —
> the right answer to the wrong question. Capacitor mass is what the spacecraft
> carries.

Per stage the search grids over capacitance, voltage and trigger position,
vectorised across the design space, subject to the constraints in §3.4.
Capacitance is matched to transit time (`C ∝ 1/v²`), so later stages
automatically grade to shorter, sharper pulses.

### 3.4 Constraints — the coilgun's analogues of chamber pressure

| Constraint | Value | Basis |
|---|---|---|
| Canister structure | 50,000 g | `σ = ρ_Re·a·L_column`; 75 MPa vs Al-6061 276 MPa, SF 3.7 |
| Coil structure | B ≤ 15 T | magnetic pressure `B²/2μ₀`; 15 T = 90 MPa needs banding |

---

## 4. Calibration procedure (GUN-01)

Burn-rate coefficient `a` and burn distance `e₁` are **not published
quantities**. They were fitted against each propellant's published reference
load, then **held fixed** for all space-launcher results.

### 4.1 Reference loads

| Cartridge | Projectile | Charge | Barrel | Published v | Published P_max |
|---|---|---|---|---|---|
| .300 Win Mag | 11.66 g | 4.60 g H1000 | 0.610 m | 900 m/s | 441 MPa |
| .308 Win | 10.89 g | 2.85 g Varget | 0.610 m | 795 m/s | 415 MPa |
| .30-06 | 10.69 g | 3.69 g IMR 4350 | 0.610 m | 850 m/s | 405 MPa |
| .300 RUM | 12.96 g | 5.30 g Retumbo | 0.660 m | 910 m/s | 441 MPa |

### 4.2 Two model defects found during calibration

A first fit matched peak pressure but **overshot velocity by a uniform +10% on
all four powders**. A systematic bias across independent cases indicates a model
error, not a coefficient error. Two causes, each contributing about half:

1. **Form function.** An idealised single-perf tubular grain (ψ = z, perfectly
   neutral) is wrong — real extruded sticks burn on their *ends* as well as the
   perforation, making them mildly degressive. `ψ = 1 − (1−z)^1.6` corrects it.
2. **Bore friction.** 8 MPa equivalent was too low for a rifled barrel driving an
   engraved jacketed bullet; **35 MPa** closes the gap.

**The 35 MPa friction deliberately does not carry to the space launcher.** A
smoothbore firing a sabot on a low-friction obturator has no rifling engraving
and uses **3 MPa**. This difference is physical, not a fitting parameter.

### 4.3 Calibrated coefficients `CALIBRATED`

| Propellant | a (m/s per MPaⁿ) | n | e₁ (m) | Geometry |
|---|---|---|---|---|
| H1000 | 5.00×10⁻³ | 0.80 | 4.20×10⁻⁴ | stick |
| Retumbo | 3.00×10⁻³ | 0.80 | 3.00×10⁻⁴ | stick |
| IMR 4350 | 4.50×10⁻³ | 0.82 | 3.80×10⁻⁴ | stick |
| Varget | 3.50×10⁻³ | 0.83 | 2.60×10⁻⁴ | stick |

---

## 5. Validation protocol

Following the v5 convention, **validation precedes results** and is reported
whether or not it flatters the model.

| Test | Method | Acceptance |
|---|---|---|
| GUN-01 external | 4 published cartridges, velocity and peak pressure | ±5% |
| GUN-01 numerical | Time-step convergence, dt = 2×10⁻⁶ → 5×10⁻⁸ s | <0.1% |
| GUN-01 physical | Thermodynamic efficiency vs. small-arms range | 25–35% |
| COIL-01 internal | Winding mean radius from COIL-03 stack vs. COIL-01 assumption | ±2% |
| COIL-01 external | Overall efficiency vs. published induction launchers | 15–50% |
| CHK-01 | Velocity choice against all v5/v6 engagement gates | pass/fail per gate |
| CART-01 | Charge law vs. INT-3/INT-4 outputs | exact |

---

## 6. Assumptions register

### 6.1 Chemical launcher

| # | Assumption | Tag |
|---|---|---|
| A1 | 0-D lumped parameters, uniform chamber pressure (Lagrange factor as first-order correction) | `DERIVED` |
| A2 | Heat loss fixed at 20% of released chemical energy | `UNVALIDATED` |
| A3 | Perfect obturation, zero leakage | `UNVALIDATED` |
| A4 | No barrel erosion or thermal soak within a shot | `UNVALIDATED` |
| A5 | γ = 1.23 constant | `SOURCED` |
| A6 | Covolume b = 1.0×10⁻³ m³/kg constant | `SOURCED` |
| A7 | Friction lumped into constant resistive pressure | `CALIBRATED` |
| A8 | Vacuum ahead of projectile (P_atm = 0) | `DERIVED` |
| A9 | Ignition modelled only as a 2% primer conversion | `UNVALIDATED` |

### 6.2 Electromagnetic launcher

| # | Assumption | Tag |
|---|---|---|
| B1 | Lumped L-R per circuit, no skin-depth redistribution | `UNVALIDATED` |
| B2 | Dipole M(x), peak coupling k = 0.60 | `UNVALIDATED` |
| B3 | Rigid armature, no melt or deformation | `UNVALIDATED` |
| B4 | Ideal switching, no crowbar, current allowed to ring | `UNVALIDATED` |
| B5 | Room-temperature resistivity, constant | `UNVALIDATED` |
| B6 | No stage-to-stage magnetic interaction | `UNVALIDATED` |
| B7 | Capacitor energy density 1500 J/kg | `UNVALIDATED` |

### 6.3 Non-linear effects **not** modelled

**Chemical:** erosive burning (raises peak pressure), pressure-wave/ignition
dynamics (the classic cause of detonation, invisible to any 0-D model), grain
fracture at high loading density, shot-to-shot thermal feedback.

**Electromagnetic:** skin effect in the armature (raises effective R_a — **so the
reported efficiency is optimistic**), armature heating over a burst, flux
diffusion through the canister wall.

Both sets of omissions bias in the **unfavourable** direction, so reported margins
should be treated as upper bounds.

---

## 7. Reproducibility

| | |
|---|---|
| Language | Python 3.11.9 |
| Dependencies | numpy 2.0.2, pandas 2.2.3, matplotlib 3.10.0 |
| Random seed | 20260723 (fixed) |
| Modules | `sims/gun01_ballistics.py`, `gun02_barrel.py`, `coil01_design.py`, `coil02_structure.py`, `coil03_stackup.py`, `chk01_velocity.py`, `cart01_charge_law.py` |
| Outputs | `outputs/gun01_*`, `gun02_*`, `coil01_*`, `coil02_*`, `coil03_*`, `chk01_*`, `cart01_*` |
| Execution | Each module runs standalone: `python sims/<module>.py` |

`gun01_ballistics.py` and `coil01_design.py` are registered in `run_all.py`; the
later additive modules are not, so the primary v5/v6 pipeline outputs are
unaffected by this study.

**Independence.** The launcher modules import `fg01/` read-only and write only
their own namespaced outputs. No v5 or v6 simulation, output or document was
modified.
