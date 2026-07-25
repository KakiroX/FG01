# GUN-01 — Internal Ballistics of a Chemical-Propellant Launcher

**Code:** [`sims/gun01_ballistics.py`](../sims/gun01_ballistics.py) (standalone —
only numpy/pandas/matplotlib). Run with `python sims/gun01_ballistics.py`.
**Outputs:** `outputs/gun01_*.csv`, `figures/gun01_*.svg`.

---

## 1. Executive summary

**400 m/s for a 100 g projectile is achievable inside a 100 MPa chamber-pressure
limit — but not at the 20 mm bore assumed at the outset.** At 20 mm the pressure
cap binds hard: the most that can be extracted is **315 m/s**, and reaching
400 m/s requires 115–132 MPa, which violates the limit at every barrel length
from 0.6 m to 2.0 m. The fix is not more powder or a longer barrel — both are
nearly powerless against the cap — it is **a bigger bore**.

Recommended design: **25 mm bore, 0.6 m barrel, 12.0 g of Retumbo**, giving
**400 m/s at 85.4 MPa peak** with a 15% safety margin against the cap. Bore
diameter is by far the strongest lever: at a fixed 100 MPa ceiling, muzzle
velocity rises from 121 m/s at 12 mm to 984 m/s at 40 mm.

**The design is not limited by ballistics. It is limited by recoil.** Each shot
delivers **43.4 N·s**, which is **4.34 m/s of Δv to a 10 kg spacecraft** — the
launcher would push a smallsat out of its orbit faster than it deorbits debris.
A platform of **≥600 kg** is required to keep per-shot Δv under 0.1 m/s. Peak
acceleration is **41,000 g**, which the payload — not the gun — has to survive.

---

## 2. Methodology

Zero-dimensional (lumped-parameter) interior ballistics, integrated with
**4th-order Runge–Kutta**.

State vector **(x, v, z)** = travel, velocity, burnt web fraction. Pressure is
*not* integrated — it follows algebraically from the energy balance each step,
which is both more accurate and far more stable than integrating dP/dt.

| # | Equation | Form used |
|---|---|---|
| 1 | Burn rate (St. Robert / Vieille) | r = a·Pⁿ |
| 2 | Web consumption | dz/dt = r / e₁ |
| 3 | Form function (grain geometry) | ψ(z) = 1 − (1−z)^1.6 (extruded stick) |
| 4 | Energy balance (**Resal**) | P·V_free = f·m_b·(1−h) − (γ−1)·½·φ·m_p·v² |
| 5 | Free volume (Noble–Abel) | V_free = V₀ + A·x − (m_c−m_b)/ρ_s − b·m_b |
| 6 | Projectile motion | m_p·dv/dt = (P − P_atm)·A − F_friction |

Equation 4 is the first law `m_b·c_v·T = m_b·f/(γ−1) − E_kin − Q_loss` combined
with the Noble–Abel gas law `P·(V − b·m) = m·R·T`. **φ = 1 + m_c/(3·m_p)** is the
Lagrange factor, accounting for the kinetic energy of the gas column itself (gas
at the breech is at rest, gas at the projectile base moves at v).

**Ignition.** St. Robert's law gives r = 0 at P = 0, so a charge cannot light
itself — the model needs a primer. `z_ignite = 0.02` (2% of the web converted
instantaneously) reproduces a realistic ~15–20 MPa magnum-primer spike. Without
it the integration stalls at zero pressure forever.

**Time-step convergence** (reference dt = 5×10⁻⁸ s): velocity error is
**2×10⁻⁶ %** even at dt = 2×10⁻⁶ s. The integration is fully converged; dt is
not a source of uncertainty.

---

## 3. Validation

Burn-rate coefficient `a` and burn distance `e₁` are **not published
quantities**. They are **calibrated** here against each propellant's published
reference load, then held fixed for all space-launcher results.

| Reference load | Published v | Model v | Δ | Published P_max | Model P_max | Δ |
|---|---|---|---|---|---|---|
| .300 Win Mag, 180 gr, 71 gr H1000 | 900 m/s | 905 | **+0.6%** | 441 MPa | 438 | **−0.6%** |
| .308 Win, 168 gr, 44 gr Varget | 795 m/s | 767 | **−3.5%** | 415 MPa | 418 | **+0.6%** |
| .30-06, 165 gr, 57 gr IMR 4350 | 850 m/s | 845 | **−0.6%** | 405 MPa | 405 | **+0.1%** |
| .300 RUM, 200 gr, 82 gr Retumbo | 910 m/s | 923 | **+1.4%** | 441 MPa | 461 | **+4.5%** |

All four agree within **±3.5% on velocity and ±4.5% on pressure**, and
thermodynamic efficiency lands at 23–26%, squarely in the 25–35% band expected
for small arms.

**Two model defects were found and fixed during validation, not papered over.**
A first calibration matched peak pressure but overshot velocity by a uniform
**+10%** across all four powders — a systematic bias, so a model error rather
than a coefficient error. Two causes, each contributing about half:

1. **Form function.** Idealised single-perf tubular grain (ψ = z, perfectly
   neutral) is wrong: real extruded sticks burn on their *ends* as well as the
   perforation, making them mildly degressive. ψ = 1 − (1−z)^1.6 fixed it.
2. **Bore friction.** 8 MPa equivalent was too low for a rifled barrel driving
   an engraved jacketed bullet; 35 MPa closes the gap.

**The 35 MPa figure deliberately does not carry over to the space launcher.**
A smoothbore firing a sabot on a low-friction obturator has no rifling
engraving, so it uses 3 MPa. That difference is physical, not a fudge factor.

> **Note on the stated validation target.** The brief suggested "10 g of H1000 in
> a 0.6 m barrel → 800–900 m/s for an 8–10 g bullet." That load is not
> physically realisable: 10 g of H1000 will not fit in a .30-calibre case (~6 cm³
> holds ~4.6 g), and if it did, chamber pressure would far exceed proof. The
> calibration therefore uses the real .300 Win Mag magnum load (4.6 g H1000,
> 11.66 g bullet, 900 m/s), which is the load the suggested figure appears to
> describe.

---

## 4. Base case — 100 g @ 20 mm, 10 g H1000, 1.0 m barrel

| Parameter | Unit | Value |
|---|---|---|
| Muzzle velocity | m/s | **501** |
| Maximum chamber pressure | MPa | **152.9** ⚠ over the 100 MPa cap |
| Peak pressure location | m | **0.022** (2.2% of barrel travel) |
| Travel time down barrel | ms | **2.93** |
| Kinetic energy at muzzle | J | **12,549** |
| Propellant energy efficiency | % | **28.9** |
| Peak acceleration | g | **48,031** |
| Chamber volume | cm³ | 16.7 |
| Charge burnt at muzzle | — | 0.87 |
| Recoil impulse | N·s | 52.9 |

The base case **exceeds the structural limit by 53%**. Note also that peak
pressure occurs after only 22 mm of travel — the pressure curve is sharply
peaked, which is exactly what wastes structural margin.

*Figure:* `gun01_base_case.svg` — pressure, velocity and position vs. time.

---

## 5. Sensitivity analyses

### Sweep 1 — barrel length (10 g fixed, 20 mm bore)

| L (m) | v (m/s) | P_max (MPa) | t (ms) | Efficiency (%) | Δv per +0.1 m |
|---|---|---|---|---|---|
| 0.5 | 443.6 | 152.9 | 1.88 | 22.6 | — |
| 0.7 | 473.0 | 152.9 | 2.31 | 25.7 | 13.1 |
| **0.9** | **493.1** | 152.9 | 2.73 | 28.0 | **9.2 ← knee** |
| 1.2 | 513.9 | 152.9 | 3.32 | 30.4 | 6.0 |
| 1.5 | 528.1 | 152.9 | 3.90 | 32.1 | 4.2 |
| 2.0 | 543.6 | 152.9 | 4.83 | 34.0 | 2.5 |

**Peak pressure is completely independent of barrel length** (152.93 MPa at every
length) because the peak occurs at 22 mm of travel — long before the muzzle.
Lengthening the barrel buys velocity but does nothing for the binding
constraint. Diminishing returns set in at **0.9 m** (90.7% of the 2.0 m
velocity); beyond ~1.2 m you are carrying barrel mass for a few m/s.

### Sweep 2 — powder charge (0.9 m barrel, 20 mm bore)

| Charge (g) | v (m/s) | P_max (MPa) | Efficiency (%) | Under 100 MPa? |
|---|---|---|---|---|
| 5 | 322.3 | 101.9 | 23.9 | **no (just)** |
| 7 | 399.1 | 123.7 | 26.2 | no |
| 10 | 493.1 | 152.9 | 28.0 | no |
| 15 | 615.8 | 195.8 | 29.1 | no |
| 20 | 711.8 | 233.7 | 29.1 | no |

Velocity scales roughly as √(charge) while pressure scales nearly **linearly** —
so adding powder buys velocity at a steadily worsening pressure cost. At 20 mm,
*no charge in the 5–20 g range* satisfies the pressure limit while making
400 m/s. Efficiency also saturates near 29%: past ~15 g the extra propellant is
mostly unburnt or vents at the muzzle.

### Sweep 3 — combined optimisation, constraint P_max < 100 MPa

*Figure:* `gun01_contour.svg` — velocity contours with the 100 MPa boundary in
red and the 400 m/s requirement dashed white. **The two lines do not intersect
at 20 mm bore.**

Exact charge to hit 400 m/s (bisection), 20 mm bore:

| L (m) | Charge (g) | P_max (MPa) | Peak accel (g) | Under cap? |
|---|---|---|---|---|
| 0.6 | 7.82 | 132.1 | 41,342 | **FAIL** |
| 1.0 | 6.86 | 122.3 | 38,213 | **FAIL** |
| 1.5 | 6.39 | 117.3 | 36,628 | **FAIL** |
| 2.0 | 6.19 | 115.1 | 35,916 | **FAIL** |

Even a 2 m barrel — quadrupling the length — only drops peak pressure from 132
to 115 MPa. **Barrel length is not a usable lever against the pressure cap.**

### The lever that does work — bore diameter

Largest charge that stays under 100 MPa, 0.9 m barrel:

| Bore (mm) | Max charge (g) | v (m/s) | Peak accel (g) | Recoil (N·s) |
|---|---|---|---|---|
| 12 | 0.62 | 121 | 11,187 | 15.1 |
| 16 | 1.98 | 211 | 19,888 | 24.3 |
| 20 | 4.83 | 315 | 31,074 | 35.1 |
| **25** | **11.84** | **463** | 48,554 | 50.4 |
| 30 | 24.73 | 626 | 69,917 | 67.3 |
| 40 | 80.44 | 984 | 124,297 | 104.1 |

**400 m/s becomes reachable under 100 MPa at bore ≥ 25 mm.** The reason is
straightforward: at fixed pressure, force scales with bore *area* (d²), so the
projectile accelerates faster, the chamber volume grows faster, and the pressure
peak is self-limiting. Bore is the only parameter with quadratic leverage.

---

## 6. Recommended design

At 25 mm bore, choosing the propellant that delivers the most **velocity per
unit of peak pressure** — since pressure, not energy, is the binding constraint:

| Propellant | Charge (g) | P_max (MPa) | v per MPa | Burnt at muzzle | Under cap? |
|---|---|---|---|---|---|
| **Retumbo** | 12.0 | **79.9** | **5.01** | 0.47 | **yes** |
| H1000 | 9.2 | 87.3 | 4.58 | 0.56 | yes |
| IMR 4350 | 8.9 | 90.4 | 4.42 | 0.59 | yes |
| Varget | 7.6 | 101.8 | 3.93 | 0.68 | **no** |

The ranking is exactly inverse to burn speed: **slower powder ⇒ flatter pressure
curve ⇒ more velocity per unit peak pressure.** Varget, the fastest, is the only
one that fails the cap.

### Design card

| Parameter | Unit | Value |
|---|---|---|
| Bore diameter | mm | **25** |
| Barrel length | m | **0.6** |
| Powder charge | g | **12.0 (Retumbo)** |
| Chamber volume | cm³ | 20.0 |
| **Muzzle velocity** | m/s | **400** |
| **Maximum chamber pressure** | MPa | **85.4** (15% margin) |
| Peak pressure location | m | 0.015 |
| Travel time down barrel | ms | **2.28** |
| Kinetic energy at muzzle | J | **8,000** |
| Propellant energy efficiency | % | **15.0** |
| Peak acceleration | g | **41,250** |
| Projectile base stress | MPa | 85.4 |
| Inertial compressive stress | MPa | 155.7 |
| Recoil impulse | N·s | **43.4** |

Efficiency drops to 15% because only 43% of the slow Retumbo charge burns before
the projectile exits. **That inefficiency is deliberate and correct** — it is the
price of a flat pressure curve, and propellant is cheap while structural mass
is not.

---

## 7. Structural analysis

Hoop stress σ = P·D/(2t), at 85.4 MPa peak and 25 mm bore:

| Material | Yield (MPa) | SF | Wall (mm) | Barrel mass (kg) |
|---|---|---|---|---|
| 4340 steel | 500 | 2.0 | 4.27 | 1.85 |
| **Ti-6Al-4V** | 900 | 2.0 | **2.37** | **0.54** |
| Inconel 718 | 1030 | 2.0 | 2.07 | 0.87 |
| **CF overwrap (T1000)** | 3000 | 2.0 | **0.71** | **0.055** |

**Best strength-to-weight: a thin Inconel or steel liner (for erosion and heat)
overwrapped with filament-wound carbon fibre.** CF alone gives a 0.71 mm wall at
55 g — 34× lighter than steel — but cannot take the bore heat or erosion, so the
practical answer is a hybrid: ~1 mm Inconel liner + CF overwrap, ~0.4 kg total.

**Projectile stress.** Base pressure is 85.4 MPa, but the *inertial* compressive
stress in a 20 mm tungsten slug at 41,250 g is **155.7 MPa** — nearly double, and
the real design driver. Tungsten (compressive strength >1 GPa) is comfortable;
a fragile payload is not.

### Recoil — the binding constraint

| Spacecraft mass | Δv per shot |
|---|---|
| 10 kg | **4.34 m/s** ⚠ |
| 100 kg | 0.434 m/s |
| **600 kg** | **0.072 m/s** |

43.4 N·s per shot is a large impulse. On a 10 kg platform the launcher imparts
more Δv to *itself* than the 51.8 m/s kick it is trying to deliver to a debris
object — after ~12 shots the platform has changed its own orbit more than the
target's. **A ≥600 kg platform is required**, which aligns with the 602 kg
reference vehicle in SIM-11 and must be corrected by propellant either way.

### Thermal — the multi-shot limit

- Heat into barrel per shot: **10.7 kJ** (20% of chemical energy)
- Ti barrel (0.54 kg): **+37.8 K per shot**
- Radiative rejection at 350 K, ε = 0.8: **71 W**
- **Time to reject one shot's heat: 150 s**

In vacuum there is no convection — radiation only. **Sustained fire is limited to
roughly one shot per 2.5 minutes**, or ~24 shots/hour, before the barrel walks
up in temperature. A 10-shot burst raises the barrel 378 K, which softens
titanium (Ti-6Al-4V loses ~40% of yield strength by 700 K) and would require
re-derating the pressure limit. **This is the strongest argument for a
high-temperature liner** — and it is the practical constraint on engagement
cadence, independent of everything in ECO-1.

---

## 8. Discussion

### Bottlenecks, in order

1. **Recoil (43.4 N·s/shot).** The hard limit. It forces a platform ≥600 kg and
   consumes propellant to correct. Momentum is conserved; there is no design
   trick that avoids it — only a heavier platform or an actively cancelled
   counter-mass.
2. **Barrel heating in vacuum (150 s/shot).** Sets the firing cadence. Radiation
   is the only sink.
3. **Peak acceleration (41,250 g).** The payload must survive it. For a dispersed
   grain cloud in a sabot this is likely fine; for anything with structure it is
   not.
4. **Pressure cap (100 MPa).** Real but solvable — by bore diameter, not by
   powder or barrel length.
5. **Propellant efficiency (15%).** Least important. Propellant is cheap and
   light; do not trade structural margin for it.

### Recommendations

- **Go to 25–30 mm bore.** This is the single highest-leverage change and it is
  what makes the requirement feasible at all.
- **Use the slowest available propellant** with the largest web. Flat pressure
  curve = more velocity per MPa. Retumbo over H1000 over Varget.
- **Shorten the barrel to 0.6 m.** Beyond ~0.9 m you buy a few m/s per 0.1 m at
  real mass cost, and pressure does not improve.
- **Inconel liner + carbon-fibre overwrap.** 0.4 kg vs 1.85 kg for steel, and the
  liner handles the heat the CF cannot.
- **Do not fly this on a smallsat.** Below ~600 kg the recoil dominates the
  mission.

### Comparison with the alternatives already studied

| | Chemical gun (this study) | Coilgun (SIM-6) | Cold gas |
|---|---|---|---|
| Muzzle energy, design shot | 8,000 J | 410 J | ~10 J |
| Prime hardware vs. demonstrated | inside (20 mm/25 mm guns are routine) | **80,491× inside** | inside |
| Consumable per shot | 12 g propellant + case | **electricity only (~$4×10⁻⁵)** | gas |
| Recoil per shot | **43.4 N·s** | 1.3 N·s | low |
| Shots per reload | limited by propellant mass | **11,152 per 30 kg** | limited |
| Cadence in vacuum | **1 per 150 s (thermal)** | 0.26 s charge at 5 kW | fast |
| Barrel/rail wear | erosion per shot | negligible | none |

**For FG01's actual design shot (2.14 g at 619 m/s = 410 J), the coilgun is
clearly superior** — 80,491× inside demonstrated hardware, no consumable, no
erosion, 33× less recoil, and a cadence limited by charging rather than by
radiative cooling. A chemical gun only becomes attractive if the required muzzle
energy rises by ~20× (to the 8 kJ studied here), where coilgun capacitor mass
starts to dominate.

That regime — 100 g projectiles at several hundred m/s — corresponds to engaging
**large** objects, which is precisely the window ECO-5 identified (100 kg to
3.7 t targets). **The chemical gun is the right launcher for the wrong half of
the FG01 concept, and the coilgun for the right half.**

---

## 9. Engineering questions

**Which propellant?** **Retumbo**, or slower. Quantitatively it delivers 5.01 m/s
of muzzle velocity per MPa of peak pressure against H1000's 4.58 and Varget's
3.93 — a 27% structural-mass advantage over Varget for the same performance.
Secondary but important: Retumbo and H1000 are Hodgdon *Extreme* series, which
are **temperature-insensitive** — decisive in orbit, where a barrel can swing
from 200 K in eclipse to 350 K in sun and a temperature-sensitive powder like
IMR 4350 would change both velocity and peak pressure shot to shot. Caveat:
nitrocellulose propellants contain volatiles and absorbed moisture that
**outgas in vacuum**, changing burn rate over months. Hermetic cartridge sealing
is mandatory, and this is an unaddressed reliability risk.

**Which barrel material?** **Inconel 718 liner (1–2 mm) + T1000 carbon-fibre
overwrap.** CF has ~30× the specific strength of steel (1.9 MJ/kg vs 0.064) and
gets the barrel to 55 g of pressure-bearing structure, but it cannot take bore
heat or erosion. The liner does that. Monolithic Ti-6Al-4V (0.54 kg) is the
simple fallback but loses ~40% of its yield strength by 700 K, which the thermal
analysis says is reachable in a 10-shot burst.

**Sabot design?** A discarding sabot is **disqualified**: it creates exactly the
debris this system exists to remove. Use a **non-discarding pusher/obturator**
that stays with the payload, or a captured design retained by the barrel. Best
material is a PEEK or polyimide obturating band on an aluminium carrier — PTFE
has the lowest friction but cold-flows and outgasses. In vacuum, sliding friction
is *higher* than in air for many polymers (no adsorbed moisture film), so the
3 MPa friction figure used here is optimistic and should be measured.

**Multi-shot heating?** Covered in §7: 10.7 kJ per shot, 71 W radiative
rejection, **150 s per shot** to hold temperature. Ten rapid shots put +378 K
into a titanium barrel and force a pressure derate. Mitigations: high-temperature
liner, a dedicated radiator, or accepting a low duty cycle. Barrel erosion also
scales with flame temperature, so a cooler-burning propellant helps twice.

**Safe pressure limit?** **100 MPa is a reasonable and slightly conservative MEOP
for a lightweight space-qualified barrel.** With a CF overwrap at SF 2.0 the wall
is 0.71 mm and mass is negligible, so the limit could be raised to 200–300 MPa
for a modest mass penalty — but there is no reason to: at 25 mm bore, 85 MPa
already meets the requirement, and margin against a fatigue-and-erosion-driven
failure mode is worth more than velocity the mission does not need.

---

## 10. Assumptions and non-linear effects

**Stated assumptions:** 0-D lumped parameters with uniform chamber pressure
(Lagrange factor applied as the first-order correction); heat loss a fixed 20% of
released chemical energy; perfect obturation, zero leakage; no barrel erosion or
thermal soak within a shot; γ = 1.23 constant; covolume 1.0×10⁻³ m³/kg constant;
friction and engraving lumped into a constant resistive pressure; vacuum ahead of
the projectile; ignition transient not modelled beyond a 2% primer.

**Non-linear effects present in the model:** the P^n burn-rate feedback loop
(pressure raises burn rate raises pressure) is the dominant non-linearity and is
what makes peak pressure so sensitive to charge and web thickness — it is why
pressure scales nearly linearly with charge while velocity scales as √(charge).

**Non-linear effects NOT modelled, and where they would bite:**
- **Erosive burning** — gas flow over grains raises local burn rate at high
  loading density. Would *raise* peak pressure above the prediction; the 15%
  margin at the design point partly covers this.
- **Pressure-wave / ignition dynamics** — longitudinal waves in the chamber can
  produce local pressures well above the lumped mean. The classic cause of
  detonation in real guns and not captured by any 0-D model.
- **Grain fracture** at high loading density — suddenly increases burning surface
  and can spike pressure.
- **Shot-to-shot thermal feedback** — a hot barrel raises propellant temperature
  and burn rate, walking pressure up over a burst. Given the 37.8 K/shot heating,
  this is a real effect for anything beyond single shots.

**Model uncertainty to carry forward:** ±3.5% velocity and ±4.5% pressure against
the validation set, plus the unquantified erosive-burning and pressure-wave
margins above. The 85.4 MPa design point should be read as **85 ± 10 MPa**, which
still clears the 100 MPa cap but not by much.
