# COIL-01 — Coilgun Design for the FG01 3 cm Rhenium Canister

**Code:** [`sims/coil01_design.py`](../sims/coil01_design.py) (standalone).
**Outputs:** `outputs/coil01_*.csv`, `figures/coil01_*.svg`.
**Given:** 3 cm bore, **100 g total accelerated mass**, rhenium payload.

---

## 1. Answers up front

| Question | Answer |
|---|---|
| **Is a coilgun optimal for FG01?** | **On performance yes, on mass no — at FG01's realised shot count the chemical gun is lighter.** See §6. |
| **Optimal length** | **1.36 m** at a 50,000 g acceleration limit — or **0.66 m** if the canister is designed for 100,000 g |
| **Capacitor energy to store** | **≈42 kJ** for 400 m/s · **≈95 kJ** for 619 m/s |
| **Capacitor mass** | 28 kg (400 m/s) · 64 kg (619 m/s) at 1500 J/kg |

**The single most useful finding:** stored energy is **essentially invariant at ~40–44 kJ** across a
10× range of acceleration limit. Only the *length* changes. So the capacitor bank is fixed by the
physics and the length is purely a payload-structure decision.

---

## 2. Two facts that constrain the design before anything else

**Rhenium cannot be driven by a coilgun, at all.**

- It is **paramagnetic** (χ ≈ +7×10⁻⁵), so a *reluctance* coilgun has nothing to grip. Using an iron
  armature instead caps magnetic pressure at saturation: B ≤ 2.1 T → B²/2μ₀ = **1.75 MPa**, roughly
  6× too low for this design, before counting the iron's dead mass.
- Its conductivity is poor (5.42×10⁶ S/m, **11× worse than copper**), so it cannot serve as its own
  *induction* armature either.

So the design is an **induction launcher** with a **thin aluminium canister that is simultaneously
the armature and the payload container**, and which flies with the payload. No discarding sabot —
that would create exactly the debris FG01 exists to remove.

**Mass budget, from the 100 g total:**

| Item | Mass |
|---|---|
| Aluminium canister/armature (3 cm dia, 4 cm long, 2 mm wall) | **19.0 g (19%)** |
| **Rhenium payload** | **81.0 g** |

The 19 g armature is **not optional** — it is the thing the field pushes on. Ask for 100 g of
rhenium accelerated and the launch mass becomes 123 g, needing ~23% more stored energy.

---

## 3. Method

Per stage, a capacitor discharges into a drive coil magnetically coupled to the armature, solving
two coupled circuits plus motion:

```
[L_c  M ] [dI_c/dt]   [V_c − R_c·I_c − v·(dM/dx)·I_a]
[M   L_a] [dI_a/dt] = [     −R_a·I_a − v·(dM/dx)·I_c]

m·dv/dt = I_c·I_a·dM/dx        dV_c/dt = −I_c/C
```

with `M(x) = M_peak/(1+(x/R)²)^{3/2}` and `M_peak = k·√(L_c·L_a)`, k = 0.6. RK4, vectorised across
the design grid. The armature sits *ahead* of the coil so dM/dx < 0; Lenz makes I_a oppose I_c, and
the product `I_c·I_a·dM/dx` is positive — repelled forward. Standard coaxial induction geometry.

**Computed armature:** L_a = 19.3 nH, R_a = 29.1 µΩ, **L/R = 0.66 ms**. The decay time is ~7× the
pulse length, so induced current persists through the whole transit — this is what makes the
efficiency respectable.

**Objective: maximise stage *efficiency*, not ΔKE.** This matters. An early version maximised
ΔKE and returned 326,000 µF capacitors storing 168 kJ to deliver 376 J — the right answer to the
wrong question. Capacitor mass is what the spacecraft actually carries.

### Constraints — the coilgun's analogues of chamber pressure

1. **Canister structure.** The base carries the rhenium column: σ = ρ_Re·a·L_column. 81 g in a 3 cm
   bore is a **7.3 mm column**, so at 50,000 g the base sees **75 MPa** against Al-6061's 276 MPa
   yield (**SF 3.7**).
2. **Coil structure.** Magnetic pressure B²/2μ₀ bursts the drive coil exactly as chamber pressure
   bursts a barrel. At 15 T that is **90 MPa** — bare copper (yield 70–200 MPa) cannot hold it, so
   the coil needs steel or composite banding. Above ~20 T this becomes pulsed-magnet engineering
   and is excluded.

---

## 4. Design points

| Target | Achieved | Stages | **Length** | **Energy stored** | Efficiency | Cap mass | Cu mass | Recoil | Charge @5 kW |
|---|---|---|---|---|---|---|---|---|---|
| 300 m/s | 308 | 18 | 0.79 m | **25.1 kJ** | 18.9% | 16.7 kg | 3.2 kg | 30.8 N·s | 5.6 s |
| **400 m/s** | **407** | **31** | **1.36 m** | **42.3 kJ** | **19.6%** | **28.2 kg** | 5.5 kg | 40.7 N·s | 9.4 s |
| 500 m/s | 500 | 44 | 1.94 m | **62.4 kJ** | 20.0% | 41.6 kg | 7.8 kg | 50.0 N·s | 13.9 s |
| **619 m/s** | **623** | **67** | **2.95 m** | **95.3 kJ** | **20.4%** | **63.6 kg** | 11.9 kg | 62.3 N·s | 21.2 s |
| 700 m/s | 701 | 83 | 3.65 m | 119.6 kJ | 20.5% | 79.7 kg | 14.8 kg | 70.1 N·s | 26.6 s |

Overall efficiency **~20%**, which is credible for a multistage induction launcher at this scale
(published large-bore systems reach 20–50%; small ones 15–25%). The first stage is only 11%
efficient because the armature starts from rest — normal, and why real coilguns use an injector.

**619 m/s is expensive.** Going from 400 to 619 m/s costs **2.3× the energy and 2.2× the length**.
Since INT-3 already recommends staying at or below 619 m/s and consumables are negligible, **400 m/s
is the better operating point** unless the mission specifically needs the higher intercept velocity.

---

## 5. The acceleration trade — the key design insight

Target 400 m/s throughout:

| Accel limit | Stages | **Length** | Energy stored | Cap mass | B peak | Base stress | **Canister SF** |
|---|---|---|---|---|---|---|---|
| 20,000 g | 90 | 3.96 m | 40.3 kJ | 26.9 kg | 6.2 T | 30 MPa | 9.2 |
| **50,000 g** | 31 | **1.36 m** | 42.3 kJ | 28.2 kg | 9.9 T | 75 MPa | **3.7** |
| **100,000 g** | 15 | **0.66 m** | 40.5 kJ | 27.0 kg | 13.3 T | 150 MPa | **1.8** |
| 200,000 g | 13 | 0.57 m | 44.0 kJ | 29.3 kg | 14.8 T | 299 MPa | **0.9 ✗ fails** |

> **Stored energy is invariant at 40–44 kJ across a 10× range of acceleration limit.**
> Only the length responds.

Because efficiency stays ~20% regardless, the capacitor bank is fixed by physics. The length is
purely a **payload-structure** decision:

- **50,000 g → 1.36 m, SF 3.7** — comfortable, recommended baseline.
- **100,000 g → 0.66 m, SF 1.8** — halves the launcher, still positive margin, needs a thickened
  canister base (the 150 MPa is a *pressure* on the base plate, so the base must be sized in
  bending, not just compression).
- **200,000 g fails** — the canister base yields.

**Coil turns barely matter:** 15→60 turns moves energy only 39.8→42.3 kJ and length 1.41→1.23 m.
The real effect is drive voltage: **N=22 at 7 kV** is meaningfully easier to switch than N=30 at
12 kV or N=45 at 18 kV, for a ~4% energy penalty. **Recommend N ≈ 22, ~7 kV.**

---

## 6. Is the coilgun optimal? — the honest answer

Against the chemical gun from [GUN_BALLISTICS.md](GUN_BALLISTICS.md), both delivering 100 g at
400 m/s:

| Metric | Coilgun | Chemical gun | Winner |
|---|---|---|---|
| Length | 1.36 m | **0.60 m** | chemical (2.3×) |
| Prime hardware mass | 43.8 kg | **5.4 kg** | **chemical (8×)** |
| Consumable per shot | **0 g** | 12 g propellant + case | coilgun |
| Recoil | **40.7 N·s** | 43.4 N·s | coilgun (6%) |
| Waste heat per shot | **8.3 kJ** | 10.7 kJ | coilgun |
| Cadence limit | **9.4 s** (charging) | 150 s (radiative cooling) | **coilgun (16×)** |
| Barrel erosion | **none** | yes | coilgun |
| Velocity tunable per shot | **yes** | no (fixed charge) | coilgun |
| Moving parts | **none** | breech, autoloader, extractor | coilgun |

**Mass crossover: 1,599 shots.** Below that the chemical gun is lighter; above it the coilgun is.

**And FG01 does not get near that number.** ECO-1 found the platform achieves **~24 engagements per
decade** — the Δv-per-engagement invariant, not the launcher, limits the mission. Even at several
shots per engagement, FG01 fires *tens to hundreds* of shots in its life, **1–2 orders of magnitude
short of the crossover**.

### So: not optimal on mass, but probably still the right choice

The coilgun loses the mass comparison by 38 kg and wins on everything else. Three of those wins are
worth more than 38 kg on a ~600 kg platform:

1. **Velocity tunable shot-to-shot.** INT-3 showed the safe intercept velocity depends on β and
   E_s,c, both uncertain, and on target size. A coilgun dials velocity by choosing how many stages
   to fire and at what voltage. A chemical gun is fixed at load time. Given that **β is the widest
   remaining uncertainty in the whole project** (1.17–1.68), being able to re-tune in flight is a
   genuine risk retirement.
2. **No consumable, no moving parts, no erosion.** A 10-year mission with a breech mechanism,
   autoloader and propellant that outgasses in vacuum is a reliability problem; a coilgun is
   capacitors and switches.
3. **16× better cadence** — 9.4 s versus 150 s. Only matters if engagements ever cluster.

**Recommendation: coilgun, at 400 m/s, ~42 kJ, N≈22 turns at ~7 kV.** Choose length from the
canister's acceleration rating — **1.36 m at 50,000 g** for margin, or **0.66 m at 100,000 g** if
the shorter launcher is worth designing a stiffer canister base.

*The one condition that would flip this:* if FG01's shot count ever exceeds ~1,600 (a much larger
fleet, or many shots per engagement), the chemical gun becomes the lighter system and its
disadvantages would need re-weighing.

---

## 7. Design card

| Parameter | Value |
|---|---|
| Architecture | Multistage **induction** coilgun (not reluctance — rhenium is non-magnetic) |
| Bore | 30 mm |
| Total accelerated mass | 100 g (**19 g Al canister/armature + 81 g rhenium**) |
| Muzzle velocity | 400 m/s |
| **Stored capacitor energy** | **42.3 kJ** |
| **Length** | **1.36 m** (31 stages @ 44 mm pitch) — or 0.66 m at 100,000 g |
| Coil | ~22 turns/stage, 40 mm long, AWG 10, ~7 kV |
| Capacitance | 1,451 µF (stage 1) grading to 13 µF (stage 31) |
| Efficiency | 19.6% |
| Peak acceleration | 49,700 g |
| Peak field | 9.9 T (90 MPa magnetic pressure — coil needs banding) |
| Capacitor mass | 28.2 kg @ 1500 J/kg |
| Copper mass | 5.5 kg |
| Recoil impulse | 40.7 N·s |
| Charge time @ 5 kW | 9.4 s |
| Waste heat per shot | 8.3 kJ |

**Capacitors grade strongly:** 1,451 µF at stage 1 down to 13 µF at stage 31. Pulse length must
match transit time (~2R/v), so as the armature speeds up the pulses must get shorter and sharper.
This is not optional — a uniform bank would waste most of its energy.

---

## 8. Limitations

- **k = 0.6 is assumed, and efficiency is quadratic-ish in it.** The theoretical ceiling for
  transformer-coupled energy transfer is k²/(1+√(1−k²))² — 11% at k=0.6, 25% at k=0.8. That the
  model returns ~20% overall reflects multi-stage re-use and persistent armature current, but a
  measured k is the single most valuable input this design lacks. `UNVALIDATED`
- **Skin effect not modelled.** Current crowds into the armature surface during a fast pulse,
  raising effective R_a and lowering efficiency. **The 20% figure is therefore optimistic**; treat
  it as 15–20%, i.e. capacitor energy 42–56 kJ for 400 m/s.
- **No switch losses, no capacitor ESR, no bus inductance.** Real systems lose another 5–15% here.
- **Ideal triggering.** Each stage fires at its optimal position; real timing jitter costs a few
  percent per stage.
- **Coil heating over a burst not modelled** — 8.3 kJ per shot into copper, with the same
  radiative-only rejection problem the chemical gun has.
- **Capacitor energy density 1500 J/kg** is good-quality pulse film. Space-qualified parts may be
  heavier; at 1000 J/kg the bank is 42 kg rather than 28 kg, which pushes the crossover to ~2,400
  shots and strengthens the chemical gun's mass case.
