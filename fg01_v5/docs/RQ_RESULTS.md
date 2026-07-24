# FG01 / KIDR — Results by Research Question

Organised in the order of *FG01 techs.md*. Every answer cites the SIM/INT/ECO
module that produced it and carries a traceability tag. Ranges, not point values,
wherever the result depends on β, price, launch cost, τ, Δh or platform
assumptions.

> **v6 update.** Answers below were written from the v5 run. Where v6's resolved
> interaction model or engagement economics changed an answer, the change is
> marked **[v6]** inline. The three material changes are: **β = 1.20** rather
> than 1.5 (raising required mass 1.25×), **cost per object $0.9M–4.2M** rather
> than ~$5,000 (because a platform achieves ~24 engagements, not 10⁴), and the
> overall verdict, which is now **not viable for small debris but competitive
> against large tracked objects**. See [V6_RESULTS.md](V6_RESULTS.md).

---

## 1. Prove that the momentum rhenium cloud is viable for space-debris removal

**Answer: the momentum-kick cloud is viable as physics and not viable as a
small-debris business.** All four physics gates pass at the selected operating
point with probability ≈1. P(viable) = **0.666** at the optimum (1 cm debris,
900 km, Δh = 200 km, v_rel = 619 m/s), rising to **0.997** once one platform
services ≥10,000 objects and falling to **0** below ~3,000.
`SIM-0, SIM-12, SIM-13` `DERIVED`

**[v6] The engagement rate has now been derived, and it is ~24 per
platform-decade, not 10⁴.** ECO-1's invariant Δv/engagement = π·v/(3.5·Ω̇·T)
gives ~107 m/s per engagement and cannot be tuned — widening the altitude band
raises opportunities and hop cost identically. Cost per 1 cm object is therefore
**$0.9M–4.2M**, and the small-debris architecture **FAILS**. The same kick
applied to a large tracked object costs **$5.6M–20M** against $10M–100M for an
ADR tug and is **competitive**. `ECO-1, ECO-4` `DERIVED`

### 1.1 Dust cloud or grains? How small?

**Discrete grains of 0.1–1 mg (0.36–0.45 mm rhenium).** Not dust, not a single
slug. Bounded above by delivery statistics — grain mass must be small enough that
≥100 grains strike the target, giving Poisson CV ≤10% (1 mg → 161 grains on a 1 cm
target, CV 7.9%). Bounded below by conjunction exposure: hazard during the
released mass's descent scales linearly with grain count, so 0.01 mg dust is 100×
worse on that axis. `SIM-3, SIM-10` `DERIVED`

**[v6] Confirmed by the resolved model, with the cloud size now specified:**
450 µm grains through a **5 cm cloud** put 197 grains on a 1 cm target from a
1.75 g shot. The Δv distribution across shots is dominated by aim scatter, not by
Poisson lumpiness — at 200+ grains the swarm averaging is effectively
deterministic. `INT-1, INT-3` `DERIVED`

### 1.1.2 If a cloud, how to predict its behaviour?

Exact Clohessy–Wiltshire relative motion, evaluated with no assumed cone angle
(AM-14). **The finding is that CW is not the operative physics here:** at
engagement times-of-flight (≤0.25 s for R ≤ 100 m), n·t ≤ 2.8×10⁻⁴ and the exact
CW solution differs from straight-line ballistic spreading by 2.7×10⁻⁸. The
pattern size is set by launcher transverse dispersion, σ = δv_T·R/v_rel — a design
variable, not an orbital-mechanics outcome. `SIM-3` `DERIVED`

### 1.2 Can it slow debris enough to reenter, and by how much?

**Yes, and the honest measure is the lifetime cut, not the reentry.** A 200 km
kick imparts 48.7–57.7 m/s across 400–1200 km and reduces orbital lifetime by a
factor of **16.7× at 400 km falling to 2.2× at 1100 km**. For 1 cm debris at
900 km: **63.5 yr → 18.2 yr**, clearing the 25-yr rule. `SIM-0, SIM-2`
`DERIVED` + `SOURCED` (NRLMSISE-00)

### 1.3 How high an orbit can we target?

**Single kick, 25-yr rule, nominal solar:** 900 km (0.5–1 cm), 800 km (2 cm),
700 km (5–10 cm). **5-yr rule:** 800/700/700/600/600 km.
**Staged kicks lift the ceiling to 1200 km** at 2–3 kicks, each independently
sub-catastrophic — so the ceiling is economic, not physical.
**But note the lower bound too:** below ~700 km the kick has no remediation value
because the debris decays inside 25 years unaided. The useful band is a diagonal:
larger debris → lower optimal altitude. `SIM-2, SIM-12` `SOURCED`

### 1.4 At what angle should the launch occur?

**Pure retrograde in the target frame.** FG01 leads the target in-plane and fires
rearward along the velocity vector. Off-axis error θ costs cos θ of useful impulse
— a 10° error costs 1.5%, negligible beside the areal-density effect of the same
error. `SIM-4` `DERIVED`

### 1.4.1 What guarantees the accuracy of the path?

**Nothing guarantees it; the shot pattern is sized so that accuracy is not
required to be perfect.** Error budget at the design point: 5.775 mm total,
composed of range-knowledge timing (5.41 mm), launcher pointing (2.00 mm),
muzzle-speed dispersion (0.162 mm), optical position knowledge (0.055 mm).
Concrete requirement: ranging bandwidth **≥100 MHz** (the plan's 1 MHz gives
σ_R = 33.5 m → 0.54 m of miss, inadequate). `SIM-4` `DERIVED` + `UNVALIDATED`
(pointing, muzzle repeatability are engineering estimates)

### 1.5 Estimate the rhenium required, f(x) where x = height

**f(h) = γ·Δv_shift(h)/(β·v_rel), independent of target mass.** At v_rel = 619 m/s,
β ∈ [1.0, 2.5]:

| h (km) | g Re on target per g Al |
|---|---|
| 400 | 0.075 – 0.187 |
| 600 | 0.071 – 0.178 |
| 800 | 0.068 – 0.171 |
| 1000 | 0.066 – 0.164 |
| 1200 | 0.063 – 0.157 |

**Launched mass is larger** by the areal penalty π·e·d_miss²/A_target — 16.7× at
the design point, giving **2.26 g launched to put 135 mg on a 1 cm target**.
`SIM-1, SIM-3` `DERIVED`

**[v6] Revised with the derived β and the resolved cloud:** β = 1.201 rather than
1.5 raises mass on target to **197 mg**, while a Gaussian cloud is 3× more
mass-efficient than the uniform disc assumed, so launched mass is **1.75 g** —
the two corrections nearly cancel. At 95% per-shot reliability the sizing is
**2.14 g**. `INT-2, INT-3` `DERIVED`

### 1.6 What must the velocity and mass ratio be?

m_p/m_t = γ·Δv_shift/(β·v_rel). The two are exchangeable at constant momentum, but
**not** at constant damage: E_s = γ·Δv·v_rel/(2β) rises linearly with v_rel, so the
energy gate sets a ceiling (1,039–1,231 m/s at β=1.5, γ=2, E_s,c=40) and the
self-disposal requirement sets a floor (87–290 m/s). The design sits at 619 m/s,
inside both. `SIM-0, SIM-1` `DERIVED`

### 1.6.1 Does the ratio depend on distance?

**No.** f(h) depends on altitude (through Δv_shift) and on v_rel and β. Engagement
range R does not enter the momentum balance at all — it enters only the *miss
distance*, and therefore the launched-mass penalty (quadratically). `SIM-1, SIM-3`
`DERIVED`

### 1.7 How precise will the shot be?

**5.775 mm (1σ) at 10 m range**, dominated by ranging and pointing, and
**independent of control latency from 0.1 µs to 100 ms**. `SIM-4` `DERIVED`

### 1.7.1–1.7.2 Does the error matter, and how much does it cut efficiency?

**It matters as mass, not as failure.** Efficiency η_eff = 1/penalty =
A_target/(π·e·d_miss²) = **6.0%** at the design point. Doubling the miss distance
quadruples the launched mass. Because consumables are ~$31 of a ~$5,031 per-object
cost, this is affordable at the design point and becomes dominant only for miss
distances above ~10 cm on centimetre targets. `SIM-3, SIM-4, SIM-9` `DERIVED`

### 1.7.3 Will rhenium behaviour influence removal accuracy?

Yes, through two channels, both quantified: the pattern σ must be matched to the
miss distance (σ\* = d_miss/√2), and grain count must keep Poisson variability of
the delivered momentum below ~10%. `SIM-3` `DERIVED`

### 1.8 Will it create more space junk?

**No. Net = +1 object removed per engagement.** The impact is sub-catastrophic
(21.8 kJ/kg = 0.55× threshold), so the target is cratered, not disrupted;
energy-based cratering gives 0.167 g of ejecta (11.8% of target mass) with 5–50×
the parent's A/m, so ejecta decay faster than the parent. Applying the mandated
NASA SBM non-catastrophic branch verbatim gives 0.26 fragments >1 cm — physically
impossible from an undisrupted 1 cm target, and reported as the conservative bound
with the extrapolation flagged. `SIM-5` `SOURCED` + `DERIVED`

**[v6] Re-asked of the *resolved* cloud — 197 separate craters rather than one
idealised impactor — and the answer strengthens.** The largest single crater is
**1.39 mm**, so no ejecta fragment can approach 1 cm and the engagement cannot
produce a trackable fragment **by crater geometry alone**, without needing the
sub-catastrophic argument at all. Grains embed (450 µm Re penetrates 1.46 mm into
Al, 3.2× its own diameter), so the rebound branch is bounded at zero. Net remains
**+1**. `INT-5` `DERIVED`

### 1.8.1 Will the transporter burn up with the missed rhenium?

**No.** FG01 and the released grains are separate objects on entirely different
trajectories: the grains are on immediate-reentry ellipses, FG01 remains in its
operating orbit. `SIM-11` `DERIVED`

### 1.8.2 Deorbit FG01 or return to a hub to reload?

Direct deorbit from 800 km costs **39 kg** of propellant (600 kg dry, I_sp 300 s);
returning to a 400 km hub costs less below ~500 km. A 30 kg magazine holds
**11,152 design shots**, so reload frequency is not the driver — engagement rate
is. `SIM-11` `DERIVED`

### 1.8.3 Will missed rhenium burn up in the atmosphere?

**Yes, and it reenters within one orbit.** Fired retrograde at v_rel ≥ Δv_direct(h)
the missed mass has a perigee below the atmosphere; from 800 km at 619 m/s it
enters at 7,666 m/s and −8.77° on its first perigee pass (~45 min). Grains of
0.05–1 mm — the entire design range — ablate >99.99% and fully demise.
`SIM-2, SIM-8` `DERIVED` + `SOURCED`

---

## 2. How cost-efficient is the method?

**[v6 — this answer changed materially.]** v5 reported **$531–$50,031 per
object**, conditional on 10³–10⁵ engagements per platform lifetime. ECO-1 has
since derived that number: **~24 per platform-decade**, giving **$0.9M–4.2M per
object** for 1 cm debris. The v5 figures remain correct arithmetic conditional on
an engagement count that turns out to be unreachable.

Per kg: **$6.4×10⁸/kg** for 1 cm debris at the realised engagement count. The
areal penalty still scales as 1/A_target, which is why the same platform reaches
**$1,038–18,954/kg** against a large tracked object. `SIM-9, ECO-1, ECO-4`
`DERIVED` + `UNVALIDATED` (platform cost)

### 2.1 Comparison with competitors

| method | USD/kg | USD/object | addresses 1–10 cm? |
|---|---|---|---|
| Ground-based laser ablation | 100 – 500 (theoretical) | ~$0.14–0.71 | **yes** |
| Electrodynamic tether | 1,000 – 5,000 | — | no |
| Aerogel/foam capture | 5,000 – 20,000 | — | no |
| ADR tug (demonstrated class) | 10,000 – 100,000 | **$10M – 100M** | no |
| **FG01 vs 1 cm debris [v6]** | 6.4×10⁸ | **$0.9M – 4.2M** | yes |
| **FG01 vs large tracked object [v6]** | **1,038 – 18,954** | **$5.6M – 20M** | no |

**Honest verdict:** against small debris, the only genuine competitor is
ground-based laser ablation and FG01 is **5–6 orders of magnitude more expensive
per object**, because its hardware must be in orbit *and* must be moved to each
target. The per-kg comparison with tethers and tugs is not meaningful for small
debris — those methods cannot engage this size class at all.

**[v6] Against large tracked objects the comparison inverts and becomes
favourable:** FG01 beats the demonstrated-class ADR tug on both cost per object
and cost per kg. `SIM-9, ECO-4`

### 2.2 Electrical energy per maximum-capacity shot

**515 J muzzle / 1,288 J stored** at η = 0.4 for the 2.69 g design shot;
33.5 kJ muzzle / 83.8 kJ stored for a 175 g pattern. Electricity cost
**$3.6×10⁻⁵ per shot**. `SIM-6` `DERIVED`

### 2.2.1 Coilgun dimensions — coils, copper thickness

**1.95 m barrel, 7 stages** (10,000 g acceleration limit, 0.30 m stages),
**1.82 kg copper** (AWG 12, 200 turns/m, 5 cm bore), **0.86 kg capacitors**
(1.5 J/g pulse-discharge). `SIM-6` `DERIVED` + `UNVALIDATED` (capacitor energy
density)

### 2.2.2 Energy to charge for maximum velocity

1,288 J, **0.26 s at 5 kW**. The launcher is 3–5 orders of magnitude inside the
demonstrated EM-launch envelope (33 MJ, US Navy EMRG 2010) — the reverse of v4's
position. `SIM-6` `SOURCED`

### 2.3 USD/kg of trash deorbited

$3.56M/kg (1 cm) to $3,756/kg (10 cm) at $50M platform / 10⁴ engagements.
**Cost is 99.4% platform amortisation** — consumables are $31 per object.
`SIM-9` `DERIVED` + `UNVALIDATED`

### 2.3.1 Weight of FG01

**602 kg dry, 632 kg wet** (250 kg coilgun, 30 kg magazine, 50 kg avionics,
167 kg power at 5 kW/30 W·kg⁻¹, 30 kg thermal, 20% structure). `SIM-11`
`UNVALIDATED` (parametric subsystem fractions)

### 2.3.2 Cost/weight to send to space

$3,000/kg (dedicated) – $7,000/kg (rideshare marginal), both ends carried in every
cost figure. Starship-class $100–500/kg appears **only** as a labelled sensitivity
and in no headline number (AM-9). `SIM-9` `SOURCED, dated`

### 2.3.3 Grams of rhenium per gram of aluminium

**On target: 68–171 mg per g** (800 km, 619 m/s, β ∈ [1.0, 2.5]).
**Launched: ×16.7** at the design guidance budget → 1.14–2.85 g per g.
`SIM-1, SIM-3` `DERIVED`

### 2.4 USD wasted per missed shot

**~$31 of consumables** — the material and launch cost of one shot. Under the
delivery model a partial miss delivers partial momentum rather than nothing, so
the cost of imperfect aim appears mainly as the areal penalty already charged on
every shot (16.7×), plus a 1/0.9 shot multiplier for the 90% delivery confidence.
`SIM-9` `DERIVED`

---

## 3. Is it acceptable to deorbit this material?

### 3.1 Environmental effect

**Negligible.** 27 kg/yr of metal at 10,000 engagements/yr, ≈**7×10⁻⁷ of the
natural cosmic-dust influx** (~40,000 t/yr). Re₂O₇ melts at 297 °C and sublimes,
so ablated rhenium disperses as vapour rather than particulate. `SIM-8`
`DERIVED` + `SOURCED`

### 3.2 Is it legal? If not, how to make it so?

**Disposal law: PASS by construction.** Released material reenters within one
orbit, so FCC 47 CFR 25.283 (5 yr) and NASA ODMSP (25 yr) are satisfied at every
altitude. The debris being remediated is not the operator's object, so its
lifetime is an effectiveness benchmark, not a licence condition.

**The gap is authorisation, not disposal.** No regime permits a third party to
apply force to another state's registered object, however small. Deliberate
release of a hypervelocity cloud is a fault-liability exposure under the 1972
Liability Convention even though the material self-disposes. Mitigations:
conjunction screening of the release corridor before each shot (operationally
identical to standard collision avoidance), and an explicit international
small-debris remediation regime. `SIM-10` `SOURCED` + open

### 3.3 Can we track the trash so aircraft are not hit?

**Yes, and the hazard is nil regardless.** Release-to-reentry is ~45 minutes;
NOTAM window ~75 minutes; footprint ~2,500 × 100 km. The largest surviving
fragment is 0.45 mm carrying ~10⁻⁴ J — five orders of magnitude below the 15 J
injury threshold. `SIM-10` `DERIVED`

### 3.4 Will trash fully burn up and not hit the ground?

**Aluminium debris demises completely up to 5 cm; 10 cm is 98.9% ablated** with a
remnant carrying 5×10⁻⁴ J. Ground-casualty expectation **0**, against the NASA
1×10⁻⁴ threshold. `SIM-8` `DERIVED`

### 3.5 Will rhenium fully burn up and not hit the ground?

**Yes across the entire design grain range.** 0.05–1 mm rhenium grains ablate
>99.99% and fully demise. Larger pieces (2 mm+) partially survive but are never
released in this design. Casualty expectation **0**. `SIM-8` `DERIVED`

---

## Cross-cutting answer: is rhenium the right material?

**No — tungsten is.** Momentum p = mv is density-independent, so the founding
rationale for rhenium has no bearing on the momentum budget. Tungsten is 73–243×
cheaper, 1,037× more available, chemically better behaved, and has a 5.7% *higher*
area-to-mass ratio at equal grain mass, so missed tungsten decays slightly faster.
The single axis on which density could have favoured rhenium runs against it.
Because launch cost dominates delivered cost at these per-shot masses, the switch
changes total cost by ~1.6×, not ~100× — so the material choice does not decide
viability, but rhenium cannot be defended on physics. `SIM-7` `SOURCED`
