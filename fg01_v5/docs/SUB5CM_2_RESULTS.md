# FG01 / KIDR — Results
## Momentum-transfer remediation of sub-5 cm LEO debris

**Document 2 of 3** · [Methodology](SUB5CM_1_METHODOLOGY.md) · [Viability](SUB5CM_3_VIABILITY.md — not present in this repository)

All values regenerate from `python run_all.py`, seed `20260723`.
Tags: `DERIVED` · `SOURCED` · `UNVALIDATED`.

> **Superseded-economics notice.** Sections 8 (Joint viability) and 9 (Economics)
> below assume 10⁴ engagements per platform lifetime as a sampled/assumed
> parameter — the same convention used in the original SIM-12/13 uncertainty
> quantification. That convention is superseded elsewhere in this project's
> current manuscript: ECO-1 derives achievable engagement rate in closed form
> from orbital mechanics (Δv per engagement = π·v_orbital/(3.5·Ω̇·T_mission)) and
> finds it fixed at approximately 24–55 per platform-*decade* at the best debris
> cluster examined, not 10⁴ per lifetime. The P(viable) figures, the $27,469
> median cost/object, and the "cost per object is essentially flat across the
> sub-5 cm class" conclusion in §9 all depend on the superseded assumption and
> should not be cited as current without re-deriving them under the corrected
> cadence. Sections 1–7 (validation, disruption gate, decay gate, interaction,
> terminal guidance, net debris, reentry) do not depend on engagement count and
> are not affected by this notice. `SUB5CM_3_VIABILITY.md`, referenced above as
> the concluding document in this series, was not provided alongside Documents
> 1 and 2 and is not present in this repository.

---

## 1. Validation

Reported first, because everything downstream is conditional on it.

| # | Check | Result | Level |
|---|---|---|---|
| 1 | Δv_shift vs. first burn of a two-burn Hohmann to the same perigee | **0.0** (exact identity) | analytic |
| 2 | E_s substitution vs. direct ½m_p v²/m_t, 2,000 random draws | **4.9×10⁻¹⁶** | analytic |
| 3 | Density table vs. fresh direct NRLMSISE-00 calls | **1.9×10⁻⁴** max rel. | internal |
| 4 | Orbit-averaged vs. 3-D Cartesian RK4, 400 revolutions, J2 off | **0.01 %** | independent method |
| 5 | Cartesian step-size convergence (20 s → 10 s) | 0.157 % | numerical |
| 6 | Full Cartesian to reentry, 250 km, 1 cm | 1.9 % | independent method |
| 7 | **Starshine 1** (NORAD 25769, 387 km, i = 51.6°, 258 d observed) | **233 d nominal (−9.8 %)**, bracketed by 1011 d (solar min) and 94.6 d (solar max) | external observation |
| 8 | Closed-form vs. exact non-central-χ² delivery sizing | 0.1–0.2 % | analytic |
| 9 | Mass efficiency vs. grain-by-grain Monte Carlo | **2.16 MC standard errors** (statistically consistent) | independent method |
| 10 | β model at 6 km/s vs. published hypervelocity values | 1.96–8.27 against reported 2–3 | cross-regime |
| 11 | Closed forms vs. explicit 3-D inertial states | \|w\| **2.1×10⁻¹⁶**, cos ψ **2.8×10⁻¹¹** | analytic |
| 12 | Monte-Carlo convergence, 10³ → 10⁶ samples | 0.5760 → 0.5749 | numerical |
| 13 | Lifetime ∝ debris diameter | **exact**, as A/m ∝ 1/d requires | internal consistency |

The Starshine comparison is the only check against a real observation. The
observed lifetime falls between the nominal and solar-max cases, consistent with
the rising phase of solar cycle 23 over 1999–2000. Its weakest input is the
assumed 39.5 kg spacecraft mass, which is not independently confirmed and to
which the error scales.

---

## 2. Gate 1 — Catastrophic disruption: **PASS**

Design specific energy `E_s = γ·Δv(h)·v_rel / (2β)` is **independent of target
mass**, so a single result covers the whole sub-5 cm class. `DERIVED`

| h (km) | Δv (m/s) | E_s at β = 1.201 (kJ/kg) | Fraction of 40 kJ/kg | vs. 30 kJ/kg bound |
|---|---|---|---|---|
| 400 | 57.72 | 29.8 | 0.75 | PASS |
| 600 | 55.23 | 28.5 | 0.71 | PASS |
| 800 | 52.91 | 27.3 | 0.68 | PASS |
| **900** | **51.81** | **26.7** | **0.67** | **PASS** |
| 1200 | 48.72 | 25.1 | 0.63 | PASS |

Across the full β envelope (1.0–2.5) the design point spans **12.1–35.7 kJ/kg**,
i.e. 0.30–0.89× the threshold, and clears 40 kJ/kg at every altitude.

**Self-disposal floor.** Requiring v_rel ≥ Δv_direct(h) puts all missed mass on
an immediate-reentry trajectory. The floor is 87–290 m/s across 400–1200 km,
against a design v_rel of 619 m/s — margin **2.1–7.1×**. Satisfied everywhere.
`DERIVED`

**Per-impact energy (AM-16).** At 450 µm grains and 619 m/s the per-grain
specific energy on a 1 cm target is **136 J/kg = 0.34 % of threshold**. The bulk
figure alone already clears the gate, so the swarm route is a **reserve, not a
dependency**.

---

## 3. Gate 2 — Orbital decay: **PASS below a size-dependent ceiling**

Natural and post-kick lifetimes, nominal solar activity, Δh = 200 km, years.
`SOURCED` (NRLMSISE-00)

| h (km) | 0.5 cm nat. → post | 1 cm nat. → post | 2 cm nat. → post | 5 cm nat. → post |
|---|---|---|---|---|
| 400 | 0.03 → 0.00 | 0.06 → 0.00 | 0.12 → 0.01 | 0.31 → 0.02 |
| 600 | 0.80 → 0.12 | 1.59 → 0.25 | 3.19 → 0.50 | 7.96 → 1.24 |
| 700 | 3.19 → 0.60 | 6.37 → 1.21 | 12.75 → 2.42 | 31.87 → **6.04** |
| 800 | 11.0 → 2.52 | 22.0 → 5.05 | 44.0 → **10.10** | 110 → 25.24 |
| **900** | 31.8 → **9.09** | 63.5 → **18.18** | 127 → 36.37 | ∞ → 90.92 |
| 1000 | 75.6 → 27.49 | 151 → 54.97 | ∞ → 110 | ∞ → ∞ |
| 1200 | ∞ → 142 | ∞ → ∞ | ∞ → ∞ | ∞ → ∞ |

**Lifetime scales exactly linearly with diameter** (A/m ∝ 1/d) — an internal
consistency check that holds to machine precision.

**Acceleration factor** falls from **16.7× at 400 km to 2.2× at 1100 km**: the
kick is least effective exactly where debris lifetimes are longest.

### Altitude ceilings (nominal solar, Δh = 200 km)

| Debris | 25-yr rule | 5-yr rule |
|---|---|---|
| 0.5 cm | **900 km** | 800 km |
| 1 cm | **900 km** | 700 km |
| 2 cm | **800 km** | 700 km |
| 5 cm | **700 km** | 600 km |

Staged engagement lifts the ceiling — each kick is independently
sub-catastrophic, so two successive 200 km kicks reach 1100–1200 km. The ceiling
is therefore **economic, not physical**.

---

## 4. Interaction and delivery

### 4.1 Momentum coupling β — a correction that worsens the result

| v_imp (m/s) | β low | β central | β oblique-corrected | Regime |
|---|---|---|---|---|
| 400 | 1.055 | 1.094 | 1.063 | transitional |
| **619** | **1.170** | **1.301** | **1.201** | transitional |
| 1000 | 1.307 | 1.568 | 1.378 | hypervelocity |
| 6000 | 1.957 | 3.072 | 2.381 | hypervelocity |

**β = 1.201 at the design point, not the 1.5 previously assumed.** At 619 m/s the
impact is only **1.9× aluminium's 320 m/s strength velocity**, so almost no mass
is ejected fast enough to enhance momentum — the interaction is close to pure
embedding. Required mass on target rises **1.25×** and E_s rises from 21.4 to
26.7 kJ/kg. `DERIVED, not SOURCED`

Supporting evidence: 450 µm rhenium grains penetrate **1.46 mm** into aluminium
at 619 m/s — 3.2× their own diameter, far short of perforating — confirming they
embed.

### 4.2 Cloud efficiency — a correction that improves it

A Gaussian cloud is **exactly 3× more mass-efficient** than a uniform disc:
ε = 2.995·(d_t/d_cloud)², because a peaked profile concentrates mass where the
target is. `DERIVED`

| Cloud diameter | ε on 0.5 cm | 1 cm | 2 cm | 5 cm |
|---|---|---|---|---|
| 2 cm | 0.171 | 0.527 | 0.950 | 1.000 |
| 5 cm | 0.030 | 0.113 | 0.381 | 0.950 |
| 10 cm | 0.007 | 0.030 | 0.113 | 0.527 |

The two corrections in §4.1 and §4.2 **very nearly cancel**.

### 4.3 Worked design point — 1 cm at 900 km

| Quantity | Value |
|---|---|
| Mass on target | **197 mg** |
| Launched mass (mass-optimal cloud) | 4.28 g |
| Launched mass at 95 % per-shot reliability | **2.14 g** |
| Areal penalty | 21.7× |
| Grains on target (450 µm Re) | **197** |
| Poisson CV of delivered momentum | 7.9 % |
| Bulk E_s | 26.7 kJ/kg (0.67× threshold) |
| Per-grain E_s | 136 J/kg (0.0034× threshold) |
| Consumable | $17.95 rhenium / **$5.41 tungsten** |

**Δv distribution across shots** (aim scatter + Poisson lumpiness, γ = 2):

| Cloud | Launched | Δv p05 | p50 | p95 | Meets requirement |
|---|---|---|---|---|---|
| 2 cm | 0.37 g | 1.1 | 41.9 | 105 | 42 % |
| **5 cm** | **1.75 g** | **43.7** | **88.0** | **118** | **91 %** |
| 10 cm | 6.68 g | 81.0 | 103 | 124 | 99.98 % |

Dominated by aim scatter, not Poisson lumpiness — at 200+ grains the swarm
averaging is effectively deterministic. This is the quantitative form of the
"many grains striking the debris" picture.

### 4.4 Size scaling within the sub-5 cm class

| Question | Answer |
|---|---|
| Do bigger targets catch more grains? | **Yes — dominant favourable scaling, but it saturates at 5.4 cm** for a 5 cm cloud |
| Do they receive more Δv? | **No — the opposite.** Removed-mass-per-launched-mass caps at **7.17 kg/kg** = β·v_rel/(γ·Δv), independent of target size |
| Do they decay faster? | **No.** A/m ∝ 1/d, so a 5 cm fragment takes 10× as long as a 0.5 cm fragment from the same orbit |

**Sub-5 cm is precisely the regime where interception has *not* saturated**: a
1 cm target catches only 11 % of a 5 cm cloud and removes 0.8 kg per kg launched,
against the 7.17 cap. The three scalings act in opposition, which is why the
viable envelope is a narrow diagonal rather than a region (§8).

---

## 5. Terminal guidance

The legacy formulation is reproduced exactly before being replaced.

| τ (ms) | M1 miss (legacy) | M1 P_hit | M3 miss (design) |
|---|---|---|---|
| 0.1 | 0.062 m | 0.987 | 5.775 mm |
| 1 | 0.619 m | 0.265 | 5.775 mm |
| 10 | 6.19 m | **2.4×10⁻⁵⁸** | **5.775 mm** |
| 100 | 61.9 m | 0 | **5.775 mm** |

**The M3 total is independent of control latency from 0.1 µs to 100 ms.**
Keplerian relative motion at 10 m separation is trivially predictable — the
unmodelled relative-acceleration term is **11 nanometres** — so latency enters
only through estimation error.

Error budget at the design point (R = 10 m, u = 10 m/s, 200 µrad, 100 MHz ranging):

| Term | Contribution |
|---|---|
| Range-knowledge timing | 5.41 mm |
| Launcher pointing | 2.00 mm |
| Muzzle-speed dispersion | 0.162 mm |
| Optical position knowledge | 0.055 mm |
| Unmodelled relative acceleration | 1.1×10⁻⁸ m |
| Velocity-extrapolation error | 5.0×10⁻⁹ m |
| **Total** | **5.775 mm** |

This is a consequence of the engagement geometry the energy gate already forces
(§ methodology 5.1), not an assumption of better hardware. Its price is that the
platform must operate inside the target's orbital plane.

**A concrete requirement:** ranging bandwidth **≥ 100 MHz**. At 1 MHz, σ_R =
33.5 m contributes 0.54 m of miss — inadequate.

---

## 6. Net debris: **PASS, +1 per engagement**

| Threshold | Item | Count |
|---|---|---|
| >1 cm | Target removed | **+1** |
| >1 cm | Ejecta (conservative bound: fragment = full crater) | **0, identically** |
| >1 cm | Embedded grains, missed grains | 0 |
| **>1 cm** | **NET** | **+1** |
| >1 mm | Ejecta, conservative / expected | −17.4 / 0 |

The 197 craters of a design engagement excavate 370 mg (26 % of a 1 cm target's
mass), but the **largest single crater is 1.39 mm** — 7× below the trackable
threshold. Production terminates at the crater diameter, so the engagement
**cannot** create a trackable fragment by this route. This is stronger than an
argument from the sub-catastrophic regime; it follows from crater geometry.

Ejecta are sub-millimetre with 10–1000× the parent's A/m: 10 µm ejecta decay in
0.02 yr against 17.6 yr for the parent.

The mandated SBM non-catastrophic branch returns 0.26 fragments > 1 cm from an
undisrupted 1.41 g target — physically impossible, and reported as the
conservative model artefact it is.

---

## 7. Reentry and environment: **PASS**

| Particle | Ablated | Outcome |
|---|---|---|
| 0.05–1 mm Re grains (**design range**) | > 99.99 % | **fully demise** |
| 2 mm Re | 83 % | remnant, 2×10⁻³ J |
| Al debris ≤ 5 cm | > 99.99 % | **fully demise** |

**Ground-casualty expectation: 0** at every population density tested, against
the NASA 1×10⁻⁴ threshold. Deposition at 10⁴ engagements/yr is 27 kg/yr,
≈7×10⁻⁷ of the natural cosmic-dust influx.

Released mass reenters within **one orbital period** (~45 min) at any
v_rel ≥ Δv_direct, so FCC 5-yr and ODMSP 25-yr disposal rules are satisfied
**by construction** at every altitude.

---

## 8. Joint viability

### 8.1 Global UQ — design variables sampled

P(all constraints satisfied), nominal solar, 10⁶ samples. Two cost conventions
and the remediation-value gate, reported together (AM-18).

**With remediation-value gate — the operative result:**

| h (km) | 0.5 cm | 1 cm | 2 cm | 5 cm |
|---|---|---|---|---|
| 400–600 | 0 | 0 | 0 | 0 |
| 700 | 0 | 0 | 0 | **0.400** |
| 800 | 0 | 0 | **0.508** | 0.168 |
| **900** | **0.582** | **0.511** | 0.021 | 0 |
| 1000 | 0.193 | 0 | 0 | 0 |
| 1100–1200 | 0 | 0 | 0 | 0 |

**Without the remediation gate** (per-object cost), values at 400–800 km rise to
0.36–0.57 across all sizes — those configurations pass every physical gate while
engaging debris that would have decayed unaided.

**Under the per-kilogram criterion** ($100,000/kg), P(viable) is **identically
zero for all sub-2 cm debris at every altitude**. This is a property of the
metric, not the design: a 1 cm sphere masses 1.414 g, so the bar is $141/object.

### 8.2 Sensitivity

| Parameter | S₁ | S_T |
|---|---|---|
| Engagements per platform life | 0.223 | **0.519** |
| Altitude | 0.240 | **0.497** |
| Platform cost | 0.022 | 0.177 |
| β | 0.017 | 0.169 |
| v_rel | 0.012 | 0.153 |
| Δh | −0.005 | 0.142 |
| **Control latency τ** | — | **absent from the top eight** |

The dominant terms are **programmatic** (engagement count, platform cost) and
**environmental** (altitude), not physical.

### 8.3 Selected operating point — design variables chosen

Epistemic-only Monte Carlo with design variables fixed (R = 10 m, 200 µrad,
u = 10 m/s, τ = 10 ms).

**Optimum: 1 cm debris at 900 km, Δh = 200 km, v_rel = 619 m/s → P(viable) = 0.662**
(63.5 yr natural → 18.2 yr post-kick).

| Constraint | P(satisfied) |
|---|---|
| **Cost < $50 k/object** | **0.663 ← binds** |
| Sub-catastrophic E_s | 0.998 |
| Net debris ≥ 0 | 0.998 |
| Shot mass ≤ 1 kg | 1.000 |
| Post-kick lifetime < 25 yr | 1.000 |
| Remediation value | 1.000 |

The two numbers 0.662 (chosen design) and 0.511 (sampled design, §8.1) differ
because they answer different questions; both are reported.

**Conditional on the dominant term:**

| Engagements per platform lifetime | P(viable) | Median cost/object |
|---|---|---|
| 1,000 | 0.000 | $274,446 |
| 3,000 | 0.222 | $91,497 |
| **10,000** | **0.997** | **$27,469** |
| 100,000 | 0.998 | $2,772 |

### 8.4 Robustness

P(viable) is flat to within 0.01 across control latency (0.1 µs – 100 ms),
engagement range (5–100 m), pointing (50–1000 µrad) and muzzle dispersion
(0.1–1 %). Only the platform–target relative speed shows any effect, and weakly
(0.661 → 0.651 from 1 to 100 m/s).

---

## 9. Economics

Cost per object at 900 km, design guidance budget, 10⁴ engagements, $50 M platform:

| Debris | Shot mass | Kicks | **Cost/object** | Cost/kg |
|---|---|---|---|---|
| 1 cm | 2.66 g | 1 | **$5,030** | $3.56 M |
| 2 cm | 5.31 g | 1 | **$5,061** | $447 k |
| 5 cm | 13.29 g | 1 | **$5,152** | $29.2 k |

**Cost per object is essentially flat across the sub-5 cm class** (2 % from 1 cm
to 5 cm) because it is **99.4 % platform amortisation** — consumables are ~$31.
Cost per kilogram falls 122× over the same range purely because the mass removed
rises.

Consequently: cost/object ≈ C_platform / N_engagements, and nothing about the
projectile, material, launcher or grain size materially affects it.

---

## 10. Where the gates fail

Reported in the same document as the successes (AM-18).

- **Above the size-dependent ceiling** (900 km for ≤1 cm, 800 km for 2 cm,
  700 km for 5 cm) a single kick does not achieve 25-yr compliance.
- **Below ~700 km** the kick has **no remediation value** — the debris decays
  unaided within 25 years.
- **Sub-2 cm debris fails the $100,000/kg criterion at every altitude**, by 1–2
  orders of magnitude.
- **Below ~3,000 engagements per platform lifetime**, P(viable) = 0 everywhere.
- **At β = 1.0 with E_s,c = 30 kJ/kg**, the margin narrows substantially though
  the window never closes.
- **No legal regime authorises the operation** on another state's registered
  object, however small.

---

## 11. Model corrections adopted

Three corrections were found during validation and adopted rather than tuned
around. Two worsen the result.

| Correction | Direction | Effect |
|---|---|---|
| β derived (1.201) rather than assumed (1.5) | **worse** | +25 % mass on target, E_s 21.4 → 26.7 kJ/kg |
| Gaussian rather than uniform-disc cloud | better | ε 3× higher |
| Remediation-value gate added | **worse** | Eliminates all configurations below ~700 km |

An earlier delivery-optimisation objective (maximise ΔKE) was also found to be
the right answer to the wrong question and replaced with stage/shot efficiency.

---

*Document 2 of 3. Conclusions and the viability assessment follow in
[SUB5CM_3_VIABILITY.md](SUB5CM_3_VIABILITY.md).*
