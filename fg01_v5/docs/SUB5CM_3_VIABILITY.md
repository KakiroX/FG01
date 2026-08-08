# FG01 / KIDR — Viability Assessment and Conclusions
## Momentum-transfer remediation of sub-5 cm LEO debris

**Document 3 of 3** · [Methodology](SUB5CM_1_METHODOLOGY.md) · [Results](SUB5CM_2_RESULTS.md)

---

> **Superseded-economics notice.** The P(viable) figures in Sections 1 and 2 —
> the headline P(viable) = 0.662 / 0.582, the per-gate probabilities such as the
> cost gate "binds at P = 0.663", and the altitude×size probability table in §2 —
> are outputs of the original SIM-12/13 joint uncertainty quantification, which
> samples *engagements per platform lifetime* as a free parameter (the same
> sampled/assumed convention, spanning up to the 10⁴-per-lifetime case, used in
> the companion Documents 1 and 2). That convention is superseded by ECO-1, which
> derives achievable engagement rate in closed form from orbital mechanics
> (Δv per engagement = π·v_orbital/(3.5·Ω̇·T_mission)) and fixes it at
> approximately 24–55 per platform-*decade* — not a sampled quantity — at the best
> debris cluster examined. Two things follow. First, the specific P(viable)
> probabilities above should not be cited as current without re-deriving them
> under the fixed cadence, and the Sobol total-effect term "engagements per
> platform lifetime (0.519)" in §3 is a sensitivity to that superseded sampled
> parameter. Second — and unlike the companion documents — the *verdict itself
> does not rest on the old numbers*: §3 already presents the ECO-1
> Δv-per-engagement invariant as the binding constraint and derives the same
> "73× short" shortfall from it directly, so the physics findings (§5), the
> viable-envelope geometry (§2), the value comparison (§4) and the recommendation
> not to deploy (§8) stand under the corrected cadence. It is the sampled
> P(viable) values, not the conclusion, that are superseded.

---

## 1. Verdict

> **The physics works. The economics do not.**
>
> For debris ≤ 5 cm, every physical gate passes at the selected operating point
> with probability ≈ 1. The concept fails on cost per object, by a factor of
> **73× at best** against its own stated target — and the binding constraint is
> neither the kick physics, the debris population, nor the platform cost, but a
> closed-form orbital-mechanical invariant governing how many objects one
> platform can reach.

| Gate | Verdict | Margin |
|---|---|---|
| 1. Sub-catastrophic specific energy | **PASS** | 0.63–0.75× threshold; P = 0.998 |
| 2. Launched mass within magazine | **PASS** | P = 1.000 |
| 3. Post-kick lifetime < 25 yr | **PASS** below a size-dependent ceiling | P = 1.000 within envelope |
| 4. Net debris ≥ 0 | **PASS** | +1 per engagement; P = 0.998 |
| 5. Cost per object | **FAIL** | binds at P = 0.663; 73× short of target |
| 6. Remediation value | **PASS** only in a narrow diagonal band | see §2 |

**P(viable) = 0.662** at the selected optimum (1 cm at 900 km), and **0.582** as
the maximum over the sampled-design global UQ (0.5 cm at 900 km). Neither
approaches the ~100 % that would justify deployment.

---

## 2. The viable envelope

Viability for sub-5 cm debris is **not a region but a narrow diagonal band**,
because three size-scalings act in opposition:

| Scaling | Direction | Effect of increasing target size |
|---|---|---|
| Interception efficiency ε | **favourable** | Bigger target catches more of the cloud (saturates at 5.4 cm) |
| Removed mass per launched mass | neutral | Caps at 7.17 kg/kg — **independent of target size** |
| Area-to-mass ratio | **unfavourable** | A/m ∝ 1/d, so bigger debris decays proportionally slower |

The net result, from the joint UQ with the remediation-value gate applied:

| Altitude | 0.5 cm | 1 cm | 2 cm | 5 cm |
|---|---|---|---|---|
| ≤ 600 km | 0 | 0 | 0 | 0 |
| 700 km | 0 | 0 | 0 | **0.400** |
| 800 km | 0 | 0 | **0.508** | 0.168 |
| **900 km** | **0.582** | **0.511** | 0.021 | 0 |
| 1000 km | 0.193 | 0 | 0 | 0 |
| ≥ 1100 km | 0 | 0 | 0 | 0 |

**Each debris size has essentially one workable altitude**, and it moves *down*
as size increases. The envelope is bounded on all four sides:

- **Below ~700 km** — the debris decays unaided inside 25 years; the engagement
  achieves nothing (remediation-value gate).
- **Above the ceiling** (900 km for ≤ 1 cm, 800 km for 2 cm, 700 km for 5 cm) —
  a single kick cannot bring the object inside the 25-yr rule.
- **Below 2 cm under a per-kilogram cost criterion** — P(viable) is identically
  zero at every altitude, because a 1 cm sphere masses 1.414 g and a
  $100,000/kg bar is $141 per object.
- **Above 5 cm** — outside this study's scope; interception saturates and the
  economics change character.

This narrowness is itself a finding. **The sub-5 cm class is the regime in which
interception has not yet saturated**: a 1 cm target catches only 11 % of a 5 cm
cloud and removes 0.8 kg per kg launched, against a 7.17 kg/kg cap.

---

## 3. Why it fails: the Δv-per-engagement invariant

Cost per object is **99.4 % platform amortisation** — consumables are ~$31, and
cost/object ≈ C_platform / N_engagements. So viability reduces to a single
question: *how many objects can one platform service?*

That number is governed by a closed form, derived and validated to 1.14 %:

$$\boxed{\;\Delta v_{\text{per engagement}}=\frac{\pi\,v_{\text{orbital}}}{3.5\,\dot\Omega\,T_{\text{mission}}}\;}$$

**This cannot be tuned.** Widening the altitude band raises engagement
opportunities linearly (through differential nodal drift) and raises the hop Δv
linearly (through the Hohmann transfer) — the two cancel exactly. Only orbital
velocity, nodal regression rate and mission duration enter, none of which the
designer controls once the target population is chosen.

At Fengyun-1C over ten years this is **107 m/s per engagement**. Consequences:

| Δv budget | Engagements per decade | Cost/object at $20 M platform |
|---|---|---|
| 500 m/s | 16 | $1.25 M |
| 2,000 m/s | 24 | $833 k |
| 10,000 m/s | 55 | **$364 k** |

Against the programme's own **$5,000/object** niche target, the best case is
**73× short**. Reaching 10,000 engagements would require ~1,070 km/s of tour Δv,
which no propulsion system can deliver.

**Two further findings compound this.** Every breakup cloud studied has spread to
a full 360° in RAAN within a few years, so the compact co-orbital cluster the
architecture assumed does not exist — what survives is a shared *inclination
band*. And because θ² ≈ (Δi)² + (sin i·ΔΩ)², a fragment whose inclination differs
from the platform's by more than the firing cone can **never** be engaged, at any
RAAN, at any time, for any Δv expenditure.

Sensitivity analysis confirms the diagnosis: the dominant Sobol total-effect
terms are **engagements per platform lifetime (0.519)** and **altitude (0.497)**,
both programmatic or environmental. Control latency — the variable that dominated
an earlier version of this analysis — does not appear in the top eight.

---

## 4. Value comparison

Removing one randomly-selected ≥ 1 cm fragment avoids **$67–$20,124** of expected
asset loss across the swept assumptions, and **at most ~$67** in the subset
consistent with the observed on-orbit loss record.

| | Cost per object | Avoided loss | Gap |
|---|---|---|---|
| Most generous assumptions | $364 k | $20,124 | **18×** |
| Credible assumptions | $364 k | $67 | **~5,400×** |

This is not a failing peculiar to FG01. It is why operational active debris
removal targets **large** objects, whose removal prevents the breakups that
create fragments in the first place, rather than the fragments themselves. **The
case for sub-5 cm removal cannot be made on per-object value; it requires a
systemic cascade-prevention argument that this study does not model.**

### Against ground-based laser ablation

The only competing method that genuinely addresses this size class.

| | FG01 | Ground laser |
|---|---|---|
| Cost per 1 cm object | $364 k – 4.2 M | **$381 – 5,702** |
| Prime hardware vs. demonstrated | **80,491× inside** | **46× beyond** |
| Pointing requirement | 876 µrad | 0.21 µrad |
| Margin over own diffraction limit | **131×** | 1.64× |
| Aiming solutions per object | **1** | 294 consecutive |
| Failure mode on aiming error | **graceful** (45 % delivered at 2× error) | threshold cliff |

The laser retains a **159–10,951×** cost advantage. Two honest qualifications:
that advantage is contingent on hardware **46× beyond demonstrated average
power**, and the published $100–500/kg figure it rests on is a *marginal-energy*
cost that omits the facility — correcting it moved the gap from six orders of
magnitude to two-to-four.

On terminal-accuracy tolerance specifically the comparison runs the other way,
because FG01's short-range areal-capture engagement degrades smoothly where a
long-standoff laser's point-hit requirement is binary. **That is a geometry
consequence, not a cost offset**, and it does not change the verdict.

---

## 5. What is genuinely established

Findings that hold independently of the economic verdict, and that constitute the
study's positive contribution:

1. **The disruption falsification is lifted.** Design specific energy is
   **12.1–35.7 kJ/kg** across the full β envelope, below the 40 kJ/kg threshold
   at every altitude — and it is **mass-independent**, so this covers the entire
   sub-5 cm class in one result.

2. **Missed mass is not debris.** Fired retrograde above Δv_direct(h), missed
   projectile mass has sub-atmospheric perigee and reenters within one orbital
   period. FCC and ODMSP disposal rules are satisfied **by construction** at
   every altitude. This eliminates what had been considered a primary risk.

3. **Net debris is +1 per engagement, by crater geometry.** The largest crater
   from a design engagement is **1.39 mm** — 7× below the trackable threshold —
   so the engagement *cannot* produce a trackable fragment. This is stronger than
   an argument from the sub-catastrophic regime.

4. **The guidance wall was a sweep-rate artefact.** The corrected co-orbital
   budget gives a **5.775 mm** miss that is independent of control latency from
   0.1 µs to 100 ms. The co-orbital geometry this relies on is *forced* by the
   energy gate (a crossing engagement exceeds the catastrophic threshold at every
   plane angle), not assumed for convenience.

5. **β is now derived rather than assumed**, moving from a factor-2.5 envelope to
   1.170–1.681, with the model validated against the hypervelocity regime where
   published data exists.

6. **Complete demise on reentry.** Design-range grains ablate > 99.99 %;
   ground-casualty expectation is **0** against the NASA 1×10⁻⁴ threshold.

---

## 6. Limitations bearing on this verdict

Ordered by their capacity to change the conclusion.

| # | Limitation | Could it overturn the verdict? |
|---|---|---|
| **L1** | Engagement rate rests on modelled cluster populations, not a catalogue (only ~3 % of the ≥1 cm population is tracked) | **No** — the Δv invariant is orbital-mechanical and independent of how many fragments exist |
| **L2** | Target acquisition and terminal tracking assumed, not designed | **No** — would worsen it |
| **L9** | Platform cost parametric ($20–500 M) | **Partially** — the $364 k best case already assumes the $20 M floor |
| **L5** | β derived, not measured | **No** — enters mass, and consumables are 0.6 % of cost |
| **L10** | Fixed F10.7 per case, not an 11-year cycle | **Partially** — no operating point is simultaneously remediation-relevant at permanent solar max and compliant at permanent solar min; a cycle-averaged propagation would sharpen the envelope |
| **L14** | Solid-sphere debris (AM-2) | **Narrows** the envelope further — real fragments have higher A/m and shorter natural lifetimes |
| **L17** | Interim on-orbit hazard argued from the lifetime ratio, not integrated | **No** — direction robust |
| — | Value model spans three orders of magnitude, and its upper end implies ~80 debris losses/year against an observed record of essentially none | **No** — the gap is ≥ 18× everywhere in the sweep |

**No single limitation, if resolved favourably, moves cost per object from
$364 k to $5,000.** The shortfall is structural.

---

## 7. What would change the answer

In descending order of impact:

1. **A fundamentally different visiting mode.** The Δv invariant is the entire
   result. Anything that removes the need to match each target's altitude
   individually would reopen the economics. INT-0 shows this cannot be done with
   a momentum kick — a retrograde impulse requires travelling anti-parallel to
   the target's velocity — so it would require a different removal mechanism, not
   a different FG01.

2. **A systemic cascade-value model.** Per-object value is what condemns sub-5 cm
   removal. Cascade prevention is the only framing under which it could pay, and
   it is not modelled here.

3. **Retargeting at larger objects.** Outside this document's scope, but the same
   launcher and the same sub-catastrophic energy apply to any target, and the
   areal penalty collapses to 1 above ~5 cm.

4. **A direct β measurement** at Al-on-Al, cm-scale, 0.5–1 km/s — the widest
   remaining physics uncertainty, though it does not affect the verdict.

---

## 8. Recommendations

**Do not deploy FG01 against sub-5 cm debris.** The concept is sound physics
attached to an unviable mission. Specifically:

- **Do not** pursue the small-debris tour architecture. It fails by 73× at best
  and the binding constraint cannot be engineered away.
- **Do** preserve the kick physics as a validated result. Items 1–6 of §5 are
  robust, reproducible and independently useful.
- **Do** report the negative result. A 0 % viability claim that a reviewer can
  dismantle is worth less than an honestly bounded failure with a named cause,
  and the Δv invariant is a general constraint on *any* co-orbital
  momentum-transfer remediation scheme, not a defect specific to this design.
- **Consider** re-scoping to intermediate-mass targets, where the areal penalty
  vanishes while the visiting cost is unchanged.
- **If sub-5 cm removal is pursued regardless**, the only defensible operating
  points are those in the §2 diagonal — and each must be justified against a
  cascade-prevention value model that does not yet exist.

---

## 9. Statement of compliance

| Requirement | Status |
|---|---|
| Traceability tag on every number | Applied (`DERIVED`/`SOURCED`/`UNVALIDATED`) |
| Anti-mistake contract AM-1 … AM-18 | Discharged per-module; mapping in Methodology §2.2 |
| Ranges rather than point estimates (AM-10) | Applied throughout |
| Failure envelope reported with successes (AM-18) | Results §10, Viability §2 |
| Validation at analytic, independent-method, external-observation and cross-regime levels | 13 checks, Results §1 |
| Limitations consolidated, not scattered | Methodology §9, Viability §6 |
| Model corrections adopted rather than tuned around | Three recorded, two adverse (Results §11) |
| Reproducibility: fixed seed, pinned dependencies, single entry point | Methodology §10 |
| Provenance separating verified from inherited sources | `CITATION_LOG.md` |

**Declared adverse findings.** Three corrections found during this work made the
result worse and were adopted: the derived β (+25 % mass), the remediation-value
gate (eliminates all configurations below ~700 km), and the areal-capture
delivery model replacing a hit-probability formulation. A fourth — that the
optimiser's original objective maximised the wrong quantity — was corrected
before it influenced any reported result.

---

*End of document 3 of 3.*
