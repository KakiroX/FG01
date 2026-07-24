# FG01 / KIDR v6 — Cloud–Debris Interaction Model & Economically-Viable Engagement Architecture

**Builds on:** `FG01_Agent_Context_and_Simulation_Plan_v5.md` (context, anti-mistake contract, parameter table) and the completed v5 simulation run.
**Two jobs of this document:** (I) replace the single-impactor stand-in with a real, resolved **cloud–debris interaction model**; (II) **fine-tune the mission architecture for economic viability**, since the completed run showed cost — not physics — is the only binding constraint.
**Audience:** simulation engineers and reviewers continuing the FG01 work.

> **Framing note (read this).** The completed v5 run found all four *physics* constraints pass with probability ≈1 at the optimum; the sole blocker is cost, which is 99.4% platform amortization. Economic viability is therefore **not** a matter of tuning the grain physics — the rhenium is ~$22–88 per engagement, a rounding error. It is entirely a matter of **how many objects one platform can service**. This plan is built around that fact: Part I gets the interaction physics right so the per-engagement result is trustworthy; Part II reorganizes the *mission* so the platform is amortized over enough targets to be economical. Where a viability figure depends on an assumption the simulations must still confirm, it is tagged `UNVALIDATED` — a fine-tuned architecture is only as good as the target population it assumes.

---

## PART 0 — Findings from the completed v5 run, carried in as premises

These are inherited results this plan builds on. Each is retained with its validation status; the ones tagged **[verify]** are load-bearing and must survive independent scrutiny before any economic claim rests on them.

1. **The disruption gate is cleared.** Design-point specific energy 12.1–35.7 kJ/kg across the β envelope, below 40 kJ/kg at every altitude (Hohmann cross-check exact; E_s substitution to 5×10⁻¹⁶). *Solid.*
2. **Missed grains are not debris.** Fired retrograde at v_rel > Δv_direct(h), the missed cloud's perigee is sub-atmospheric; it reenters on the first perigee pass (~45 min). Defuses the missed-mass persistence risk and makes FCC/ODMSP compliance automatic. *Geometrically sound; keep the v_rel > Δv_direct(h) condition explicit as a firing constraint.*
3. **The engagement is near-co-orbital, low-relative-velocity. [verify]** The energy gate forces v_rel low; a co-orbital target moves across the aim point at u ≈ 10 m/s, not at the 619 m/s closing speed (which comes from the launcher). This is what turns v4's latency "cliff" into a latency-independent ~5.8 mm miss. **This single reframing reverses v4's dominant finding and is the most attackable claim in the whole project.** It is a premise of everything below; it must be independently validated (Part I, INT-0).
4. **Hit/miss is the wrong variable.** A target intercepts only the mass incident on its cross-section, so aim error costs *launched mass* (∝ d_miss²/A_target), a graceful quadratic — exactly the multi-grain picture the PI described. Formalized in Part I.
5. **Cost is 99.4% platform amortization** ($5,031/object; consumables ~$31). The dominant Sobol terms are engagements-per-platform (0.52) and altitude (0.50); control latency does not appear. *This is the entire economic problem.*
6. **Rhenium is not materially justified;** tungsten wins on every differentiating axis (including the missed-grain-decay axis, by 5.7%). Carried into Part II as a material decision, not re-litigated.
7. **P(viable) = 0.67 at the optimum** (1 cm target, 900 km, Δh 200 km, 63.5 yr → 18.2 yr), rising to 0.997 only once one platform services ≥10,000 objects — a number the run could not derive because no acquisition-cadence model existed. **Part II builds that model.**

---

# PART I — The cloud–debris interaction model

**Purpose.** Replace the single-macroscopic-impactor momentum balance (v5 SIM-1) and the point-hit probability (v5 SIM-4) with one resolved model of a dispersed rhenium cloud coupling momentum to a target it partially overlaps. This closes AM-15/AM-16 by modeling the swarm explicitly instead of assuming equivalence.

The interaction runs INT-0 → INT-1 → INT-2 → INT-3 → INT-4 → INT-5 in sequence; INT-3 exports the per-engagement Δv distribution that Part II's economics consume.

---

### INT-0 — Engagement geometry, frames & the co-orbital premise (validate the load-bearing claim)

**Purpose.** Pin down the geometry every later module assumes, and **independently test premise #3** — that the target's transverse sweep across the aim point is ~u (co-orbital drift), not v_rel.

**Questions answered.** From what relative state does FG01 fire? What is the true transverse velocity of the target relative to the intended intercept point, resolved into along-track / cross-track / radial components? Is the 5.8 mm, latency-independent miss real or an artefact of collapsing 3-D relative motion onto one axis?

**Governing equations.** Clohessy–Wiltshire relative state between FG01 and target in the Hill frame; decompose the closing velocity into the line-of-sight (closing) component and the transverse component. Miss during control latency τ: `d_miss,transverse = v_transverse · τ`, with `v_transverse` the CW cross-line-of-sight rate — **not** the scalar closing speed. `DERIVED`

**Inputs.** FG01–target relative state: separation R = 1–100 m; relative velocity magnitude = v_rel (500–1134 m/s from the SIM-0 safe zone); geometry swept over the full range of approach angles, *including the worst-case (near-transverse) approach*, not just the head-on retrograde best case. τ = 0.1–100 ms.

**Outputs.**
1. Table: v_transverse and d_miss vs. approach geometry and τ — showing the *range*, so the best case (near-zero transverse) and worst case (near-transverse pass) are both visible.
2. Explicit statement: "The latency-independent miss holds only for approach geometries within X° of pure retrograde; outside that cone, transverse rate grows to Y m/s and the latency term reappears." This is the honest boundary of premise #3.
3. The firing-geometry constraint that Part II's cadence model must respect (you may only engage a target when the approach is within the validated cone).

**Anti-mistake checks.** AM-5, AM-14, AM-18. **Validation.** Reproduce with a full 6-DOF relative-motion propagator (GMAT or equivalent); confirm the 1-D and 3-D miss agree only inside the stated cone. **If they do not, premise #3 fails and the guidance problem returns — report that outcome plainly.**

---

### INT-1 — Cloud generation & areal density field

**Purpose.** Describe the dispersed powder as a resolved areal number-density field at the intercept plane, `n(x, y)` [grains/m²], from launch conditions — the quantitative version of "many grains touching the debris."

**Questions answered.** Dust vs. grain mode; grain size; how the cloud's spatial profile evolves to the intercept plane; what fraction of launched mass lands within the target silhouette.

**Governing equations.** Launch dispersion → CW propagation (v5 SIM-3 Eq. 7) → 2-D areal density at time-of-flight R/v_rel. Model the transverse profile as (near-)Gaussian with σ_cloud(t) from CW; areal mass density `μ_A(x,y) = m_grain · n(x,y)`. `DERIVED` (grain-swarm formation from a solid slug remains `UNVALIDATED` — flag per AM-15).

**Inputs.** Total launched mass m_L = 0.5 mg–1 g; grain diameter d_g = 50 µm–2 mm (ρ_Re → grain mass); dispersion half-angle / δv_launch giving cloud diameters d_cloud = 2, 5, 10, 20 cm at intercept; N_grains = m_L/m_grain.

**Outputs.**
1. Table: d_cloud, N_grains, peak areal density, and **mass efficiency ε = (mass within target silhouette)/(mass launched)**. For a 1 cm target: d_cloud = 2/5/10/20 cm ⇒ area ratio ≈ 4/25/100/400× ⇒ ε ≈ 0.25/0.04/0.01/0.0025.
2. Plot: ε vs. d_cloud, with grains-per-target-silhouette annotated (must be ≫ 1 for the "many grains" averaging to hold — e.g. keep grains-on-target ≥ 100).
3. Recommendation on grain size and cloud diameter that keeps grains-on-target ≥ 100 *and* per-grain energy sub-catastrophic (INT-2) *and* launched mass affordable.

**Anti-mistake checks.** AM-1, AM-15. **Validation.** CW dispersion against analytic solution; ε against a direct Monte-Carlo grain-by-grain overlap count.

---

### INT-2 — Single grain → target coupling at sub-km/s (pin down β from physics)

**Purpose.** Determine the momentum-coupling coefficient β for one rhenium grain striking aluminum at the engagement speed, replacing the wide `U(1.0, 2.5)` envelope with a physically-bounded value for *this* regime.

**Questions answered.** At 0.5–1.1 km/s, does a rhenium grain embed, rebound, or crater-with-ejecta? What β does that imply? Is the impact even in the hypervelocity regime (it is not — see below), and does that make β closer to 1?

**Governing equations & regime.** Impact speed 0.5–1.1 km/s is **below the hypervelocity threshold** (~3 km/s / target sound-speed regime). In this regime momentum enhancement from ejecta is small: β is expected near 1.0–1.3, i.e. mostly plastic embedding/cratering, *not* the ejecta-dominated β up to 2.5 assumed for hypervelocity. Use a cratering/penetration model (e.g. modified Poncelet or an empirical low-velocity penetration relation) plus, for validation, a hydrocode spot-check (e.g. iSALE/Autodyn) at 2–3 design points. `SOURCED` + `DERIVED`

**Inputs.** Grain d_g = 50 µm–2 mm; v_rel = 500–1134 m/s; Re grain onto Al-6061 target; normal and oblique incidence (the cloud hits a curved target, so most grains strike obliquely).

**Outputs.**
1. β(v_rel, incidence) table with a **narrowed, defensible central value and range** for the design regime.
2. **Critical honesty flag:** if β ≈ 1.0–1.3 rather than 1.5, required incident mass rises ≈15–50% — recompute INT-3 and the consumable cost accordingly. This *darkens* the mass budget slightly and must be reported, not hidden by keeping the optimistic 1.5.
3. Per-grain crater volume / ejecta mass (feeds INT-5).
4. Rebound fraction (rebounding grains transfer up to 2× momentum but risk becoming new debris — quantify both effects).

**Anti-mistake checks.** AM-4a (β now moves from UNVALIDATED toward SOURCED *only* with a real low-velocity impact source or hydrocode run — otherwise stays UNVALIDATED), AM-16. **Validation.** Cross-check against published low-velocity metal-on-metal penetration data and the hydrocode spot-checks; note that the field's β data is overwhelmingly hypervelocity, so extrapolation into this regime is itself a stated limitation.

---

### INT-3 — Aggregate momentum delivery, mass efficiency & the imparted-Δv distribution

**Purpose.** Sum the grain-level coupling over the cloud to get the actual Δv imparted to the target per engagement, and its distribution (not a single number) — the export Part II's economics need.

**Questions answered.** How much Δv does one dispersed shot deliver? How much launched rhenium per engagement, including the ε penalty? How does Δv vary shot-to-shot given aim scatter and cloud non-uniformity?

**Governing equations.** `Δv_imparted = (Σ_hit m_grain β v_rel)/m_t = ε · m_L · β̄ · v_rel / m_t`, with ε from INT-1 and β̄ from INT-2. Launched mass to guarantee the design shift: `m_L = γ · m_t · Δv_target / (ε · β̄ · v_rel)`. `DERIVED`

**Inputs.** Δv_target = 53 m/s (200 km shift at 800–900 km) and the altitude-dependent values from v5 SIM-0; ε from INT-1; β̄ from INT-2; v_rel across the safe zone; γ = 2 margin on *launched* mass (not on per-grain energy — keeps the disruption gate intact, AM-16).

**Outputs.**
1. Table: h, v_rel, d_cloud, ε, β̄, **m_L (launched mass), incident mass, Re cost/shot**. Worked point: 1 cm @ 900 km, 5 cm cloud, ε = 0.04, β̄ = 1.15, v_rel = 619 ⇒ m_L ≈ 3 g, cost ≈ $22–30.
2. **Distribution of Δv_imparted** across many simulated shots (aim scatter + cloud lumpiness), feeding a *distribution* of orbit shifts and hence decay times into v5 SIM-2 — replacing the deterministic 53 m/s.
3. Confirmation that per-grain E_s stays < 40 kJ/kg while the *bulk* Δv is delivered (the AM-16 result, now earned rather than assumed).

**Anti-mistake checks.** AM-1, AM-2, AM-16, AM-18. **Validation.** Momentum conservation closes to machine precision; total launched-mass cost reconciles with the v5 consumable figure.

---

### INT-4 — Aim tolerance as a graceful mass penalty (the PI's point, formalized)

**Purpose.** Replace P_hit with the correct decision variable: how much *extra launched mass* buys a target reliability, given aim error.

**Governing equations.** For a Gaussian cloud of 1σ radius σ_cloud centered with aim error d_miss on a target of area A_t, the fraction of cloud mass on target is `ε(d_miss) = ε₀ · exp(−d_miss²/2σ_cloud²)` (small-target limit). To hold delivered Δv constant, launched mass scales as `m_L(d_miss) = m_L(0) · exp(+d_miss²/2σ_cloud²)` — smooth and finite, no cliff. `DERIVED`

**Inputs.** d_miss from INT-0 (≈5.8 mm central, plus the worst-case-cone value); σ_cloud from INT-1; reliability target (e.g. deliver ≥ design Δv on ≥ 95% of shots).

**Outputs.**
1. Plot: launched mass (and Re $/shot) vs. aim error, for several σ_cloud — the "graceful quadratic."
2. Optimal cloud diameter that minimizes launched mass at the required reliability (bigger cloud = more aim-tolerant but lower ε; there is a minimum).
3. Per-shot reliability at the tuned operating point, exported to Part II's yield model.

**Anti-mistake checks.** AM-5, AM-18. **Validation.** Compare against a grain-by-grain Monte-Carlo hit count.

---

### INT-5 — Secondary ejecta & net-debris check for the dispersed cloud

**Purpose.** Even sub-catastrophic, thousands of grain impacts crater the target and produce ejecta. Confirm net debris ≥ 0 for the *cloud* engagement, not just the single-impactor idealization.

**Governing equations.** Per-grain ejecta mass/size from INT-2 crater volume; aggregate over N_grains-on-target; classify ejecta by size vs. the 1 mm / 1 cm trackable thresholds; net = 1 (removed, eventually) − (trackable ejecta) − (rebounding grains that stay in orbit). `SOURCED`+`DERIVED`

**Inputs.** N_grains-on-target and crater model from INT-1/INT-2; target Al-6061; v_rel range.

**Outputs.**
1. Table: ejecta mass, N_ejecta(>1 mm), N_ejecta(>1 cm), rebounding-grain count, net-debris tally per engagement.
2. **Critical finding:** whether the many-small-impacts approach is genuinely net-positive, and the ejecta size below which fragments self-demise quickly (INT feeds v5 SIM-8/SIM-2).
3. Interim-hazard note (v5 R7): the target sits on a lowered but eccentric orbit for the decay interval — debit that on-orbit time in the net accounting.

**Anti-mistake checks.** AM-3, AM-8, AM-13. **Validation.** Cross-check aggregate ejecta against NASA SBM non-catastrophic (cratering) branch and any available multi-impact erosion data.

---

# PART II — The economically-viable engagement architecture (fine-tuning)

**The problem in one line:** cost/object = C_platform / N_engagements + consumables(~$30). Consumables are negligible; **the only lever is N_engagements per platform.** Everything here maximizes N at low Δv.

### The cost surface (from the completed run + INT costs)

| Platform cost | N=10³ | N=10⁴ | N=10⁵ | N=10⁶ |
|---|---|---|---|---|
| $50M | $50,031 | **$5,031** | **$531** | **$81** |
| $150M | $150,031 | $15,031 | $1,531 | $181 |
| $500M | $500,031 | $50,031 | $5,031 | $531 |

Benchmarks: ground-laser (theoretical) ~$1–5/object-equivalent for this size class; the "do nothing" alternative carries a collision-risk cost. **Reading:** a *defensible niche* (~$5k/object) needs N≈10⁴ on a $50M platform; *laser-competitive* economics (~$50–500/object) needs N≈10⁵–10⁶. The whole question is whether that many low-Δv targets are reachable.

---

### The reformulation that makes it possible: **cluster targeting**

Chasing random singletons is fatal: each new target on a different plane costs a plane-change Δv of hundreds–thousands of m/s, so N stays tiny and cost stays at $50k+/object. **The only economical architecture is to service co-orbital clusters** — debris that already shares an orbital plane and altitude, so (a) the platform matches the cluster mean orbit *once* and every member is automatically low-v_rel (satisfying the energy gate *and* the INT-0 co-orbital premise for free), and (b) moving between members costs only small intra-cluster Δv.

Real small-debris populations are dominated by exactly such clusters — major breakup clouds (e.g. Fengyun-1C, Cosmos-2251/Iridium-33, Cosmos-1408) each contribute tens of thousands to >100,000 untracked ≥1 cm fragments sharing a plane. **The dense 600–1000 km debris the PI flagged is largely this clustered breakup debris.** This reframes "10,000 objects per platform" from an arbitrary hope into a checkable question about specific clouds.

---

### ECO-1 (SIM-14) — Intra-cluster engagement cadence & Δv budget

**Purpose.** Compute how many cluster members one platform can engage per year, and the Δv it costs — the number that sets economic viability.

**Questions answered.** Given a debris cloud of known spread (in altitude, eccentricity, and RAAN drift from differential J2), how many members drift within engagement range of a station-keeping platform per day? What Δv/year to stay optimally placed and to occasionally reposition? How many engagements over a 10-year life?

**Governing equations.** Differential nodal precession `Δ(dΩ/dt) ≈ −(3/2)n J2 (R_E/a)² [cos i spread]` drives relative RAAN drift → members sweep past the platform's plane over months. Encounter rate = cluster spatial density × relative drift velocity × engagement cross-section. Repositioning Δv via low-thrust within the cluster. `DERIVED`

**Inputs.** Cluster models from ECO-2; engagement range R (INT-0 cone); platform propulsion (electric I_sp = 1,500–3,000 s); firing-geometry constraint from INT-0 (only engage within the validated retrograde cone at low v_transverse).

**Outputs.**
1. Encounters/day and engageable-fraction (those within the INT-0 cone at safe v_rel) vs. cluster density and platform station strategy.
2. N_engagements over 10 years, with Δv/year budget — the input to the cost surface.
3. Recommendation: station-keep-and-let-them-come vs. actively tour the cluster, whichever maximizes N per unit Δv.

**Anti-mistake checks.** AM-9, AM-11, AM-18. **Validation.** Propagate a synthetic cluster (differential J2 + drag) and count engagement windows directly; cross-check against published breakup-cloud evolution.

> **The central tension this SIM must resolve honestly:** the energy gate wants *low* v_rel (co-orbital), but a station that lets targets "come to it" sees them at *high* random relative velocity. The resolution is that within a *single breakup cloud* the members' velocity dispersion is small (tens–hundreds of m/s), so co-orbital low-v_rel engagement of many members is possible *inside one cloud* — but **not** across clouds. If ECO-2 shows individual clouds don't hold ≥10⁴ engageable members within the cone, the niche economics fail and this must be reported.

---

### ECO-2 (SIM-15) — Cluster availability in the real debris population

**Purpose.** Determine, from real catalog + statistical debris data, how many co-orbital ≥1 cm targets are actually reachable — the population that N is drawn from.

**Questions answered.** How many untracked ≥1 cm fragments does each major cluster hold within a band a platform can service at low Δv? How many clusters exist? Does any single cluster reach the 10⁴ (niche) or 10⁵–10⁶ (laser-competitive) thresholds?

**Governing equations / data.** ESA MASTER / NASA ORDEM spatial-density models for the ≥1 cm population by altitude/inclination; catalog breakup provenance for plane assignment. `SOURCED`

**Inputs.** MASTER/ORDEM density fields; major-breakup inventory; the engageable-band width from ECO-1.

**Outputs.**
1. Table: cluster, altitude, inclination, estimated ≥1 cm members within the engageable band, decay-benefit of a 200 km shift at that altitude (from v5 SIM-2).
2. **The viability-defining number:** max engageable co-orbital members per platform deployment, with uncertainty. Tag `UNVALIDATED` where it rests on statistical (not cataloged) fragment counts — most ≥1 cm objects are *not* individually tracked, so this number is itself a modeled estimate.
3. Cross-check: does the highest-density cluster also sit where the 200 km shift gives a real decay benefit (v5 R2)? Report clusters that are dense-but-too-high (shift useless) separately from dense-and-low (shift effective).

**Anti-mistake checks.** AM-11, AM-17. **Validation.** MASTER vs. ORDEM density cross-comparison; sanity-check total against published untracked-population estimates.

---

### ECO-3 (SIM-16) — Platform cost & consumable logistics

**Purpose.** Drive C_platform down and confirm the consumable/reload chain supports the required N.

**Questions answered.** What is the cheapest credible platform that carries the coilgun, sensors, and propulsion? How much consumable per platform-life, and does supply/reload support it? Rhenium or tungsten?

**Governing equations.** Mass budget (v5 SIM-11) → launch cost (v5 parameter table) + build cost; consumable mass = N × m_L (INT-3); reload cadence vs. payload capacity. `DERIVED`

**Inputs.** m_L per shot from INT-3 (~3 g); N from ECO-1/2; payload 10–100 kg per reload; Re vs. W (SIM-7 result).

**Outputs.**
1. C_platform range for a minimized design (target: push toward the $50M end — smaller platform, COTS avionics, electric propulsion, mass-produced if a fleet).
2. Consumable mass per platform-life: N=10⁴ ⇒ ~30 kg; N=10⁵ ⇒ ~300 kg; N=10⁶ ⇒ ~3,000 kg. **Decision:** at N ≥ 10⁵ the consumable mass and the 60–81 t/yr rhenium supply ceiling force a switch to **tungsten** (84,000 t/yr, 1/50–1/100 the price); at N=10⁶ even one platform's 3 t of rhenium is ~4% of world production — untenable for a fleet. Tungsten is the economically-required material, confirming SIM-7.
3. Reloads required over life and their Δv/cost (hub rendezvous vs. expendable-platform-per-cluster).

**Anti-mistake checks.** AM-9. **Validation.** Mass/cost cross-check vs. comparable smallsats; consumable supply vs. USGS production.

---

### ECO-4 (SIM-17) — Economic optimization & the honest viability gate

**Purpose.** Combine INT + ECO into cost/object, find the architecture that minimizes it, and state plainly whether it reaches viability — and against which benchmark.

**Governing equations.** `Cost/object = C_platform(ECO-3)/N(ECO-1,ECO-2) + m_L·P_material + energy + reload`, optimized over {cluster choice, altitude, station strategy, cloud diameter, v_rel, material}. `DERIVED`

**Outputs.**
1. **Tuned design card** (the deliverable the PI asked for): target cluster, altitude, Δh, v_rel, cloud diameter, material, m_L/shot, N over life, C_platform, **cost/object**, and the decay benefit achieved.
2. Cost/object vs. cluster density curve, with the $5k (niche) and $500 (laser-adjacent) lines drawn — showing which real clusters clear each.
3. **Two honest headline numbers, both reported:**
   - *Niche viability:* achievable if a single dense, low-altitude cluster supplies ≥10⁴ engageable members — cost/object ≈ $5k on a $50M platform, justified against avoided-collision value, **not** against lasers.
   - *Laser-competitive viability:* requires ≥10⁵–10⁶ members and ~27–274 engagements/day for a decade with tungsten consumable logistics — report whether ECO-2 finds any cluster that supports it (current expectation: **unlikely for a single platform; possibly reachable for a mass-produced fleet servicing the largest clouds**).
4. Binding-constraint ledger and the explicit list of `UNVALIDATED` conditions the chosen headline rests on (chiefly: the INT-0 co-orbital premise, the ECO-2 engageable-member count, and β̄ from INT-2).

**Anti-mistake checks.** AM-8, AM-9, AM-18. **Validation.** End-to-end arithmetic reconciliation with the v5 cost module; sensitivity to every `UNVALIDATED` input.

---

## PART III — The fine-tuned operating point (target design card)

To be *confirmed* by the SIMs above, not asserted — but this is the architecture the plan is tuned toward, and the concrete target:

| Lever | Tuned choice | Why | Confidence |
|---|---|---|---|
| Mission mode | **Service one dense breakup cluster** per deployment | Only way to keep every engagement co-orbital & low-Δv | `DERIVED` |
| Target cluster | Densest ≥1 cm cloud with real decay benefit (low-to-mid altitude, ~600–800 km) | Balances member count (ECO-2) vs. 200 km-shift decay benefit (v5 R2) | `UNVALIDATED` (ECO-2) |
| Altitude | ~600–800 km (not the 900 km v5 optimum, nor a "does-nothing" 500 km) | High cluster density *and* meaningful lifetime cut | `UNVALIDATED` (ECO-2) |
| v_rel | Upper safe-zone (~1,000–1,100 m/s) | Less launched mass per shot; still sub-catastrophic; check aim (INT-4) | `DERIVED` |
| Cloud diameter | Minimum that holds ≥95% shot reliability (INT-4) & ≥100 grains on target | Balances aim tolerance vs. mass efficiency ε | `DERIVED` |
| Material | **Tungsten** | Economically required at N≥10⁵; no momentum penalty | `SOURCED` |
| Platform | Minimized ~$50M, electric propulsion, COTS | Amortization is 99.4% of cost — this is the lever | `DERIVED` |
| N target | ≥10⁴ (niche) → ≥10⁵ (laser-adjacent) | Sets cost/object | `UNVALIDATED` (ECO-1/2) |

**Honest ceiling.** The physics is in hand; the economics live or die on one number — how many engageable co-orbital members a real cluster holds (ECO-2). If that is ≥10⁴, the concept is a **viable niche system** (~$5k/object) worth deploying against high-value-orbit collision risk. If it must reach 10⁵–10⁶ to beat lasers, that is a **fleet-scale, tungsten, mass-produced** proposition whose feasibility ECO-1/2 must decide — and which this plan does not pre-judge. Either way the result must state the number it rests on.

---

## PART IV — What must be true for economic viability (checkable conditions)

Report each as PASS/FAIL when the SIMs run; the concept is economically viable **iff all hold**:

1. **INT-0 co-orbital premise holds** within a usable approach cone (else the guidance wall and v_rel·τ cliff return, and per-shot mass/cost balloon).
2. **β̄ ≥ ~1.1** at engagement speed (INT-2) — if grains rebound-dominate or under-couple, mass/cost rise and possibly rebound-debris appears.
3. **A real cluster holds ≥10⁴ engageable ≥1 cm members** within a low-Δv, decay-beneficial band (ECO-1/ECO-2) — the make-or-break number.
4. **Intra-cluster Δv/year is affordable** on a minimized platform over 10 years (ECO-1/ECO-3).
5. **Net debris ≥ 0** including cloud ejecta and rebounds (INT-5), and interim on-orbit hazard (R7) is acceptable.
6. **Material = tungsten** at N≥10⁵ (Re supply ceiling), accepting rhenium only for a small niche demo.

**Risks that still darken it (carry from v5 + new):** the INT-0 reversal is the most attackable claim in the project (R1′); ECO-2's member count is a modeled, not cataloged, estimate (R2′); benchmarking against theoretical lasers is unfair *both* ways — they aren't deployed either, so "avoided-collision value" is the fairer frame but is itself uncertain (R8′); and calling it "rhenium" while the economics force tungsten is a naming/marketing gap the manuscript must own (R3′).

---

## PART V — Deliverables

1. **Interaction module** (INT-0…INT-5): resolved cloud–target coupling code, replacing v5 SIM-1/SIM-4, with the β̄ study and the ε/mass-penalty curves.
2. **Economic module** (ECO-1…ECO-4 / SIM-14…17): cadence, cluster-availability, platform-cost, and the optimization producing the design card and cost/object surface.
3. **Updated master + RQ-indexed results** folding these into the v5 pipeline (single source of truth; `run_all.py` regenerates).
4. **Two-number viability statement** (niche vs. laser-competitive), each with its binding `UNVALIDATED` condition named.
5. **Updated citation log** adding the low-velocity β source(s)/hydrocode runs (INT-2), the MASTER/ORDEM cluster data (ECO-2), and the dated tungsten pricing/production (ECO-3).

---

## Appendix — Module → question map

| Question | Module |
|---|---|
| How do the grains actually couple to the debris? | INT-1…INT-3 |
| Is the co-orbital / no-latency-cliff claim real? | INT-0 (validation gate) |
| What β applies at these speeds? | INT-2 |
| How much rhenium/tungsten per shot, and how forgiving is aim? | INT-3, INT-4 |
| Does the cloud create net debris? | INT-5 |
| How many objects can one platform service? | ECO-1, ECO-2 |
| What does it cost per object, and is that viable? | ECO-3, ECO-4 |
| What is the tuned, best-economics design? | Part III design card |
| What has to be true for it to work? | Part IV conditions |
