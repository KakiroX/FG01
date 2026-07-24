# FG01 / KIDR — Agent Context & Accelerated-Reentry Simulation Plan (v5)

**Project:** Rhenium-based momentum-transfer remediation of small (1–10 cm) LEO debris
**Platform:** FG01 (electromagnetic-launch "transporter" spacecraft)
**This revision (v5):** Pivot from *direct reentry* (perigee → 100 km, **falsified in v4**) to *accelerated reentry* (single targeted momentum kick that lowers the orbit ≈200 km, after which natural atmospheric drag completes the deorbit).
**Audience:** Every Claude agent that touches this project — simulation engineers, writers, reviewers.
**Supersedes:** "prompt 2" (SIM execution prompt v4) and the v4 falsification manuscript, which are retained as reference and history, not as the current design basis.

> **Read this first.** This document is the single source of truth for project context, physics, parameters, the anti-mistake contract, and the simulation plan. If anything you are asked to do conflicts with the *Anti-Mistake Contract* (Part 2) or the *Methodological Risks* register (Part 6), stop and flag it. The prior version of this work was invalidated by exactly the kinds of shortcuts those two sections forbid.

---

## PART 0 — Orientation for agents (read in 3 minutes)

- **What FG01 does.** An orbital platform carries an electromagnetic launcher (coilgun) and a rhenium payload. To remediate a tracked piece of small debris, FG01 launches a concentrated, briefly-dispersing cloud of rhenium grains onto a converging trajectory with the target. The hypervelocity encounter transfers momentum ("a kick"), lowering the target's orbit.
- **The mechanism this revision models (confirmed with the PI): *targeted momentum kick*, not a persistent global dust layer.** FG01 engages individual tracked objects; it does not seed a planet-wide drag cloud. This distinction is load-bearing — see Part 1.5 (why it matters for both the physics and the SOTA comparison).
- **The one-sentence pivot.** v4 tried to impart enough Δv to drop perigee straight to 100 km (immediate reentry) and was falsified four independent ways. v5 imparts only enough Δv to lower the orbit ≈200 km and lets atmospheric drag finish the job over months–years. The smaller Δv is the whole point: it is what lets the concept escape the prior study's *primary* failure.
- **The honest goal.** The PI wants the plan to *maximize* demonstrated viability and, ideally, identify an operating envelope where P(viable) → ~100%. That is a legitimate objective **only if pursued through falsification-first rigor**: we find the envelope where the gates genuinely pass, characterize it, and state plainly which assumptions a high-viability result rests on. A ~100% number that a reviewer can dismantle is worth less than an honestly bounded 60% number. Part 6 lists every issue that can pull the number down; do not bury them.
- **Traceability tags are mandatory on every number.** `DERIVED` (from stated physics), `SOURCED` (specific citation/DOI/page), or `UNVALIDATED` (assumption with no direct source at the relevant scale). No bare point values in anything destined for the manuscript.

---

## PART 1 — Project context

### 1.1 The concept and the hardware

Small orbital debris (≈1–10 cm) is below the ~10 cm ground-tracking threshold in LEO yet carries enough mass and closing velocity (7–15 km/s in random encounters) to destroy an operational spacecraft. Conventional ADR (robotic capture, nets, harpoons, tethers) needs a tracked, cooperative-enough target and does not scale to this size class. FG01's premise is momentum transfer from a launched projectile mass, so the object need only be *detected and briefly tracked terminally*, not cataloged and rendezvoused.

Rhenium was originally chosen for its density (21.01 g/cm³, fourth-densest stable element). **This justification is physically weak** and is a standing vulnerability (Part 6): momentum p = mv is independent of density at fixed mass; density changes only grain size, not momentum delivered.

### 1.2 What changed from v4 → v5: direct vs. accelerated reentry

| | v4 — Direct reentry (FALSIFIED) | v5 — Accelerated reentry (this plan) |
|---|---|---|
| Deorbit target | Perigee lowered to 100 km → immediate reentry | Orbit lowered ≈200 km → natural decay finishes |
| Δv required (800 km) | ≈194 m/s `DERIVED` | ≈53 m/s `DERIVED` |
| Design-point specific energy E_s (β=1.5, γ=2, v_rel=619 m/s) | ≈80 kJ/kg (2× catastrophic threshold) | ≈22 kJ/kg (**below** 40 kJ/kg threshold) |
| Reentry timescale | Minutes | Months–years (drag-limited) |
| Success criterion | Perigee < 100 km | Post-kick orbital lifetime < 25 yr (stretch: < 5 yr) |

The success criterion moves from "did it reenter immediately?" to "did the ≈200 km shift cut the natural orbital lifetime below the regulatory limit?" This makes **orbital-decay modeling (SIM-2) the new center of gravity** of the whole analysis, not a downstream check.

### 1.3 The prior study's four independent falsifications (do not silently reintroduce)

v4 concluded 0.0000% joint viability (10⁶-sample Monte Carlo). The four independent bases were:

1. **Structural margin–disruption finding (primary).** Because both deorbit impulse and disruption risk scale with delivered momentum, E_s = Δv_deorbit · v_rel / (2β). Any finite mass safety margin (γ=2) pushed the design point to exactly 2× the NASA Standard Breakup Model (SBM) 40 kJ/kg catastrophic threshold. A design point with both adequate margin *and* sub-catastrophic energy did not exist.
2. **Terminal guidance collapse.** Coupling the Clohessy–Wiltshire cloud geometry into the hit-probability model, P_hit fell from 0.268 at 1 ms control lag to ≈6×10⁻⁵⁸ at 10 ms — because positional drift v_rel·τ (≈6.2 m) dwarfed the cloud cross-section (σ_cloud ≈ 0.38 m). Control lag drove 68% of total output variance (Sobol).
3. **Net debris increase.** In the catastrophic regime the kick generated ≈443 trackable fragments per engagement: net −442 objects.
4. **Cost.** Best-case gross ≈$163,500/kg engaged, exceeding every competing method before the (formally undefined) fragmentation penalty.

**How the v5 pivot maps onto these:** the 200 km shift *directly defeats #1 and #3* by cutting Δv (and therefore E_s) below the catastrophic threshold — see 1.4. It **does not by itself defeat #2 or #4**, which are largely independent of Δv. That asymmetry is the honest core of this plan and dictates where the simulation effort must concentrate.

### 1.4 Why the pivot works, quantitatively (the viability lever)

Δv to lower perigee by 200 km via a single apogee burn (vis-viva), and resulting design-point E_s at v_rel = 619 m/s, β = 1.5, γ = 2 `DERIVED`:

| Altitude (km) | Δv direct→100 km (m/s) | Δv shift −200 km (m/s) | E_s design, direct (kJ/kg) | E_s design, shift (kJ/kg) | v_max,safe shift, β=1.5 (m/s) |
|---|---|---|---|---|---|
| 400 | 87.4 | 57.7 | 36.1 | 23.8 | 2,079 |
| 600 | 142.0 | 55.2 | 58.6 | 22.8 | 2,173 |
| 800 | 193.8 | 52.9 | 80.0 | 21.8 | 2,268 |
| 1000 | 243.1 | 50.7 | 100.3 | 20.9 | 2,365 |
| 1200 | 290.0 | 48.7 | 119.7 | 20.1 | 2,463 |

**Read-out:** across the entire LEO band the 200 km-shift design point sits at ≈20–24 kJ/kg — comfortably below the 40 kJ/kg catastrophic gate — and the safe intercept-velocity ceiling rises from the v4 range of ~410–620 m/s to ~1,400–2,500 m/s (β = 1.0–1.5). **The structural falsification of v4 is genuinely lifted by the pivot.** This is the plan's strongest, most defensible result and should be established first (SIM-0).

**But note the built-in tension the plan must expose, not hide:** a 200 km drop is a *fixed altitude decrement*, so its effect on *decay time* is altitude-dependent and highly nonlinear. From 600 km → 400 km lands in the years-to-months decay regime (density rises ~10× per ~150–200 km); from 1000 km → 800 km leaves the object in a multi-decade-to-century regime that likely still fails the 25-yr rule. Debris is densest at 600–1000 km (PI note), which is exactly where the decrement is least sufficient. The "optimal altitude" is therefore a real optimization, not a given — see SIM-2 and SIM-13.

### 1.5 SOTA landscape and where FG01 sits (2024–2026)

- **The drag-dust-cloud prior art is tungsten, and it was not adopted.** Ganguli et al. (US Naval Research Laboratory, 2011) proposed injecting ~20–40 t of **tungsten** micro-dust to raise drag above ~900 km and decay untracked debris over ~25 yr. It stalled on three objections that FG01 must be measured against: indiscriminate drag also decays *operational* satellites; the dust itself is new debris and can clump; and a metal layer degrades astronomical observation at some wavelengths. FG01's *targeted* kick avoids the indiscriminate-drag and global-astronomy objections — but re-inherits the targeting/guidance problem that the global cloud sidesteps. This is the fundamental trade: **untargeted → indiscriminate harm; targeted → the guidance wall (falsification #2).**
- **Operational ADR today** targets large, trackable objects: Astroscale ADRAS-J (launched 2024, characterized a rocket upper stage at ~50 m station-keeping); ESA/ClearSpace-1 (robotic-arm capture). Neither addresses the untracked 1–10 cm class FG01 targets — that gap is FG01's actual niche and should be stated as the motivation.
- **Competing small-debris concepts & published costs** (for SIM-9): ground-based laser ablation ~$100–500/kg (theoretical); electrodynamic tether ~$1,000–5,000/kg (demonstrated); aerogel capture ~$5,000–20,000/kg; ADR tug ~$10,000–100,000/kg. FG01 must be benchmarked against these as ranges, with sources.
- **Models that are the field's accepted currency** (use these, not hand-fits): NASA Standard Breakup Model (Johnson et al., 2001) for fragmentation and the 40 J/g catastrophic criterion; NRLMSISE-00 for atmospheric density/decay; Clohessy–Wiltshire relative motion for cloud dispersion; Fay–Riddell stagnation heating for reentry ablation; USGS Mineral Commodity Summaries for material prices.

### 1.6 Research-aims map (from *FG01 techs.md*)

Every FG01 techs.md aim is assigned a primary SIM in Part 8. The three top-level aims are: **(1)** prove the momentum-kick rhenium cloud can deorbit small debris (mode, slow-down, altitude ceiling, launch angle, mass estimate f(h), mass/velocity ratio, shot precision, and whether it creates more junk); **(2)** cost-efficiency vs. competitors (energy per shot, coilgun design, $/kg deorbited, cost of misses); **(3)** is it acceptable to deorbit — environment, legality, aviation-strike tracking, and whether trash and rhenium both fully ablate.

---

## PART 2 — Anti-Mistake Contract (binding)

Breaking any rule reproduces an error that already invalidated a prior version. These carry over from v4 and add the v5-specific gates.

**AM-1 — No assumed efficiency.** Never insert an unsourced interception-efficiency η or momentum-coupling factor to make a number work. Required projectile mass comes from momentum balance; hit probability comes from the coupled cloud/guidance model.
**AM-2 — Solid-sphere target mass.** Debris mass = solid Al-6061 sphere (ρ = 2.70 g/cm³). State it; do not quietly use a hollow or low-density shortcut that flatters mass ratios.
**AM-3 — SBM gate before design.** Evaluate the catastrophic-disruption boundary before any projectile/guidance/platform sizing. If no design point satisfies it, that is a reportable result.
**AM-4a — β is UNVALIDATED.** No confirmed Al-on-Al, cm-scale, sub-orbital-velocity momentum-enhancement source exists. Use a flat U(1.0, 2.5) envelope; run every value; never collapse to a point. A regression may reappear *only* with a real DOI/figure. A secondary context run up to β≈12 (strengthless-target bound) may be reported, clearly labeled "context, not design basis."
**AM-5 — No hidden substitutions in guidance.** Report P_hit as a probability, never as a guaranteed hit.
**AM-6 — Real NRLMSISE-00.** Use an actual implementation with stated version; no fitted density curve. State F10.7 / Ap driver assumptions.
**AM-7 — Rate-limited reentry heating.** Fay–Riddell q_w = k√(ρ_∞/R_n)·v_∞³, not a total-energy ablation budget.
**AM-8 — Misses and failures are in the cost.** Cost per kg removed includes missed shots and the fragmentation penalty; if net removal ≤ 0, report cost-per-kg-removed as undefined, not a large finite number.
**AM-9 — No aspirational tech/costs as current.** Demonstrated launcher benchmarks only; Starship-class launch prices only as explicitly labeled sensitivity.
**AM-10 — Report envelopes, not point estimates,** wherever Part 3 marks a range.
**AM-11 — Single source of truth.** Downstream SIMs read decay/fragment/lifetime numbers from upstream SIM outputs at runtime, not re-typed constants.
**AM-12 — Honest UQ.** Monte Carlo over joint distributions (not one-at-a-time sweeps); report convergence and Sobol indices.
**AM-13 — SBM mass term.** Verify the effective-mass definition against Johnson et al. (2001) directly (total colliding mass vs. KE-based) before coding fragment counts.
**AM-14 — Couple cloud → guidance.** P_hit uses SIM-3's computed cloud cross-section σ_cloud(t), not an independent assumed cone angle.
**AM-15 — Swarm ≠ single impactor (flag every time).** Treating N micro-grain impacts as one macroscopic impactor of equal total mass is an UNVALIDATED assumption. In v5 this cuts both ways (see AM-16); label it in every output that depends on it.
**AM-16 — [v5] Per-impact vs. bulk energy must be reported separately.** The catastrophic-disruption gate is a *per-impact specific-energy* criterion. A dispersed grain swarm may deliver the *bulk* momentum needed for the 200 km shift while keeping *each* grain's specific energy sub-catastrophic. This is potentially the concept's second escape hatch — **and** an unvalidated claim. Never report the bulk-momentum success without simultaneously reporting the per-impact energy and flagging AM-15/AM-16.
**AM-17 — [v5] Decay-time is a hard gate, not a footnote.** "Accelerated reentry" only counts as success if the post-kick lifetime clears the regulatory bar. Report T_decay for the *shifted* orbit at every altitude and never assume a 200 km drop is sufficient; prove it per altitude.
**AM-18 — [v5] "Maximize viability" ≠ cherry-pick.** Selecting the best operating point is legitimate; hiding the altitudes/velocities where gates fail is not. Every headline viability figure must state the envelope it holds within and the constraints binding outside it.

---

## PART 3 — Honest parameter table (2026-current, with traceability)

Run across the full range wherever a range is given; report the envelope.

| Parameter | Value / range | Units | Traceability |
|---|---|---|---|
| μ (Earth GM) | 398,600.4418 | km³/s² | DERIVED |
| R_E (mean radius) | 6,371 | km | DERIVED |
| ρ_Al (6061) | 2.70 | g/cm³ | SOURCED (USGS/matls) |
| ρ_Re | 21.01 | g/cm³ | SOURCED |
| ρ_W (tungsten, comparator) | 19.25 | g/cm³ | SOURCED |
| E_s,c (SBM catastrophic threshold) | 40 (30–50) | kJ/kg (= J/g) | SOURCED (Johnson 2001) |
| β (momentum enhancement) | 1.0–2.5 (flat) | — | **UNVALIDATED** |
| γ (mass safety margin) | 2.0 | — | DERIVED |
| Orbit shift target Δh | 200 (150–250) | km | DESIGN VARIABLE |
| Altitude sweep h | 400, 500, 600, 700, 800, 900, 1000, 1100, 1200 | km | DESIGN VARIABLE |
| v_rel (intercept) | constrained to SIM-0 safe zone per h | m/s | DERIVED bound |
| P_Re (rhenium, 2026 spot) | ≈7,283 (2026-07-22 spot); 6,389 (2026-Q1 end); 2,486 (2025 open) | USD/kg | SOURCED (market, dated) |
| Re production | ≈60–81 | t/yr | SOURCED |
| P_W (tungsten) | 30–100 | USD/kg | SOURCED |
| W production | ≈84,000 | t/yr | SOURCED |
| C_launch | 3,000 (dedicated-full) – 7,000 (rideshare-marginal) | USD/kg | SOURCED, dated; Starship-class 100–500 = labeled sensitivity only |
| P_elec | 0.10 | USD/kWh | SOURCED |
| Control-loop latency τ | 0.1, 1, 5, 10, 50, 100 | ms | DESIGN VARIABLE (dominant risk) |
| F10.7 (solar) | 70 / 150 / 250 (min/nom/max) | sfu | SOURCED |
| H_eff ablation Re / Al | 4.5 (3.5–5.5) / 12 (10–14) | MJ/kg | SOURCED (range) |
| Casualty threshold (NASA) | 1×10⁻⁴ per event | — | SOURCED |
| FCC disposal rule | 5 | yr | SOURCED (47 CFR 25.283) |
| NASA ODMSP | 25 | yr | SOURCED |

> **Note on rhenium price:** it roughly tripled from 2025 to mid-2026. Always report which date/edition a figure is from; do not present $2,600 (2025 USGS average) as "current."

---

## PART 4 — Simulation plan

Run in the order given. **SIM-0, SIM-2, and SIM-4 are the three gates that decide viability; everything else characterizes or costs the survivors.** For each SIM: *Purpose · Questions answered · Governing equations · Inputs · Outputs · Anti-mistake checks · Validation.* This structure mirrors the v4 execution prompt so results drop straight into the manuscript.

The pipeline dependency is: **SIM-0 → SIM-1 → SIM-2 (gate) → SIM-3 → SIM-4 (gate) → {SIM-5 … SIM-11} → SIM-12 (UQ) → SIM-13 (operating-point selection).**

---

### SIM-0 — Catastrophic-disruption boundary for the 200 km-shift kick (the re-opened gate)

**Purpose.** Redraw the v4 falsification gate with Δv_shift (≈200 km drop) instead of Δv_deorbit (→100 km), and show explicitly where the design point now sits relative to the 40 kJ/kg catastrophic threshold. This is the result that lifts v4's primary falsification — establish it first and cleanly.

**Questions answered.** "Prove the momentum-kick cloud is viable" → does a non-catastrophic (h, v_rel) region exist with adequate margin? "Can it slow debris enough?" → does the safe-zone Δv meet the 200 km-shift requirement (it does by construction; verify). "How high an orbit can we target?" → the *energy* gate is permissive at all altitudes; the *binding* ceiling comes from SIM-2, not here — state this hand-off.

**Governing equations.**
- Δv for a 200 km perigee drop via apogee burn (vis-viva): `Δv_shift(h) = √(μ/(R_E+h)) − √{μ[2/(R_E+h) − 1/a']}`, with a' = R_E + h − Δh/2 for the post-burn ellipse (apogee R_E+h, perigee R_E+h−Δh). `DERIVED`
- Momentum balance: m_t·Δv_shift = β·m_p·v_rel. `DERIVED`
- Design specific energy: `E_s,design = γ·Δv_shift(h)·v_rel / (2β)`. `DERIVED`
- Margin-free safe velocity ceiling (same convention as v4 Eq. 4): `v_max,safe(h,β) = 2·E_s,c·β / Δv_shift(h)` — this is what the Part 1.4 table reports. The margin-carrying design (γ=2) remains sub-catastrophic for v_rel up to **v_max,safe/γ**; at 800 km, β=1.5 that is ≈1,134 m/s, still well above the ≈619 m/s used in the E_s,design column. `DERIVED`

**Inputs.** h = 400…1200 km (step 100); Δh = 150, 200, 250 km; β = full flat 1.0–2.5; E_s,c = 30, 40, 50 kJ/kg; γ = 2; v_rel = candidate 500–2,500 m/s.

**Outputs.**
1. Table: h, Δv_shift, E_s,design, v_max,safe for each (β, E_s,c, Δh).
2. Plot: v_max,safe(h) with β as parameter, SAFE/CATASTROPHIC shading, one panel per E_s,c — overlaid with the demonstrated EM-launch envelope (SIM-6) and the GNC feasibility limit (SIM-4).
3. Plot: E_s,design(h) for direct-reentry (v4) vs. 200 km-shift (v5), same axes, threshold line at 40 kJ/kg — this single figure is the headline of the pivot.
4. Explicit conclusion as a range: "For β ∈ [1.0, 2.5] and Δh = 200 km, the non-catastrophic regime exists at **every** tested altitude with v_rel up to X–Y m/s; the design point sits at 20–24 kJ/kg, i.e. 0.5–0.6× the catastrophic threshold, reversing the v4 structural finding."

**Anti-mistake checks.** AM-3, AM-4a, AM-10, AM-17, AM-18. **Validation.** Verify Δv_shift against a Hohmann two-burn cross-check and against GMAT; verify the E_s substitution algebraically (show it, don't assert it).

---

### SIM-1 — Momentum transfer & mass scaling f(h) (no assumed efficiency)

**Purpose.** Compute rhenium mass required per kg of debris for the 200 km shift, from physics only.

**Questions answered.** "Velocity/mass ratio between Re and debris?" "Does the ratio depend on distance?" → depends on h (via Δv_shift) and v_rel (via β), **not** on engagement range R (R only affects whether you hit — SIM-4's job). "Grams of Re per gram of Al?" → f(h) = m_p,design/m_t. "Estimate Re on a given orbit, f(x), x = height?" → tabulate f(h).

**Governing equations.** `m_p,design = γ·Δv_shift(h)·m_t/(β·v_rel)`; `f(h) = m_p,design/m_t = γ·Δv_shift(h)/(β·v_rel)` (target-mass-independent); grain count/size from ρ_Re. `DERIVED`

**Inputs.** m_t = (4/3)πr³·ρ_Al for d = 0.5, 1, 5, 10 cm; v_rel in SIM-0 SAFE zone per h; β flat 1.0–2.5 (primary); appendix β up to 12 (context only).

**Outputs.**
1. Table: h, v_rel, β, m_p,min, m_p,design, E_s (consistency-check column), f(h).
2. Plot: f(h) vs. h across the β envelope (band, not line) — compare v5 vs. v4 magnitudes (v5 should need ≈3–4× *less* Re per kg because Δv is ≈3–4× smaller).
3. Plot: minimum grain diameter vs. h (from m_p,design and ρ_Re), with the same for tungsten overlaid.
4. Range statement: "At β ∈ [1.0, 2.5], h = 800 km, Δh = 200 km, deorbiting 1 kg Al requires X–Y g Re at v_rel = Z m/s."
5. **Required AM-15/AM-16 limitation statement** in every output: single-macroscopic-impactor vs. N-grain-swarm equivalence is unvalidated.

**Anti-mistake checks.** AM-1, AM-2, AM-3, AM-4a, AM-10, AM-15, AM-16. **Validation.** Dimensional check; verify E_s < E_s,c at all reported points (consistency check).

---

### SIM-2 — Orbital decay & the accelerated-reentry success gate (NRLMSISE-00) — **new center of gravity**

**Purpose.** Determine whether a 200 km shift actually cuts the natural orbital lifetime below the regulatory bar, per altitude. In v5 this is the *primary* physical success criterion, not a downstream check.

**Questions answered.** "Slow debris enough that it reenters?" → yes iff T_decay(shifted orbit) < 25 yr (stretch < 5 yr). "How high an orbit can we target?" → the highest h at which the shifted orbit still clears the bar — **this defines the real altitude ceiling.** "Will missed rhenium decay?" → T_decay for missed grains on their un-shifted orbit. "Where is the optimal altitude?" → feeds SIM-13.

**Governing equations.** Near-circular decay `da/dt = −C_D (A/m) ρ(h)√(μa)`; full 3-D Cartesian propagation (with atmospheric co-rotation) for the eccentric post-kick orbit that sweeps orders of magnitude in density per revolution; ρ(h) from NRLMSISE-00. `SOURCED`

**Inputs.** Pre-kick: circular at h. Post-kick: apogee = h, perigee = h−Δh (Δh = 150/200/250). Debris: solid Al spheres d = 1, 2, 5, 10 cm (compute A/m). Missed Re grains: d = 0.05, 0.1, 0.5, 1, 2 mm. F10.7 = 70/150/250. NRLMSISE-00 (state implementation + version).

**Outputs.**
1. Table: h, d_debris, Δh, T_decay (yr) for each F10.7 — pre-kick vs. post-kick, and the **ratio** (the acceleration factor).
2. Table: h, d_grain, T_decay for missed rhenium (and tungsten comparator).
3. Plot: T_decay vs. h for shifted debris, with the 5-yr and 25-yr lines drawn — the crossing points are the altitude ceilings.
4. **Critical findings, stated plainly:** (a) "A 200 km shift reduces lifetime by a factor of ≈N at altitude h." (b) "The shift clears the 25-yr rule only below h ≈ X km; above that, accelerated reentry as specified does **not** achieve compliant decay." (c) "Missed Re grains ≥ d mm at h km persist Y yr — violates the FCC 5-yr rule if Y > 5."

**Anti-mistake checks.** AM-6, AM-11, AM-17. **Validation.** Cross-check decay rates against TLE-derived histories from Space-Track for comparable A/m (v4 used Starshine-1, NORAD 25769, 387 km, 258-day observed lifetime — reuse and report the % error).

---

### SIM-3 — Cloud dispersion, evolution & the cross-section that SIM-4 consumes

**Purpose.** Replace any assumed cone angle with orbital mechanics; produce σ_cloud(t) as the exported input to SIM-4. In v5 the cloud's finite cross-section is a *deliberate viability lever* for guidance — model it as such.

**Questions answered.** "Dust cloud or grains? If grains, how small? If cloud, how to predict it?" "Will cloud behavior influence removal accuracy?" → quantitatively, via σ_cloud(t) into SIM-4.

**Governing equations.** Clohessy–Wiltshire relative motion in the Hill frame:
`δz(t)=(δv/n)sin(nt)`, `δx(t)=(2δv/n)(1−cos nt)`, `δy(t)=(3δv/n)(nt−sin nt)`, with n=√(μ/a³); σ_cloud = cloud angular half-width evaluated at the time-of-flight for engagement range R. `DERIVED` (δv_launch UNVALIDATED — no basis for converting a solid slug into a uniform swarm at these scales; flag it).

**Inputs.** Dust mode: N large, m_p ≈ 0.01 mg per grain. Grain mode: N = 10³–10⁶, m_p = 1–100 mg. δv_launch = 0.1–10 mm/s. h = 400, 600, 800, 1000 km. Initial cone half-angle θ₀ = 1°, 5°, 10° (initial condition only).

**Outputs.**
1. Plot: cloud dimensions (δy, δz) vs. time, dust vs. grain mode.
2. Table: cloud number density vs. time.
3. Table: **σ_cloud(t) at the exact times-of-flight SIM-4 uses (R = 1, 5, 10, 50, 100 m)** — the mandated export.
4. **The v5 tension, quantified:** areal grain density vs. cloud cross-section — a bigger cloud raises P_hit but lowers momentum-per-transit and can drop per-grain delivery below what SIM-1 needs. Report the Pareto curve.
5. Recommendation on mode and grain size, justified against (a) persistence (SIM-2), (b) interception probability (SIM-4), (c) coilgun feasibility (SIM-6), (d) per-impact energy staying sub-catastrophic (SIM-5/AM-16).

**Anti-mistake checks.** AM-1, AM-14, AM-15. **Validation.** CW equations checked against analytical solutions.

---

### SIM-4 — Terminal GNC, launch angle & shot precision (the gate the pivot does NOT auto-solve)

**Purpose.** Honestly assess detection, tracking, and interception of cm-scale debris using SIM-3's real σ_cloud(t). This was v4's dominant failure (68% of variance); the pivot does not remove it, so **this is where viability is won or lost in v5.**

**Questions answered.** "What launch angle?" → retrograde θ that maximizes v_rel within the SIM-0 safe zone (which is now much wider). "What guarantees accuracy?" → nothing; report P_hit as a probability. "How precise?" → d_miss in metres. "Does the error matter, and how much does it cut efficiency?" → η_eff = P_hit, propagated to SIM-9/SIM-12.

**Governing equations.** `σ_pos(R)=(λ/D)·R`; `d_miss=√(σ_pos² + (v_rel·τ)²)`; `P_hit(R)=exp[−d_miss²/(2σ_cloud²)]`, with σ_cloud from SIM-3. `DERIVED`

**Inputs.** Target: 1 cm Al sphere (m_t = 1.41 g). Sensors: optical (D = 10 cm, λ = 550 nm), LIDAR (1 MHz BW, SNR = 10). τ = 0.1, 1, 5, 10, 50, 100 ms. R = 1, 5, 10, 50, 100 m. v_rel from SIM-0 safe zone. σ_cloud primarily from SIM-3 (AM-14); fixed-angle version reported alongside as a bounding comparison only.

**Outputs.**
1. Table: sensor, max detection range, σ_pos, TOF, d_miss, P_hit (CW-derived primary and fixed-angle secondary, side by side).
2. Plot: P_hit vs. R for each sensor and τ; **and P_hit vs. τ** with the 1 ms cliff shown.
3. **The v5 viability question answered explicitly:** "Given the wider safe-velocity envelope, is there a (v_rel, σ_cloud, τ, R) combination with P_hit > 0.5 at a *realistically achievable* τ? If the only solutions need τ < 1 ms, the concept remains guidance-limited despite the pivot." Report the minimum τ and maximum σ_cloud that yield P_hit ≥ 0.5 and ≥ 0.9.
4. Sensitivity: how much does enlarging σ_cloud (SIM-3) buy in τ-tolerance, and at what momentum-density cost (SIM-1)?

**Anti-mistake checks.** AM-1, AM-5, AM-14, AM-18. **Validation.** Compare against demonstrated rendezvous/tracking limits (OSAM, ELSA-d, ADRAS-J proximity ops). Use GMAT only to verify Lambert/relative-motion algebra, never to validate GNC performance.

---

### SIM-5 — Fragmentation & secondary-debris hazard (now mostly a confirmation, plus the swarm question)

**Purpose.** Confirm that in the v5 sub-catastrophic regime the kick does **not** generate net debris, and quantify the residual per-grain cratering/ejecta that even a sub-catastrophic swarm produces.

**Questions answered.** "Will it create more space junk?" → Net = removed − fragments − missed Re. In v5 this should flip positive at the sub-catastrophic design point — verify, don't assume.

**Governing equations.** NASA SBM cumulative count `N(L_c)=0.1·M^0.75·L_c^−1.71 (L_c>1 mm)`; effective mass M per Johnson 2001 regime rules (AM-13). For sub-catastrophic impacts use the cratering/ejecta (non-catastrophic) branch, **not** the total-mass catastrophic branch — this correction is central to v5. `SOURCED`

**Inputs.** E_s from SIM-1 per design point across β; m_t solid sphere; fragment sizes L_c = 1 mm, 1 cm, 5 cm, 10 cm. Explicitly include the AM-16 swarm case: N grains each at sub-catastrophic per-impact E_s.

**Outputs.**
1. Table: E_s, regime (with source equation noted), N(>1 cm), N(>1 mm), ejecta mass, net debris flux.
2. Plot: net debris flux vs. h for different v_rel — v5 expected positive below the SIM-2 ceiling.
3. **Critical finding:** "At the v5 sub-catastrophic design point the kick removes 1 object and generates F fragments ≥ 1 cm; Net = +/− Y." Plus: "Even sub-catastrophic, each hypervelocity grain produces ejecta; total secondary mass = Z g/engagement."

**Anti-mistake checks.** AM-3, AM-8, AM-11, AM-13, AM-15, AM-16. **Validation.** Cross-check against NASA EVOLVE 4.0 / ESA MASTER for known breakup events; verify the non-catastrophic branch is the correct one at v5 energies.

---

### SIM-6 — Electromagnetic launch technology assessment

**Purpose.** Determine the achievable (m_p, v) envelope from demonstrated tech; do not assume beyond it.

**Questions answered.** "Electrical energy per max-capacity shot?" "Exact coilgun dimensions — coils, copper thickness?" "Energy to charge to reach max velocity?"

**Governing equations.** E_kin = ½m_p v²; E_elec = E_kin/η; standard coilgun stage/inductance relations. `DERIVED`

**Inputs.** m_p = 0.01, 0.1, 1, 10, 100 g; v = 500, 1000, 2000, 3000, 5000 m/s; η = 0.3, 0.4, 0.5; L' = 10–100 µH/m; L_stage = 0.2–0.5 m. Demonstrated benchmarks (cite, dated): US Navy/ONR 3.2 kg @ 2,520–2,579 m/s (10.16–10.64 MJ, 2008); 3 g @ 5,900 m/s; 10 g @ ≈2,000 m/s.

**Outputs.**
1. Design table: m_p, v, E_elec, L_b, N_stage, wire gauge (AWG), m_cap, m_Cu.
2. Gap analysis: "Required (from SIM-1): X. Demonstrated: Y. Gap: Z×." (v5 velocities are *lower* than v4's, so this should be comfortably inside the envelope — state it.)
3. Plot: achievable (m_p, v) envelope vs. required from SIM-1.

**Anti-mistake checks.** AM-9. **Validation.** Cross-check against Turman (1996) and McNab (2003).

---

### SIM-7 — Material trade study (rhenium vs. tungsten — a standing vulnerability)

**Purpose.** Determine whether rhenium is actually justified, or whether tungsten (the NRL prior-art choice) dominates it. Report the honest answer even though the project is named for rhenium.

**Questions answered.** "Why rhenium?" → compare on the metrics that matter, not density alone.

**Metrics.** (1) Momentum efficiency p = mv — identical across materials at fixed mass (density is irrelevant here; this corrects the concept's founding justification). (2) Cost USD/kg. (3) Availability t/yr. (4) Oxide toxicity/volatility (Re₂O₇ volatile vs. WO₃ vs. Ta₂O₅). (5) Orbital persistence of missed grains — A/m = 3/(4ρr); **higher density → smaller grain for a given mass → potentially faster decay of missed material.** This is the *one* axis where rhenium may legitimately beat tungsten; quantify it. (6) Ablation energy H_eff (higher = better reentry survival = *worse* for missed material).

**Inputs.** Re (21.01 g/cm³, $2,600–7,283/kg, 60–81 t/yr, Re₂O₇ volatile); W (19.25, $30–100/kg, 84,000 t/yr, WO₃); Ta; steel; Cu. Prices from USGS Mineral Commodity Summaries (state edition) plus dated 2026 market spot.

**Outputs.**
1. Ranking table across all six metrics.
2. Weighted score (weights stated) or Pareto frontier.
3. **Conclusion, stated honestly:** "Rhenium is/is not justified. On momentum, cost, and availability, tungsten dominates by ≈50–100× on cost and ≈1,000× on availability with a 2.4% density penalty. Rhenium's only defensible edge is [smaller/denser missed-grain → faster decay], worth X if and only if that margin survives SIM-2." If tungsten dominates, say so — a reviewer will, and the v4 study already did.

**Anti-mistake checks.** AM-9, AM-10. **Validation.** Prices from USGS 2025/2026 (state which edition each figure comes from) + dated market data.

---

### SIM-8 — Reentry ablation, survival & environmental impact

**Purpose.** Determine whether debris and rhenium fully ablate, and assess deposition.

**Questions answered.** "Will trash fully burn up?" "Will rhenium fully burn up?" "Environmental effect?" → deposition mass flux. "Is it okay to deorbit?" → casualty expectation vs. NASA threshold.

**Governing equations.** Rate-limited Fay–Riddell q_w = k√(ρ_∞/R_n)·v_∞³ (AM-7); ablated fraction, surviving ballistic coefficient, ground-casualty expectation. `SOURCED`

**Inputs.** v_∞ = 7.5 km/s; d = 0.1, 0.5, 1, 2, 5 cm for Re and Al; H_eff,Re = 4.5 (3.5–5.5), H_eff,Al = 12 (10–14) MJ/kg; NRLMSISE-00 density during reentry; ρ_pop = 0 (ocean) to 10,000/km² (urban).

**Outputs.**
1. Table: d, material, ballistic coefficient B, f_ablate, survives? (Y/N), ground-impact probability.
2. Environmental assessment: annual Re/Al deposition mass flux from missed material vs. natural cosmic-dust background (v4 estimate ≈1 kg/yr Re₂O₇, orders of magnitude below background — reconfirm).
3. Note on Re₂O₇ volatility (m.p. 297 °C, sublimes) — assess oxidation fraction during reentry.
4. Casualty expectation per deorbit event vs. NASA 1×10⁻⁴ threshold.

**Anti-mistake checks.** AM-7, AM-11. **Validation.** Compare against NASA DAS and ORSAT.

---

### SIM-9 — Cost-economic model & competitor comparison

**Purpose.** Honest cost per kg removed, including misses, using a realistic launch-cost band.

**Questions answered.** "How cost-efficient?" "Vs. competitors — who and what cost?" "Energy per shot / coilgun / charging?" "USD/kg deorbited? FG01 weight? Cost/kg to orbit? Re per Al?" "USD wasted per missed shot?"

**Governing equations.** `C_gross = C_shot/[m_t·max(P_hit,10⁻⁶)]`, C_shot = Re material + launch + amortized platform; if Net removal ≤ 0, cost-per-kg-removed is undefined (AM-8). `DERIVED`

**Inputs.** P_Re = 2,600 and 7,283 (run both, dated); C_launch = 3,000–7,000 (both ends, labeled); P_elec = 0.10 USD/kWh; C_platform = 50M–500M amortized over 10 yr; P_hit from SIM-4 (CW primary); m_p,design from SIM-1 (≈3–4× smaller than v4 → material cost per shot drops). Competitor baselines (cite each): ground laser $100–500/kg; EDT $1,000–5,000/kg; aerogel $5,000–20,000/kg; ADR tug $10,000–100,000/kg.

**Outputs.**
1. Cost breakdown per shot.
2. C_rhenium per kg debris vs. h, across β and both launch endpoints.
3. Sensitivity: cost vs. P_hit, P_Re, C_launch.
4. Wasted USD per shot vs. P_hit.
5. Side-by-side table vs. laser / EDT / aerogel / ADR tug.
6. **Range conclusion:** "The rhenium method costs $X_low–X_high/kg (launch-profile + β range). Break-even vs. the cheapest viable competitor (EDT) requires $Z improvement — achievable only if P_hit ≥ P* and material switches to tungsten." Report the tungsten-substituted cost as a parallel column.

**Anti-mistake checks.** AM-1, AM-8, AM-9. **Validation.** Sensitivity analysis; arithmetic cross-check.

---

### SIM-10 — Regulatory, legal & aviation-safety compliance

**Purpose.** Check binding regulations against the full post-engagement population.

**Questions answered.** "Legal? If not, how to make it so?" "Okay to deorbit?" "Track debris so aircraft aren't hit?" "Environmental effect?"

**Checks.** FCC 47 CFR 25.283 (5-yr): does the entire launched mass (missed Re + any fragments) decay within 5 yr? NASA ODMSP (25-yr). IADC net-debris-flux. Liability Convention (1972) state liability. NEPA-style environmental for Re₂O₇. Aviation: NOTAM window + aircraft-encounter probability.

**Inputs.** Post-engagement population from SIM-5; lifetimes from SIM-2; reentry survival from SIM-8; FAA corridor density.

**Outputs.**
1. Compliance matrix (PASS/FAIL) for FCC, NASA, IADC, Liability, Environmental — **per altitude**, since SIM-2 makes compliance altitude-dependent.
2. Required mitigations where FAIL.
3. NOTAM window (time + geographic bounds) per deorbit event.
4. Aircraft-encounter probability per event.
5. Conclusion: "Compliant below h ≈ X km; non-compliant above, driven by missed-grain persistence." Note: v5 should PASS FCC/NASA at low altitudes where v4 failed — quantify the improvement honestly, including where it still fails.

**Anti-mistake checks.** AM-11, AM-17. **Validation.** Regulatory-text cross-check.

---

### SIM-11 — FG01 platform design, mass budget & end-of-life

**Purpose.** Analyze the transporter spacecraft.

**Questions answered.** "FG01 weight?" → design parameter, optimize. "Deorbit itself or return to hub to reload?" "Does the transporter burn up with missed rhenium?" → No; separate objects; compute separately.

**Governing equations.** Subsystem mass-fraction estimation; Tsiolkovsky for propellant. `DERIVED`

**Mass-budget template.** Structure (20% dry), avionics (50 kg COTS), power (solar+battery, 5 kW @ 30 W/kg), propulsion (chemical I_sp = 300 s, or electric 3,000 s as sensitivity), coilgun (100–500 kg, SIM-6), rhenium payload (10–100 kg design variable), thermal (30 kg). Dry = Σ; wet = dry + m_p.

**Inputs.** h_hub = 400 km; h_operational = 400–1000 km; I_sp = 300 s (or 3,000 s labeled sensitivity); mission life 10 yr; shots per reload from payload/m_p,design (v5 uses less Re per shot → more shots per reload — quantify).

**Outputs.**
1. Mass budget with parametric ranges.
2. Δv and propellant for hub rendezvous vs. direct deorbit.
3. Recommendation: "Deorbit FG01 if propellant < X kg; else rendezvous."
4. Reentry analysis for FG01 itself (separate from missed rhenium).

**Anti-mistake checks.** AM-9. **Validation.** Cross-check vs. known spacecraft (Dragon, Cygnus).

---

### SIM-12 — Global uncertainty quantification & sensitivity

**Purpose.** Honest joint uncertainty propagation across the whole pipeline; produce the headline P(viable).

**Questions answered.** "Does the error matter?" → quantify. "How much does error cut efficiency?" → posterior of η_eff. "Prove viability" → P(all constraints satisfied simultaneously).

**Method.** Monte Carlo, 10⁶ samples, fixed seed; report convergence + Sobol first-order and total indices.

**Input distributions.** β ~ U(1.0, 2.5); v_rel ~ U(500, 2500) within SIM-0 safe zone; P_hit ~ N(μ,σ) from SIM-4; m_t ~ N(1.414, 0.1) g; P_Re ~ U(2600, 7283); C_launch ~ U(3000, 7000); E_s,c ~ U(30, 50); τ ~ (discrete or log-uniform over 0.1–100 ms); Δh ~ U(150, 250); h swept.

**Viability definition (5 joint constraints).** (1) E_s < E_s,c (non-catastrophic); (2) P_hit > 0.5 [v5 raises the bar from v4's 0.01 — justify]; (3) T_decay(shifted) < 25 yr (stretch: report < 5 yr separately); (4) net debris flux ≥ 0; (5) cost < $100,000/kg removed.

**Outputs.**
1. Posterior histograms: design Δv, cost/kg, T_decay.
2. **P(viable) as a function of altitude** — this is how the plan finds the ~100% envelope: report the altitude band (and τ ceiling) where P(viable) is maximized, and the value it reaches there.
3. Sobol indices — expect τ (control lag) and altitude to dominate; confirm/refute.
4. Conclusion: "P(viable) peaks at X% for h ∈ [A, B] km with τ ≤ C ms; the dominant residual uncertainty is Y." State plainly the conditions required for X → ~100% and which of them are UNVALIDATED.

**Anti-mistake checks.** AM-12, AM-18. **Validation.** Convergence tests; compare to analytical limits.

---

### SIM-13 — Operating-point selection & viability maximization (new in v5)

**Purpose.** Turn the pipeline into an explicit answer to the PI's goal: find the single operating point (h, Δh, v_rel, grain mode, τ requirement, material) that maximizes P(viable), and report the maximum honestly.

**Questions answered.** "What is the most optimal altitude?" "What is the best achievable viability, and under what conditions?"

**Method.** Over the SIM-12 sample grid, locate argmax P(viable). Report the constraint that binds first as you move away from the optimum in each direction (the "why not higher/lower" answer). Produce the design card.

**Outputs.**
1. **Design card:** optimal h*, Δh*, v_rel*, grain mode, required τ*, sensor, material, f(h*), cost/kg, P(viable at optimum).
2. Robustness: how P(viable) degrades as each input drifts from optimum.
3. **Binding-constraint ledger:** at the optimum, list which of the five viability constraints has the least margin — that is the concept's true weak point and the next thing to engineer.
4. Explicit statement of the conditions the maximum viability rests on, each tagged DERIVED / SOURCED / UNVALIDATED.

**Anti-mistake checks.** AM-17, AM-18. **Validation.** Confirm the optimum is inside the demonstrated EM-launch envelope (SIM-6) and the achievable-τ envelope (SIM-4), not an extrapolation.

---

## PART 5 — How this plan maximizes viability (honestly)

The PI's target is P(viable) → ~100%. The plan pursues that through four legitimate levers, each tied to a SIM, rather than through optimistic assumptions:

1. **The Δv cut (SIM-0/1).** The 200 km shift needs ≈¼ of the direct-reentry Δv, which drops design-point E_s from ~80 kJ/kg to ~22 kJ/kg — below the catastrophic gate at every altitude. This lever is `DERIVED` and robust; it is the reason the concept is worth re-analyzing at all.
2. **Sub-catastrophic swarm delivery (SIM-3/5, AM-16).** Splitting the momentum across many small grains keeps *per-impact* energy low while delivering the *bulk* momentum. If validated, this makes fragmentation (v4 failure #3) structurally impossible. It is currently `UNVALIDATED` — treat it as a hypothesis the plan tests, not a fact it assumes.
3. **Cloud cross-section vs. control-lag (SIM-3/4).** A deliberately larger cloud raises P_hit and relaxes the τ requirement — the direct counter to v4 failure #2. The plan finds how much τ-tolerance a given σ_cloud buys and what momentum-density it costs.
4. **Altitude/operating-point selection (SIM-2/12/13).** Restricting operations to the altitude band where the decay gate passes and P_hit is achievable is where P(viable) is actually maximized.

**The honest ceiling.** Even with all four levers, three constraints can hold P(viable) below 100% and must be reported at the optimum, not buried: the terminal-guidance τ requirement (does the optimum need τ < 1 ms?), the decay-time ceiling (does the densest debris band at 800–1000 km fall outside the compliant zone?), and the swarm-equivalence validity (does AM-16 survive contact with swarm-impact literature?). If the optimum requires τ the field cannot yet deliver, the honest headline is "viable **conditional on** sub-millisecond terminal guidance," not "viable."

---

## PART 6 — Methodological risks that darken viability or could discredit the method

**This section is mandatory reading and mandatory reporting.** These are the issues a hostile reviewer (or the project's own v4 audit) will go straight to. Surfacing them up front is what makes a high viability number credible; hiding them is what invalidated earlier versions.

**R1 — Terminal guidance is the load-bearing weakness and the pivot does not fix it.** v4's P_hit collapsed by ~50 orders of magnitude between 1 ms and 10 ms control lag, and control lag drove 68% of variance. Lowering Δv changes nothing here. If SIM-4's optimum needs τ < 1 ms, the concept is guidance-limited regardless of the fragmentation result. *Do not let the clean SIM-0 result mask this.* Report the required τ next to every viability figure.

**R2 — The decay-time gate may exclude exactly the debris you care about.** A 200 km shift is a fixed decrement; its lifetime benefit collapses with altitude. The PI notes debris is densest at 600–1000 km — precisely where a shift to 400–800 km can still leave multi-decade lifetimes that fail the 25-yr (let alone 5-yr) rule. There is a real risk the "optimal altitude" for viability sits *below* the altitude where the debris actually is. Report P(viable) and debris-density on the same altitude axis so the mismatch is visible.

**R3 — Rhenium is hard to justify against tungsten, and the prior art already chose tungsten.** Momentum is density-independent (p = mv), so rhenium's founding rationale is physically void. Tungsten is ~50–100× cheaper, ~1,000× more available (60–81 t/yr Re is a hard supply ceiling for any at-scale campaign), and was the material the NRL dust-cloud concept selected. Rhenium's *only* plausible edge — smaller/denser missed grains decaying faster — must be shown to survive SIM-2 or the material choice reads as motivated. A reviewer will raise this immediately; pre-empt it in SIM-7.

**R4 — The swarm-vs-single-impactor equivalence (AM-15/16) is doing quiet heavy lifting.** Two central v5 escapes — sub-catastrophic per-impact energy and a large interception cross-section — both depend on converting a solid slug into a well-behaved grain swarm and on N small impacts summing to the intended bulk momentum without anomalous coupling, spallation, or ricochet. There is no experimental basis for this at cm-target, sub-orbital-velocity scales. If this assumption fails, v5's fragmentation and guidance escapes fail with it. Flag it in every dependent output.

**R5 — β is unvalidated and the whole momentum budget rides on it.** No confirmed Al-on-Al, cm-scale, sub-orbital-velocity momentum-enhancement measurement exists. The flat U(1.0, 2.5) envelope is defensible but wide; a low real β inflates required mass, cost, and per-impact energy simultaneously. Never present a β-dependent result as a point value.

**R6 — Missed rhenium is a regulatory and debris liability of its own.** Every grain that misses stays on a near-original orbit. At high altitude these persist well beyond the FCC 5-yr rule; the method can be a *net* debris source through its own misses even when the kick itself is clean. Missed-mass persistence must be in the compliance matrix and the net-debris tally, not treated as negligible.

**R7 — "Accelerated reentry" defers rather than removes risk.** Debris on a lowered but still-eccentric orbit remains a collision hazard for months–years before it reenters, passing repeatedly through populated shells. IADC net-flux accounting should credit the *eventual* removal but debit the *interim* on-orbit time. A method that lowers an object into a busier shell for a decade is not unambiguously beneficial.

**R8 — Cost is unlikely to beat the cheap competitors even after the pivot.** Lower Δv reduces Re mass per shot, but launch, platform amortization, and misses dominate C_shot. Ground-laser ($100–500/kg) and EDT ($1,000–5,000/kg) set a low bar. Report cost as a range against these and state candidly where FG01 lands; a favorable cost claim that ignores misses (AM-8) is the exact error that has to be avoided.

**R9 — Provenance discipline.** v4's audit found one likely-fabricated citation. Every SOURCED number needs a real, dated DOI/report/page. The 40 J/g threshold (Johnson 2001), NRLMSISE-00 version, coilgun benchmarks (Turman 1996, McNab 2003, ONR 2008), and every price must be independently verifiable. One bad citation discredits the whole document.

**R10 — Detection of the target in the first place.** The concept assumes the 1–10 cm object can be detected and terminally tracked at engagement range. Sub-10 cm objects are below the reliable ground-catalog threshold; FG01's own terminal sensor budget (SIM-4) must close this, and if it cannot, the "no rendezvous needed" advantage over conventional ADR partly evaporates. State the detection assumption explicitly.

---

## PART 7 — Deliverables & reporting format

1. **Master Results Document** — SIM-0 … SIM-13, every table/figure reproducible from code.
2. **RQ-indexed Results Document** — organized by the FG01 techs.md question order (viability → mode → slow-down → altitude → angle → mass ratio → precision → junk creation → cost → legality/environment), not by SIM number. This is what the manuscript's Results section is written from.
3. **Code repository** — versioned (Git); `python run_all.py` regenerates every result; pinned `requirements.txt`; fixed Monte Carlo seed; stated NRLMSISE-00 implementation/version; archived with a DOI (Zenodo).
4. **Citation log** — one file listing every external number (β, prices, launch costs, thresholds, material properties) with DOI/report/page and date. This is the artifact that lets a reviewer check R3/R5/R9 in one pass.
5. **Consolidated "Limitations & Unvalidated Assumptions" section** — one place, listing every UNVALIDATED item from Part 6 (β, AM-16 swarm equivalence, δv_launch conversion, SBM mass-term choice, τ achievability). Do not scatter caveats.

**Every number** carries (a) a range/CI from SIM-12 wherever it depends on β, P_hit, price, launch cost, τ, or Δh, and (b) a DERIVED / SOURCED / UNVALIDATED tag. **Figures:** vector (SVG/PDF), one figure = one finding, uncertainty shown as shaded bands, colorblind-safe, self-contained captions. **Tables:** units in headers, ranges not bare points, formatted + CSV.

---

## PART 8 — Question → SIM map (from FG01 techs.md)

| Research aim (FG01 techs.md) | Primary SIM | Output to look for |
|---|---|---|
| Prove momentum-kick cloud viable | SIM-0 + SIM-12 + SIM-13 | Viability boundary + P(viable) + optimal operating point |
| Dust cloud vs. grains; how small | SIM-3 | Mode comparison + grain-size recommendation |
| Slow debris enough to reenter | SIM-0 (Δv) + SIM-2 (decay) | Safe-zone Δv meets shift; T_decay < 25 yr |
| How high an orbit can we target | SIM-2 | Altitude ceiling where shifted-orbit decay complies |
| Launch angle; path accuracy | SIM-4 | Retrograde θ; d_miss; P_hit |
| Estimate Re required, f(h) | SIM-1 | f(h) band vs. altitude |
| Velocity/mass ratio; distance-dependence | SIM-1 | Ratio depends on h & β, not R |
| Shot precision; does error matter | SIM-4 + SIM-12 | P_hit, η_eff, variance contribution |
| Will it create more junk | SIM-5 | Net debris flux (v5 expected ≥ 0 sub-catastrophic) |
| Cost efficiency; vs. competitors | SIM-9 | USD/kg range + side-by-side table |
| Energy per shot; coilgun dimensions; charging | SIM-6 | E_elec, L_b, N_coils, AWG, m_Cu |
| USD/kg deorbited; FG01 weight; cost/kg to orbit; Re per Al | SIM-9 + SIM-11 + SIM-1 | Direct outputs |
| USD wasted per missed shot | SIM-9 | C_waste vs. P_hit |
| Okay to deorbit; legal | SIM-10 | Compliance matrix per altitude |
| Environmental effect | SIM-8 + SIM-10 | Re₂O₇ deposition + NEPA-style check |
| Track trash so aircraft aren't hit | SIM-10 + SIM-11 | NOTAM window + encounter probability |
| Trash fully burns; rhenium fully burns | SIM-8 | f_ablate for Al and Re |
| More junk from transporter/missed Re; deorbit/reload FG01 | SIM-11 + SIM-2 | Mass budget, EOL, missed-grain decay |

---

## Appendix A — Change log v4 → v5

- **Reframed the deorbit mechanism** from direct reentry (perigee → 100 km) to accelerated reentry (≈200 km shift + natural decay). Δv, and therefore design-point E_s, drop by ≈4×, lifting v4's structural margin–disruption falsification (SIM-0, 1.4).
- **SIM-2 promoted to a primary gate** (AM-17): success is now defined by post-kick orbital lifetime clearing the regulatory bar, which is altitude-dependent and not guaranteed by a 200 km shift.
- **Added SIM-13** (operating-point selection) to answer the "most optimal altitude / maximum viability" question directly.
- **Added AM-16, AM-17, AM-18** governing per-impact vs. bulk energy, the decay-time gate, and the anti-cherry-pick rule.
- **Raised the P_hit viability bar** from 0.01 (v4) to 0.5, with the guidance question re-centered as the load-bearing risk (R1).
- **Updated rhenium pricing** to 2026 spot (≈$7,283/kg, 2026-07-22) with the 2025→2026 trajectory noted.
- **Confirmed the mechanism as a targeted kick** (not a global dust layer) with the PI, and documented the resulting targeted-vs-indiscriminate trade against the NRL prior art (1.5).

## Appendix B — Reference documents in this folder

- `FG01_KIDR_Manuscript.md` — the v4 direct-reentry falsification study (history; the four failure modes it found are in 1.3).
- `FG01 techs.md` — the research-aims list (mapped in Part 8).
- `prompt 2.md` — the v4 SIM execution prompt (this plan's structural ancestor).

