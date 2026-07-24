# FG01 / KIDR v5 — Limitations and Unvalidated Assumptions

One consolidated list, as required by the plan (Part 7.5). Caveats are **not**
scattered through the results; every one of them is here, with what it would do
to the conclusions if it turned out to be wrong.

Ordered by how much damage each would do.

---

## Tier 1 — would change the headline conclusion

### L1. Engagements per platform lifetime is not derived
The dominant Sobol term (S_T = 0.52) and the binding constraint at the optimum.
P(viable) is **0 below ~3,000 engagements and 0.997 above 10,000**, so the entire
economic conclusion rests on a quantity this study samples over a log-uniform
prior (10³–10⁵) rather than deriving. Setting it requires a target-acquisition
and approach-cadence model — how often FG01 can detect, match orbits with, and
close to 10 m of a tracked centimetre object in its own shell. **That model does
not exist here.** `UNVALIDATED`

### L2. Target detection and terminal tracking is assumed, not designed (R10)
The concept assumes each 1–10 cm object can be detected and terminally tracked at
engagement range. Sub-10 cm objects are below the reliable ground-catalogue
threshold. The error budget assumes a 10 s tracking arc at 1 kHz with a
diffraction-limited 10 cm aperture, which is internally consistent at 10–100 m
range — but *acquiring* the target in the first place is outside the modelled
scope. If FG01 needs ground cueing, the "no rendezvous required" advantage over
conventional ADR partly evaporates. `UNVALIDATED`

### L3. Launcher pointing (200 µrad) and muzzle-speed repeatability (0.1%)
These two engineering estimates set the 5.775 mm miss distance, and launched mass
scales as the **square** of that number. They are plausible for spacecraft-grade
hardware but are not sourced to a built coilgun. A 1 mrad pointing error would
raise the mass penalty from 16.7× to ~400×, moving the shot from 2.3 g to ~54 g —
still affordable, but a 5 mrad error would not be. `UNVALIDATED` (E in the
citation log)

### L4. The co-orbital engagement geometry is a requirement, not a free choice
The whole guidance result depends on the platform–target relative speed u being
~10 m/s, i.e. FG01 operating **inside the target's orbital plane**. This is a
derived consequence of the energy gate (a low-v_rel intercept from a fast flyby is
geometrically impossible), but it means one platform serves one inclination band —
plane changes cost 130 m/s per degree of inclination. Global coverage needs a
fleet. If a future variant needs to engage across planes, the v4 sweep-rate model
becomes the correct one again and the guidance wall returns. `DERIVED` constraint,
`UNVALIDATED` operational feasibility

---

## Tier 2 — would change magnitudes, not direction

### L5. β is unvalidated and the momentum budget rides on it (R5, AM-4a)
No confirmed Al-on-Al, cm-scale, sub-orbital-velocity momentum-enhancement
measurement exists. A flat U(1.0, 2.5) envelope is used and every value is run;
no result is reported as a point estimate in β. **A physical caution the plan does
not state:** at 619 m/s the impact is well below the hypervelocity regime where
ejecta-driven momentum multiplication arises, so the true β is likely near the
*bottom* of the envelope (≈1.0–1.2). Using the envelope mean flatters required
mass by ~1.75×. At β = 1.0 the design point is still sub-catastrophic, but the
margin narrows and the velocity ceiling drops to ~756 m/s (E_s,c = 40) or
~567 m/s (E_s,c = 30).

### L6. Swarm-versus-single-impactor equivalence (R4, AM-15/AM-16)
Treating N grain impacts as one impactor of equal total mass has no experimental
basis at these scales. **The design deliberately does not rely on it:** the bulk
(single-impactor) specific energy at the design point is 21.8 kJ/kg, already
sub-catastrophic, so the per-impact escape hatch (136 J/kg, 0.34% of threshold) is
a robustness reserve rather than a load-bearing assumption. It becomes
load-bearing only at β = 1.0 with E_s,c = 30 kJ/kg, where the bulk figure would
otherwise fail. A further unmodelled subtlety: grains arrive over ~10² µs, which
is ~50 shock-transit times across a 1 cm target, so impacts are neither cleanly
independent nor cleanly simultaneous. Both bounds are reported; neither is
claimed to be correct.

### L7. The NASA SBM is extrapolated below its calibration
The SBM was fitted to hypervelocity (≳1–2 km/s) data; the v5 design point is
0.5–1.1 km/s. Fragment yields are therefore over-predicted — the conservative
direction for a debris-generation claim, but it means the reported 0.26 fragments
>1 cm is a model artefact rather than a prediction. The energy-based cratering
bound (0 fragments >1 cm, 0.167 g ejecta) is the physically applicable estimate,
and it rests on a cratering resistance Y = 5×10⁸ J/m³ that is itself an
order-of-magnitude estimate. `SOURCED but extrapolated` / `UNVALIDATED`

### L8. δv_launch → grain swarm conversion
No basis exists for converting a solid slug into a controlled Gaussian shot
pattern of the specified σ at these scales. The required transverse dispersion
(σ = δv_T·R/v_rel → ~0.4 m/s at 10 m and 619 m/s, i.e. ~0.6 mrad) is modest and
has a direct terrestrial analogue in shotgun ballistics, but no space-qualified
mechanism is designed here. `UNVALIDATED`

### L9. Platform cost and subsystem mass fractions
$50M–$500M is a parametric range, not an estimate from a costed design. The mass
budget uses conventional fractions (20% structure, 30 W/kg power) rather than a
vehicle study. Since cost per object is 99.4% platform amortisation, this range
propagates directly into the headline. `UNVALIDATED`

### L10. Fixed F10.7 per case
Solar drivers are held constant for the entire multi-decade decay, per the plan's
specification. Real activity cycles through all three cases on an 11-year period.
The three cases **bracket** the answer but no single one is a realistic history —
which is why no operating point is simultaneously remediation-relevant at a
permanent solar maximum and compliant at a permanent solar minimum, and why
selection is made on the nominal case. A cycle-averaged propagation would be the
correct refinement. `SOURCED model, UNVALIDATED usage`

### L11. Ap pairing with F10.7
Ap = 4/15/45 paired with F10.7 = 70/150/250 is an assumption, not a sourced
pairing. `UNVALIDATED`

---

## Tier 3 — bounded, quantified, or checked

### L12. Orbit-averaged propagation
Validated against full 3-D Cartesian RK4 to **0.01%** over 400 revolutions and to
1.9% on a full fast-decay case. J2 is omitted from the averaged model, which is
correct — J2 has no secular effect on a or e. Third-body and solar-radiation
pressure are omitted; both are negligible for these area-to-mass ratios in LEO.

### L13. Density tabulation
A 1 km log-interpolated table built from the real NRLMSISE-00 model, not a fitted
curve (AM-6). Interpolation error vs. direct calls: **max 1.9×10⁻⁴**. Global
averaging over lat/lon/LST/season is appropriate for multi-year decay but not for
a single-pass calculation.

### L14. Debris modelled as solid Al-6061 spheres (AM-2)
Stated, not quietly assumed. Real fragments are irregular with higher and
tumbling-dependent A/m, so real lifetimes are likely **shorter** than modelled —
which makes the decay gate conservative but also *narrows* the band where the kick
has remediation value.

### L15. Spatial-density profile for transient hazard
The catalogued-object density profile used in SIM-10 is an order-of-magnitude
construction, not a MASTER/ORDEM extraction. It supports the conclusion
(~10⁻⁶ per engagement) only at order-of-magnitude level. `UNVALIDATED`

### L16. Reentry thermal model
Lumped-mass, single-node, with q_stag applied over the frontal area as a stand-in
for the sphere-averaged load and an assumed emissivity. Sensitivity to H_eff
(3.5–5.5 MJ/kg) and emissivity (0.2–0.6) is tabulated. The qualitative conclusion
— design-range grains demise, and anything surviving is harmless — is robust
across the whole sensitivity grid.

### L17. Interim on-orbit hazard (R7) is estimated, not integrated
A kicked object spends months to years on a lowered, eccentric orbit crossing
populated shells before reentry. This study reports the lifetime reduction and the
released-mass transit risk, but does **not** compute the time-integrated collision
probability before versus after the kick. Since the kick cuts lifetime by 2.2–16.7×,
the integrated exposure almost certainly falls, but that is an argument, not a
calculation.

### L18. Cost metric is a policy choice
Per-kg and per-object criteria give qualitatively different verdicts for
centimetre debris, and neither is privileged here. Both are reported everywhere.
The plan's $100,000/kg bar is $141 for a 1 cm object, which no orbital system can
meet; this is stated as a property of the metric rather than being quietly
dropped.

---

## v6 additions — interaction model and economics

### L19. The Fengyun-1C ≥1 cm member count is a model output, not an observation
Only ~3% of the ≥1 cm population is tracked (40,000 of 1.2 million). Every
per-cluster member count in ECO-2 is modelled. **This no longer drives the
verdict** — ECO-1's binding constraint is Δv per engagement, which is purely
orbital-mechanical and independent of how many fragments exist — but any claim
about *opportunity* counts rests on it. `UNVALIDATED`

### L20. β is derived from a scaling law extrapolated below its calibration
Housen–Holsapple ejecta scaling is fitted to hypervelocity cratering; INT-2
applies it at 619 m/s, where the impact is only 1.9× the target's strength
velocity. The model reproduces β ≈ 2–3 at 6 km/s (the regime where data exists),
which is the consistency check, but the sub-km/s value is an extrapolation. A
light-gas-gun shot or hydrocode run at 0.6 km/s would close this. Note the
direction: the extrapolation *lowers* β relative to v5 and therefore darkens the
mass budget, so it is not a self-serving assumption. `DERIVED, not SOURCED`

### L21. The tour model charges round trips, not a chained tour
ECO-1 charges each engagement as a round trip from a base altitude. A
well-ordered tour chaining between neighbouring targets would pay at most the
one-way cost, so the true Δv is between 0.5× and 1× of what is charged — the
result is conservative by at most a factor of two. Against a shortfall of
~400×, this does not change the verdict. `DERIVED, conservative`

### L22. The value model spans three orders of magnitude
Avoided-collision value per object ranges $67–$20,124 across the swept
assumptions, and the upper end implies ~80 debris losses per year across the
active fleet against an observed record of essentially none. The model therefore
over-predicts and the low end is credible. The *gap* between cost and value is
robust (≥10²× everywhere in the sweep); the precise ratio is not. `UNVALIDATED`

### L23. Architecture B is costed, not designed
The large-object variant is competitive on cost per object, but 197–1,166 kg of
tungsten per target implies a platform far larger than the 400–600 kg vehicle
costed in ECO-3, and that mass has not been re-budgeted. The favourable cost
figure therefore carries an unmodelled platform-scaling penalty. Treat it as a
direction worth studying, not a result. `UNVALIDATED`

### L24. Cluster synthesis uses an isotropic SBM Δv
Real breakup fragment velocity distributions are anisotropic (they depend on
impact geometry) and the SBM's Δv law is itself fitted to a small number of
events. The RAAN-spreading conclusion is robust to this — any plausible Δv
spread produces full RAAN dispersal over 5–20 years — but the inclination-spread
figures (1.4–1.6°) inherit the assumption. `SOURCED model, UNVALIDATED usage`

### L25. Engagement cadence (8 h per engagement) is assumed
Not load-bearing at the realised counts, because the tour is Δv-limited before
it is cadence-limited. It would become load-bearing in any architecture that
solved the Δv problem. `UNVALIDATED`

---

## v7 addition — laser-ablation accuracy comparison

### L26. SIM-14's laser-side inputs come from a secondary brief, not primary literature
Every laser-specific number (orbit-prediction requirement, "sub-arcsecond"
pointing, the p-N-omission residual) is taken from `laser accuracy comm.md`,
an informal brief supplied for this comparison and not independently checked
against primary sources in this study (CITATION_LOG key `laser_brief`,
status E — same tier as this study's own engineering estimates, not the
tier of a verified citation). The brief itself gives ranges and qualitative
claims rather than point estimates for pointing, so this study's own
0.1–2.0 arcsec envelope and its 500 km engagement range are assumptions
layered on top of the brief, not sourced from it. `UNVALIDATED`

### L27. SIM-14 is scoped to accuracy tolerance only, and says nothing about cost
FG01 remains 3–4 orders of magnitude more expensive per object than
ground-based laser ablation (SIM-9); SIM-14 does not reopen, revisit, or
offset that finding. A reader should not infer overall superiority from a
single-dimension result. `SCOPE, stated in sim14_conclusions.json`

### L28. The laser is modelled as a fixed, binary-hit system; no defocus tradeoff
A real laser could in principle spread its spot to trade peak fluence for
aim tolerance, an analogue of FG01's own sigma-vs-mass tradeoff (SIM-3).
This is not modelled — the brief's own comparison table frames the laser's
failure mode as strictly binary ("miss entirely; no effect"), and SIM-14
follows that framing rather than inventing an untested laser design
freedom. If real systems do trade focus for tolerance, SIM-14's laser
numbers understate its achievable tolerance to an unquantified degree.
`UNVALIDATED`

---

## What would most change the answer

*(L1 and L2 below were written against the v5 results; v6's ECO-1 supplied the
engagement-rate model L1 called for, and the answer was ~24 per platform-decade
rather than the 10⁴ v5 had to assume. They are retained as written for the
record, with that resolution noted.)*

1. ~~**A costed engagement-rate model** (L1, L2)~~ — **DELIVERED in v6 (ECO-1).**
   It collapsed the economic case: the Δv-per-engagement invariant
   π·v/(3.5·Ω̇·T) gives ~107 m/s per engagement at Fengyun-1C, so a realistic
   platform achieves tens of engagements, not 10⁴.
2. **A single β measurement** at Al-on-Al, cm scale, 0.5–1 km/s (L5, L20) — still
   the widest physics uncertainty, now narrowed from U(1.0, 2.5) to 1.17–1.68 by
   derivation rather than measurement.
3. **A mass-budget study of Architecture B** (L23) — the only architecture that
   survived the economics, and the one least studied.
4. **A systemic cascade-value model** (L22) — per-object value is what condemns
   small-debris removal; cascade prevention is the only framing under which it
   could pay, and it is not modelled here.
5. **A demonstrated launcher pointing figure** (L3) — sets launched mass
   quadratically, though at realised engagement counts the consumable cost is
   negligible either way.
