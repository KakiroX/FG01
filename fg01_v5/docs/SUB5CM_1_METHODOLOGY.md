# FG01 / KIDR — Methodology
## Momentum-transfer remediation of sub-5 cm LEO debris

**Document 1 of 3** · [Results](SUB5CM_2_RESULTS.md) · [Viability](SUB5CM_3_VIABILITY.md)

---

## 1. Scope and definitions

This study assesses whether a targeted hypervelocity grain cloud, launched from a
co-orbital platform, can remediate **aluminium debris fragments of characteristic
diameter ≤ 5 cm** in low Earth orbit (400–1200 km).

| Term | Definition used throughout |
|---|---|
| **Target class** | Al-6061 solid spheres, d = 0.5, 1, 2, 5 cm (m = 0.18, 1.41, 11.3, 177 g) |
| **Remediation** | Reducing post-engagement orbital lifetime below the 25-yr NASA ODMSP bar |
| **Accelerated reentry** | A single retrograde impulse lowering perigee ≈200 km; atmospheric drag completes the deorbit |
| **Engagement** | One platform–target encounter, comprising one or more shots |
| **Success** | Five joint constraints simultaneously satisfied (§7.3) |

**Debris mass model (AM-2).** Targets are modelled as *solid* Al-6061 spheres,
ρ = 2700 kg/m³. This is stated rather than quietly assumed because it sets both
target mass and silhouette. Real fragments are irregular with higher and
tumbling-dependent area-to-mass ratio, so real lifetimes are **shorter** than
modelled — the decay gate is therefore conservative, but the band in which the
kick has remediation value is correspondingly **narrower** (§9, L14).

**Out of scope.** Objects > 5 cm are excluded except where needed to locate the
boundary of the sub-5 cm regime. Target acquisition and cataloguing are not
designed (§9, L2).

---

## 2. Research-integrity protocol

### 2.1 Traceability

Every quantity carries one of four tags, applied consistently in code, outputs
and text:

| Tag | Meaning |
|---|---|
| `DERIVED` | Follows from stated physics or definitions within this study |
| `SOURCED` | Specific external citation, dated, recorded in the citation log |
| `UNVALIDATED` | Assumption with no direct source at the relevant scale |
| `DESIGN` / `CHOSEN` | A variable the programme controls, swept or selected |

Provenance for every external number, with a **verification status** separating
what was independently checked from what was inherited, is maintained in
`CITATION_LOG.md`.

### 2.2 Anti-mistake contract

The study operates under a binding contract whose rules exist because each
encodes an error that invalidated an earlier version of this work. Compliance is
asserted per-module in the results document.

| Rule | Requirement | How discharged |
|---|---|---|
| AM-1 | No assumed interception efficiency | ε computed from a resolved areal-density field (§5.2) |
| AM-2 | Solid-sphere target mass, stated | §1 |
| AM-3 | Disruption gate evaluated **before** sizing | SIM-0 runs first in the pipeline |
| AM-4a | β never collapsed to a point value | Envelope carried; later **derived** (§5.3) |
| AM-5 | Hit reported as probability, never as certainty | Superseded by a delivered-mass formulation (§5.4) |
| AM-6 | Real NRLMSISE-00, not a fitted curve | §4.2, fidelity measured |
| AM-7 | Rate-limited reentry heating | §6.2 |
| AM-8 | Misses included in cost; undefined if net removal ≤ 0 | §7.2 |
| AM-9 | Demonstrated technology only | Launcher benchmarks dated and sourced |
| AM-10 | Report envelopes, not point estimates | All results as ranges |
| AM-11 | Downstream modules read upstream **outputs**, not re-typed constants | Enforced at runtime |
| AM-12 | Joint Monte Carlo, with convergence and Sobol indices | §7.3 |
| AM-13 | SBM effective-mass branch verified | §6.1, both branches implemented |
| AM-14 | P_hit consumes the **computed** cloud cross-section | §5.4 |
| AM-15/16 | Swarm-vs-single-impactor equivalence flagged; per-impact and bulk energy reported together | §5.2, discharged by explicit swarm modelling |
| AM-17 | Decay time is a hard gate, proven per altitude | §4.3 |
| AM-18 | Operating-point selection must not hide failures | Failure envelope reported in the same tables |

### 2.3 Falsification-first posture

Constraints are evaluated as **gates**: a configuration that fails any gate is
reported as failing, not re-parameterised until it passes. Where a model
correction *worsens* the result it is adopted and stated (three such corrections
are recorded in the results document).

---

## 3. Orbital mechanics

**Required impulse.** For a circular orbit at altitude *h*, the retrograde Δv at
apogee that lowers perigee by Δh (vis-viva): `DERIVED`

$$\Delta v(h)=\sqrt{\frac{\mu}{r}}-\sqrt{\mu\left(\frac{2}{r}-\frac{1}{a'}\right)},\qquad a'=r-\frac{\Delta h}{2},\quad r=R_E+h$$

**Self-disposal floor.** Because the projectile is fired retrograde, the
intercept velocity also determines the fate of everything that misses. Requiring

$$v_{rel}\ \ge\ \Delta v_{\text{direct}}(h)$$

places all released mass on a trajectory whose perigee is below the atmosphere.
This constraint is carried explicitly and is `DERIVED`.

**Constants.** μ = 3.986004418×10¹⁴ m³/s², R_E = 6.371×10⁶ m, J₂ = 1.08262668×10⁻³.

---

## 4. Atmospheric decay

### 4.1 Propagator

Orbit-averaged Gauss variational equations in (a, e). Drag is evaluated at **96
Gauss–Legendre nodes in eccentric anomaly**, with the full vector co-rotation
term **v_rel = v − ω⊕ × r**, time-weighted by (1 − e cos E). Integration is RK2
(midpoint) with an adaptive step limiting Δa to min(2 km, 0.2 % of a) and Δr_p to
2 km per step. `DERIVED`

A full 3-D Cartesian RK4 propagator (two-body + J2 + drag) is retained solely as
an independent validator.

### 4.2 Atmosphere (AM-6)

Density comes from the **real NRLMSISE-00 model** (`nrlmsise00` 0.1.2, Python
wrapper of Brodowski's C port), not a fitted curve. A 1 km-resolution table is
built by direct model evaluation averaged over:

- latitude, cos-weighted (equal-area), −82.5° to +82.5° in 15° steps
- longitude, 8 points
- local solar time, 4 points
- season, 4 points (equinoxes and solstices)

Averaging is appropriate because an object decaying over months to decades
samples all longitudes, local times and a wide latitude band. Solar drivers are
held at F10.7 = 70 / 150 / 250 with Ap = 4 / 15 / 45, run as three **bracketing
cases** rather than a single history (§9, L10). `SOURCED` model,
`UNVALIDATED` Ap pairing.

### 4.3 Decay gate (AM-17)

Reentry is declared at perigee altitude < 100 km. For every (size, altitude,
Δh, solar case) the pre- and post-kick lifetimes are computed and compared
against the 25-yr (ODMSP) and 5-yr (FCC 47 CFR 25.283) bars. The **altitude
ceiling** is the highest starting altitude at which a kick brings the object
inside the bar; it is proven per altitude, never assumed.

### 4.4 Remediation-value test

A gate the original plan did not contain. An object that would decay unaided
within 25 years has not been *remediated* by being made to decay sooner.
Engagements are credited only where

$$T_{\text{natural}}\ \ge\ 25\ \text{yr}\quad\wedge\quad T_{\text{post-kick}}<25\ \text{yr}$$

Without this test the optimiser selects operating points that pass every physical
gate while achieving nothing.

---

## 5. Cloud–debris interaction

### 5.1 Engagement geometry

Only the component of the impact relative velocity **opposing the target's
velocity** lowers the orbit. For a crossing engagement of plane angle θ the
retrograde efficiency is `DERIVED`

$$\cos\psi=\sin(\theta/2)$$

so useful-Δv per unit deposited specific energy is β/v_orbital for *any* crossing
angle, against 2β/|w| co-orbital. The co-orbital advantage is 2·v_orbital/|w|.
This is a kinematic result and it **forces** the co-orbital architecture before
any guidance argument is made.

The operational case — a nearly co-orbital platform with plane offset θ — has

$$E_s\ \propto\ 1+\left(\frac{v_{\text{geom}}}{v_{\text{launch}}}\right)^2,\qquad v_{\text{geom}}=2v_{\text{orb}}\sin(\theta/2)$$

which defines the **firing cone** that the engagement-rate model must respect.

### 5.2 Resolved cloud (AM-1, AM-15/16)

The cloud is modelled as a 2-D Gaussian areal mass-density field at the intercept
plane. The fraction of launched mass intercepted by a target disc of radius r_t
offset by d from the centroid is exact (non-central χ², equivalently 1 − Marcum Q₁):

$$\varepsilon=F\!\left(\left(\tfrac{r_t}{\sigma}\right)^{2};\,df=2,\,nc=\left(\tfrac{d}{\sigma}\right)^{2}\right)
\;\xrightarrow[\ r_t\ll\sigma\ ]{}\;\frac{A_t}{2\pi\sigma^{2}}$$

Cloud σ at intercept follows from the launcher's transverse dispersion,
σ = δv_T·R/v_rel. Clohessy–Wiltshire relative motion is evaluated exactly, and
shows that at engagement times-of-flight (n·t ≤ 2.8×10⁻⁴) it is numerically
indistinguishable from straight-line spreading — so the cloud size is a **design
variable**, not an orbital-mechanics outcome.

Grain count on target and the Poisson variability of delivered momentum are
computed explicitly, discharging AM-15/16 by modelling the swarm rather than
assuming equivalence. Both **bulk** and **per-impact** specific energy are
reported together in every dependent output (AM-16).

### 5.3 Momentum coupling β (AM-4a)

No measurement of β exists for Al-on-Al at centimetre scale and sub-km/s. Rather
than carry an assumption, β is **derived** from crater-ejecta scaling
(Housen–Holsapple, strength regime):

$$\beta=1+\frac{p_{\text{ejecta}}}{p_{\text{impactor}}},\qquad
\frac{p_e}{p_i}=\frac{3\mu k}{1-3\mu}\left[1-\left(\frac{v_{\min}}{v_{imp}}\right)^{1-3\mu}\right]$$

with the ejection-speed floor set by the target's own strength velocity
v_min = √(Y/ρ) = 320 m/s for Al-6061. An oblique-incidence correction
⟨cos i⟩ = 2/3 is applied for a curved target.

**Consistency check:** the same model at 6 km/s reproduces the hypervelocity
β ≈ 2–3 reported in the literature, which is why its much lower value in the
design regime should be believed rather than the literature value extrapolated
downward. β moves from `UNVALIDATED` to `DERIVED` — it is not `SOURCED`.

### 5.4 Delivery, not hit probability (AM-5, AM-14)

The success variable is **delivered momentum**, which is continuous in the miss
distance, not a Bernoulli hit. Launched mass to deliver a required incident mass
at confidence *c* against a Rayleigh miss of per-axis σ_miss:

$$\frac{M_{\text{launch}}}{m_{\text{incident}}}=\frac{-2\pi e\,\sigma_{\text{miss}}^{2}\ln(1-c)}{A_t},
\qquad \sigma_{\text{cloud}}^{*}=\sigma_{\text{miss}}\sqrt{-\ln(1-c)}$$

The legacy Bernoulli formulation is retained and reported **alongside** the
delivered-mass result so the effect of the change is visible rather than assumed
(AM-18).

### 5.5 Terminal error budget

Three miss models are carried side by side:

| Model | Lag term | Purpose |
|---|---|---|
| M1 | v_rel·τ | Reproduces the legacy formulation exactly |
| M2 | u·τ | Co-orbital sweep rate, no state extrapolation |
| M3 | via estimation error only | Design model, Keplerian extrapolation |

M3 sums position knowledge, velocity-extrapolation error, unmodelled relative
acceleration (3n²·sep), launcher pointing, muzzle-speed dispersion and ranging
timing in quadrature. The choice between M1 and M3 is the single largest
modelling decision in the study and is reported as such.

---

## 6. Secondary effects

### 6.1 Fragmentation (AM-3, AM-13)

NASA Standard Breakup Model, `N(>L_c) = 0.1·M_eff^0.75·L_c^−1.71`, with the
regime-dependent effective mass implemented in **both** branches and the branch
used reported with every result:

- catastrophic (E_s ≥ E_s,c): M_eff = m_t + m_p
- non-catastrophic: M_eff = m_p·v[km/s]

The SBM is calibrated to hypervelocity data; the design point is 0.5–1.1 km/s, so
its use here is an **extrapolation** and is labelled as such. An energy-based
cratering bound (crater volume = KE/Y) is computed in parallel as the physically
applicable estimate, and the largest liberated fragment is bounded by the crater
diameter.

### 6.2 Reentry ablation (AM-7)

Rate-limited, not a total-energy budget. Sutton–Graves stagnation-point
convective heating `q = k√(ρ∞/R_n)·v∞³` with radiative re-emission over the whole
surface, integrated along the trajectory. Mass loss begins at the **melting**
point — metals ablate by melt-stripping, which is what makes the threshold
consistent with tabulated effective heats of ablation.

### 6.3 Regulatory accounting

A distinction the source plan did not draw: FCC 47 CFR 25.283 and NASA ODMSP bind
*the operator's own* objects — the platform and everything it releases. The
target debris is pre-existing and not the operator's object, so its post-kick
lifetime is a **mission-effectiveness benchmark**, not a licence condition. Both
readings are reported.

---

## 7. Economics and uncertainty

### 7.1 Engagement rate

Derived rather than assumed. Fragments become engageable only when the platform
lies within the INT-0 firing cone of their plane *and* within a reachable
altitude band. Cone entries are enumerated analytically from differential J2
nodal regression; each engagement is then charged the Hohmann Δv to reach the
target's altitude.

Synthetic clouds are generated from real breakup parents using the SBM ejection-
velocity distribution and aged under differential J2 to the present epoch.

### 7.2 Cost (AM-8)

$$C_{\text{object}}=\underbrace{n_{\text{shots}}\,m_{\text{shot}}(P_{\text{mat}}+C_{\text{launch}})}_{\text{consumables}}+\underbrace{\frac{C_{\text{platform}}}{N_{\text{engagements}}}}_{\text{amortisation}}$$

Misses enter through the shot multiplier *and* through the full areal-mass
penalty charged on every shot. If net removal ≤ 0 the cost is reported as
**undefined**, never as a large finite number.

Both a per-kilogram and a per-object criterion are carried. For centimetre debris
these differ by orders of magnitude and neither is privileged: a 1 cm sphere
masses 1.414 g, so a $100,000/kg bar is $141 per object.

### 7.3 Uncertainty quantification (AM-12)

Monte Carlo over **joint** distributions, 10⁶ samples, fixed seed 20260723, with
convergence testing and Saltelli-estimator Sobol first-order and total indices.

**Viability = five joint constraints:**

| # | Constraint |
|---|---|
| 1 | Sub-catastrophic specific energy, E_s < E_s,c |
| 2 | Launched mass per shot within the magazine limit |
| 3 | Post-kick lifetime < 25 yr |
| 4 | Net debris ≥ 0 |
| 5 | Cost below the stated bar |
| (6) | Remediation value (§4.4) — applied as a sixth gate |

**Two sampling modes, reported separately.** In the global UQ, design variables
are *sampled* (how likely is a randomly configured system to work). In
operating-point selection they are *chosen* and only epistemic uncertainties are
sampled (how well does a deliberately designed system work). This decomposition
is stated because it is what makes a selected-optimum probability meaningful.

---

## 8. Numerical methods and validation protocol

| Item | Method |
|---|---|
| Decay propagation | Orbit-averaged Gauss VE, RK2 adaptive, 96-node quadrature |
| Validator | 3-D Cartesian RK4, J2 optional, fixed step |
| Interaction ODEs | RK4 |
| Root finding | Bisection to 10⁻⁴ relative |
| Monte Carlo | 10⁶ samples, fixed seed, convergence reported |
| Sensitivity | Sobol S₁ and S_T, Saltelli estimator |

**Validation is required at four levels**, and every check is reported with its
numerical error in the results document:

1. **Analytic identity** — closed forms against exact algebra
2. **Independent numerical method** — orbit-averaged against full 3-D Cartesian
3. **External observation** — decay model against an observed reentry
4. **Cross-regime consistency** — a derived model must reproduce the regime where
   published data exists

---

## 9. Consolidated limitations

Recorded here rather than scattered. Ordered by consequence.

| # | Limitation | Effect if wrong |
|---|---|---|
| **L1** | Engagement rate rests on a modelled cluster population, not a catalogue | Dominant economic term |
| **L2** | Target acquisition and terminal tracking assumed, not designed | Removes the "no rendezvous" advantage |
| **L3** | Launcher pointing (200 µrad) and muzzle repeatability (0.1 %) are engineering estimates | Launched mass scales as σ² |
| **L5** | β derived, not measured | Enters mass inversely |
| **L6** | Swarm arrival is neither cleanly simultaneous nor independent (~10² µs ≈ 50 shock transits) | Both bounds reported |
| **L7** | SBM extrapolated below its calibration velocity | Over-predicts fragments (conservative) |
| **L8** | Slug-to-swarm conversion mechanism not designed | No space-qualified mechanism exists |
| **L9** | Platform cost parametric ($50–500 M) | Propagates directly to cost/object |
| **L10** | Fixed F10.7 per case, not a cycle | Brackets rather than predicts |
| **L14** | Solid-sphere debris (AM-2) | Real lifetimes shorter; useful band narrower |
| **L16** | Lumped single-node reentry thermal model | Sensitivity run over H_eff and emissivity |
| **L17** | Interim on-orbit hazard argued from lifetime ratio, not integrated | Direction robust, magnitude not |

---

## 10. Reproducibility

| Item | Value |
|---|---|
| Entry point | `python run_all.py` |
| Random seed | `20260723` (fixed) |
| Python | 3.11.9 |
| Dependencies | `numpy 2.0.2`, `scipy 1.15.1`, `pandas 2.2.3`, `matplotlib 3.10.0`, `nrlmsise00 0.1.2` (pinned) |
| Outputs | CSV + Markdown tables, JSON conclusions and validation records |
| Figures | SVG/PDF vector, colourblind-safe palette |

Every table and figure in the results document regenerates from source. Downstream
modules read upstream results from `outputs/` at runtime rather than from re-typed
constants (AM-11), so an upstream correction propagates automatically.

---

*Document 1 of 3. Results are presented in [SUB5CM_2_RESULTS.md](SUB5CM_2_RESULTS.md);
the viability assessment and conclusions in [SUB5CM_3_VIABILITY.md](SUB5CM_3_VIABILITY.md).*
