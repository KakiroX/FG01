# FG01 / KIDR v5 — Citation Log

Every external number used anywhere in the study, with its source and a
**verification status recorded honestly**:

- **V** — independently verified during this study against a primary or
  authoritative secondary source, with the URL given.
- **P** — carried forward from the v5 plan / v4 manuscript parameter table and
  **not independently re-verified here**. These are the entries a reviewer should
  check first (R9).
- **E** — engineering estimate made in this study; no external source exists at
  the relevant scale. Treated as UNVALIDATED throughout.

R9 warned that the v4 audit found one likely-fabricated citation. This log
separates what was checked from what was inherited rather than presenting all
entries as equally solid.

---

## Impact and breakup physics

| Quantity | Value | Source | Status |
|---|---|---|---|
| Catastrophic-disruption threshold | **40 J/g** (= 4.0×10⁴ J/kg) | Johnson NL, Krisko PH, Liou J-C, Anz-Meador PD. *NASA's new breakup model of EVOLVE 4.0.* Advances in Space Research 2001;28(9):1377–1384. doi:10.1016/S0273-1177(01)00423-9 — [ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S0273117701004239) | **V** — threshold and its definition ("relative kinetic energy of the smaller object divided by the mass of the larger ≥ 40 J/g") confirmed |
| Non-catastrophic regime definition | "fragmentation of projectile and crater or hole on target" | same | **V** — directly supports the SIM-5 argument that this branch cannot yield >1 cm fragments from an undisrupted 1 cm target |
| SBM cumulative count law | N(L_c) = 0.1·M^0.75·L_c^−1.71 | same | **P** |
| Effective mass, non-catastrophic branch | M_eff = m_p·v[km/s] | Krisko PH, *Proper Implementation of the 1998 NASA Breakup Model*, Orbital Debris Quarterly News 15(4), 2011 | **P** — AM-13 flagged this as error-prone; both branches implemented explicitly and the branch used is reported with every result |
| Threshold uncertainty band | 30–50 kJ/kg | spread across the impact-test literature | **P** |
| Momentum-enhancement factor β | U(1.0, 2.5) flat | **no source exists** for Al-on-Al, cm-scale, sub-orbital-velocity impacts | **E / UNVALIDATED** (AM-4a) |
| Cratering resistance of Al-6061, Y | 5×10⁸ (3×10⁸–1×10⁹) J/m³ | order-of-magnitude from ductile-metal crater-volume scaling | **E** |
| Ejecta/projectile mass ratio | 10 (5–20) per km/s | Gault-type cratering scaling | **P/E** |

## Atmosphere and orbital decay

| Quantity | Value | Source | Status |
|---|---|---|---|
| Atmosphere model | NRLMSISE-00 | Picone JM, Hedin AE, Drob DP, Aikin AC. *NRLMSISE-00 empirical model of the atmosphere.* JGR Space Physics 2002;107(A12):1468. doi:10.1029/2002JA009430 | **P** (model identity); implementation verified in-session |
| Implementation | `nrlmsise00` 0.1.2 (Python wrapper of D. Brodowski's C port) | PyPI; version confirmed by `pip show` in this session | **V** |
| Table fidelity | max 1.9×10⁻⁴ relative error vs. direct model calls | computed in this study, `outputs/sim02_atmosphere_table_verification.json` | **V** (self-check) |
| F10.7 cases | 70 / 150 / 250 sfu | conventional solar-min/nominal/max values | **P** |
| Ap pairing | 4 / 15 / 45 | **assumed pairing**, not sourced | **E** |
| Drag coefficient, free-molecular sphere | 2.2 (2.0–2.4) | King-Hele DG, *Satellite Orbits in an Atmosphere*, 1987 | **P** |
| **Starshine 1** validation target | 48 cm diameter hollow Al sphere; 387 km, i = 51.6°; deployed from STS-96; reentered **18 Feb 2000** on its 4,212th orbit | [eoPortal](https://www.eoportal.org/satellite-missions/starshine); [Gunter's Space Page](https://space.skyrocket.de/doc_sdat/starshine-1.htm); [satobs.org decay report](http://www.satobs.org/starshinedecay.html) | **V** for diameter, orbit, inclination, reentry date |
| Starshine 1 mass | 39.5 kg (used); **39 kg confirmed** in this revision | eoPortal STARSHINE mission page: "the S/C mass is 39 kg" | **V** (confirmed against an independent mission-database source in this session, resolving the prior weakest-link status). The 1.3% discrepancy (39.5 vs. 39 kg) implies the model's A/m = 4.58×10⁻³ m²/kg understates the true value by ~1.3%, i.e., true decay is marginally *faster* than modelled — negligible against the model's own stated 9.8% validation error, and conservative in direction (it would not weaken the decay-gate conclusion) |
| Observed lifetime | 258 days (deployment 5 Jun 1999 → 18 Feb 2000); cross-check 4,212 orbits × 90 min ≈ 263 d | as above | **V** (consistent to 2%) |
| μ (Earth GM) | 3.986004418×10¹⁴ m³/s² | WGS-84 / EGM96 | **P** |
| J2 | 1.08262668×10⁻³ | EGM96 | **P** |

## Materials

| Quantity | Value | Source | Status |
|---|---|---|---|
| ρ Al-6061 | 2.70 g/cm³ | standard alloy data | **P** |
| ρ Re | 21.01 g/cm³ | USGS Mineral Commodity Summaries | **P** |
| ρ W | 19.25 g/cm³ | USGS Mineral Commodity Summaries | **P** |
| Re melting point | 3,459 K | CRC Handbook | **P** |
| W melting point | 3,695 K | CRC Handbook | **P** |
| Al melting point | 933 K | CRC Handbook | **P** |
| H_eff ablation, Re | 4.5 (3.5–5.5) MJ/kg | plan parameter table | **P** — sensitivity run over the full range (`sim08_thermal_sensitivity`) |
| H_eff ablation, Al | 12 (10–14) MJ/kg | plan parameter table | **P** |
| Surface emissivity | 0.2–0.6 | assumed | **E** — sensitivity run |

## Economics (all dated)

| Quantity | Value | Date/edition | Status |
|---|---|---|---|
| Rhenium spot | $7,283/kg | 2026-07-22 | **P — market quote inherited from the plan; not re-verified.** Also run at $6,389 (2026-Q1) and $2,486 (2025) |
| Rhenium production | 60–81 t/yr | USGS MCS | **P** |
| Tungsten price | $30–100/kg | USGS MCS | **P** |
| Tungsten production | ~84,000 t/yr | USGS MCS | **P** |
| Launch cost | $3,000 (dedicated) – $7,000 (rideshare marginal) /kg | dated commercial pricing | **P** |
| Starship-class launch | $100–500/kg | **labelled sensitivity only** (AM-9); used in no headline figure | **P** |
| Electricity | $0.10/kWh | — | **P** |
| Platform cost | $50M–$500M | parametric | **E** |
| Competitor: ground laser | $100–500/kg (theoretical) | plan §1.5 | **P** |
| Competitor: electrodynamic tether | $1,000–5,000/kg | plan §1.5 | **P** |
| Competitor: aerogel capture | $5,000–20,000/kg | plan §1.5 | **P** |
| Competitor: ADR tug | $10,000–100,000/kg | plan §1.5 | **P** |

## Electromagnetic launch

| Quantity | Value | Source | Status |
|---|---|---|---|
| Largest demonstrated EM launch | 33 MJ muzzle energy | US Navy / ONR EMRG, 2010 | **P** |
| 2008 record shot | 10.16 MJ, 3.2 kg @ ~2,520 m/s | US Navy / ONR | **P** |
| Coilgun efficiency | 0.30–0.50 | Turman BN (1996); McNab IR (2003) | **P** |
| Pulse-capacitor energy density | 1.5 J/g | current metallised-film pulse capacitors | **E** |

## Prior art and regulation

| Item | Detail | Source | Status |
|---|---|---|---|
| **Tungsten dust-cloud prior art** | Ganguli G et al., NRL. Micron-scale (10–30 µm) **tungsten** dust to enhance drag on sub-10 cm debris in sun-synchronous orbit; 20–40 t total; 12 t deorbits ballistic-coefficient ≤3 objects from 1100→900 km in 25 yr | *A Concept for Elimination of Small Orbital Debris*, [arXiv:1104.1401](https://arxiv.org/pdf/1104.1401) (2011); US Patent 8,469,314 B2 | **V** — material, particle size, tonnage, altitude band and 25-yr timescale all confirmed |
| FCC disposal rule | 5 yr | 47 CFR 25.283 | **P** |
| NASA ODMSP | 25 yr | ODMSP 2019 update | **P** |
| Casualty threshold | 1×10⁻⁴ per event | NASA-STD-8719.14 | **P** |
| Serious-injury KE threshold | 15 J | conventional in debris-casualty analysis | **P** |
| Mean human projected area | 0.36 m² | conventional | **P** |
| Liability Convention | absolute liability on the surface, fault-based in space | Convention on International Liability for Damage Caused by Space Objects, 1972, Art. II–III | **P** |
| Cosmic-dust influx | ~40,000 t/yr | conventional estimate | **P** |
| Catalogued-object spatial density | 1.5–3.5×10⁻⁸ objects/km³ across 400–1200 km | order-of-magnitude profile constructed for this study | **E** — used only for a transient-risk estimate, tagged UNVALIDATED in `sim10_transient_hazard` |

## Heating correlation

| Quantity | Value | Source | Status |
|---|---|---|---|
| Stagnation-point convective heating | q = k√(ρ∞/R_n)·v∞³, k = 1.7415×10⁻⁴ (SI) | Sutton K, Graves RA. *A general stagnation-point convective-heating equation for arbitrary gas mixtures.* NASA TR R-376, 1971 | **P** — note: the plan calls this "Fay–Riddell"; the k√(ρ/R_n)v³ form given is the Sutton–Graves correlation. The functional form and the rate-limited treatment required by AM-7 are unaffected |

---

## v6 additions — interaction model and debris population

| Quantity | Value | Source | Status |
|---|---|---|---|
| Ejecta-scaling exponent μ | 0.40–0.55 | Housen KR, Holsapple KA. *Ejecta from impact craters.* Icarus 2011;211(1):856–875. doi:10.1016/j.icarus.2010.09.017 | **P** — the scaling form is standard; its application below 1 km/s on a metal target is this study's extrapolation |
| Ejecta mass coefficient k | 0.2–0.5 (strength regime) | same | **P** |
| Al-6061 yield strength | 276 MPa | standard alloy data | **P** |
| Al-6061 bulk sound speed | 5,100 m/s | standard alloy data | **P** |
| Strength velocity √(Y/ρ) | 320 m/s | **DERIVED** from the two above | — |
| β at 619 m/s | 1.17–1.68, central 1.30, oblique-corrected **1.201** | **DERIVED** in INT-2 | **DERIVED, not SOURCED** — no direct measurement exists at Al-on-Al, cm-scale, sub-km/s |
| β hypervelocity cross-check | 1.96–8.27 at 6 km/s | model output vs. literature β ≈ 2–3 | **V** (internal consistency) |
| DART momentum enhancement | β ≈ 3.6 at 6.1 km/s | Daly RT et al., *Nature* 2023 (DART deflection) | **P** — cited only as an upper bound on a rubble-pile target, explicitly NOT applied here |
| Poncelet penetration coefficient A | 1.5 | standard form | **E** |
| Mean cos(incidence) on a sphere | 2/3 | **DERIVED** (projected-area weighting) | — |
| Ejecta size-distribution slope | 2.7 | assumed | **E** — bounded both ways in INT-5 |
| Largest ejecta fragment | 0.2× crater diameter (expected); 1.0× (conservative bound) | laboratory cratering convention | **E** — both bounds reported |

### Debris population (ECO-2)

| Quantity | Value | Source | Status |
|---|---|---|---|
| Total ≥1 cm population | **1.2 million** | ESA Space Environment Report 2025 / MASTER, Aug 2024 — [ESA](https://www.esa.int/Space_Safety/Space_Debris/ESA_Space_Environment_Report_2025) | **V** |
| Total ≥10 cm | 54,000 | same | **V** |
| Tracked objects | ~40,000 (~11,000 active payloads) | same | **V** |
| Fengyun-1C | 2007 ASAT, **865 km, i = 98.8°**, 3,500+ tracked | [Orbital Radar](https://orbitalradar.com/events/fengyun-1c-2007); NASA ODQN | **V** for orbit and tracked count |
| Fengyun-1C ≥1 cm members | ~35,000–40,000 | commonly cited modelled estimate | **P / UNVALIDATED** — a model output, not an observation |
| Iridium-33 / Cosmos-2251 | 2009 collision, **789 km**, i = 86.4° / 74.0°, 628 / 1,668 catalogued | [Wikipedia 2009 collision](https://en.wikipedia.org/wiki/2009_satellite_collision); NASA ODQN | **V** for orbit and catalogued counts |
| Cosmos-1408 | 2021 ASAT, **465–490 km, i = 82.5°**, ≥1,500 tracked | [EU SST](https://www.eusst.eu/newsroom/news/eu-sst-confirms-fragmentation-space-object-cosmos-1408); Pardini & Anselmo, *Acta Astronautica* 2023;210:465 | **V** |
| Per-cluster ≥1 cm counts other than FY-1C | scaled from tracked counts | **E** — this study's scaling, not a published figure |
| SBM ejection-velocity distribution | log₁₀(Δv) ~ N(0.9χ+2.9, 0.4) | Johnson et al. 2001 (as above) | **P** |
| SL-8 / SL-16 upper-stage masses | 1.4 t / 8.3 t | standard launch-vehicle data | **P** |
| Satellite build cost | $50k–400k per kg dry | smallsat-to-exquisite benchmark range | **E** — see v9 addition below: the $50k/kg figure used as this study's primary (not sensitivity-only) platform cost is below every real benchmark located in this revision's literature search |
| Fleet learning curve | 85% | standard aerospace convention | **E** |
| Active LEO fleet / mean asset value | 11,000 / $2M–30M swept | **E** — parametric; see R12 |

**Nothing in the v6 additions is presented as SOURCED when it is derived or
estimated.** The two numbers a reviewer should attack first are the Fengyun-1C
≥1 cm member count (modelled) and β at 619 m/s (derived from a scaling law
extrapolated below its calibration) — though note that neither drives the v6
verdict, which rests on the Δv-per-engagement invariant, a purely orbital-
mechanical result.

---

## Laser comparison (LAS-1…3)

| Quantity | Value | Source | Status |
|---|---|---|---|
| Momentum coupling C_m, peak, ns pulses on metal | 1×10⁻⁴ N·s/J | Phipps CR et al., laser-propulsion / ORION literature | **P** |
| Ablation threshold fluence, Al | 1 J/cm² (range 0.5–5) | laser-ablation literature | **P** — the cliff's *position* moves with it, its *existence* does not |
| Optimum-coupling fluence | 5 J/cm² | same | **P** |
| C_m roll-off above optimum | ∝ Φ^(−1/2) | same | **P** |
| Spot diameter, 2.44·λ·R/D | — | standard diffraction | **DERIVED** |
| Kasten–Young airmass | — | standard | **P** |
| Fried parameter r₀, zenith, 500 nm | 15 cm (good site) | standard site statistics | **P** |
| r₀ airmass scaling, X^(−3/5); λ^1.2 | — | standard turbulence theory | **P** |
| AO fitting error 0.34(d_act/r₀)^(5/3) | — | standard AO error budget | **P** |
| Actuator count | 1,000 across the aperture | facility-class AO benchmark | **E** — chosen *generously* to the laser |
| Atmospheric extinction τ | 0.15 at zenith | good astronomical site | **P** |
| LLNL Mercury | 100 J @ 10 Hz | demonstrated | **P** |
| DiPOLE-100 (STFC) | 100 J @ 10 Hz | demonstrated | **P** |
| HAPLS / ELI-Beamlines | 200 J @ 10 Hz | demonstrated | **P** |
| **Demonstrated max average power, repetitive ns** | **2 kW** | the three above | **P** — sets the readiness benchmark (AM-9) |
| Solar constant / solar V magnitude | 1361 W/m² / −26.74 | standard | **P** |
| Debris albedo | 0.1 | conventional assumption | **E** |
| Laser facility cost | $157M–1.5B | parametric (telescope + $50/W laser) | **E** |
| Terminator window / weather | 3–4 h/night, 60–70% clear | site statistics | **P** |
| Published laser cost being corrected | $100–500/kg | v5 plan §1.5 | **P** — shown here to be a marginal-energy figure omitting the facility |

The laser model was parameterised **in the laser's favour** wherever a choice
existed (generous AO actuator count, good site, optimum-coupling fluence,
best-case aperture). The 46× technology gap and the 536–8,031× cost correction
are therefore lower bounds on the laser's difficulty.
## v7 addition — laser-ablation accuracy comparison (SIM-14)

| Quantity | Value | Source | Status |
|---|---|---|---|
| Debris relative speed | up to 15 km/s | `laser accuracy comm.md` (informal brief provided for this study; not independently verified against primary literature) | **E** |
| Orbit-prediction requirement | "better than 1 m" | same | **E** |
| Orbit-prediction optimistic claim | "<1 cm within 10 s of detection" | same, attributed by the brief to unnamed "research" | **E** |
| Pointing requirement | "sub-arcsecond" (no exact figure given) | same | **E** — this study assigns a 0.1–2.0 arcsec envelope, its own estimate |
| p-N correction, if omitted | bias "on the order of the size of the debris objects themselves" | same | **E** — implemented literally as a residual equal to target radius |
| Engagement stand-off range | 500 km (fixed across all SIM-14 scenarios) | **not given in the brief**; this study's own representative assumption for a space-based laser-ablation concept | **E** |

**Key caveat, stated plainly:** `laser_accuracy_comm.md` is a secondary brief, not a peer-reviewed source, and gives ranges/qualitative claims rather than point estimates for several quantities (pointing, in particular). SIM-14 uses its stated claims at face value and flags every figure this study had to supply itself (range, exact pointing value) as its own assumption, not the brief's. SIM-14's scope is pointing/prediction error tolerance only; it does not reopen SIM-9's cost comparison, where FG01 remains 3–4 orders of magnitude more expensive per object than ground-based laser ablation.

---

## v8 additions — launcher structural detail (COIL-02, COIL-03, GUN-02)

| Quantity | Value | Source | Status |
|---|---|---|---|
| Coilgun structural wall mass share | 5% of launcher mass (1.7 kg of 34.1 kg, 31-stage COIL-01 design point) | **DERIVED** in COIL-02 from magnetic-pressure sizing (34.2 MPa) and standard laminate/gauge design rules | — |
| Coil banding thickness | 0.93 mm S-glass/epoxy, SF 2 | **DERIVED** in COIL-02 | — |
| Bore-tube insulation requirement | must be non-conductive (G10/FR4/PEEK); a metal tube shorts the armature coupling | standard induction-coilgun design constraint | **P** |
| Radial envelope, 30 mm bore | 52 mm OD, winding is 47% of wall thickness (2 layers, AWG 10) | **DERIVED** in COIL-03 | — |
| Chemical-gun barrel, 100 g @ 619 m/s, 30 mm bore | 89.1 MPa peak pressure; tapered Ti-6Al-4V barrel 0.72 kg, 47% lighter than constant-wall | **DERIVED** in GUN-02, standard interior-ballistics/thin-wall-pressure-vessel sizing | — |

These are structural-sizing refinements of the launcher case already reported in Section 3.8 (GUN-01/COIL-01/CHK-01); they do not change any cost, mass-crossover, or feasibility figure reported elsewhere, and are not a re-derivation of the 400 m/s mass-crossover shot count.

## v9 additions — engineering-assumption precedent search (platform cost, pointing, muzzle-speed repeatability)

Located and checked in response to a direct request to justify L3 and L9 of `LIMITATIONS.md` more thoroughly. Findings cut both ways — one assumption gained real support, one did not — and both are reported regardless of which direction they point.

| Quantity | Value | Source | Status |
|---|---|---|---|
| Structure mass fraction, 20% of dry mass | Matches this study's own convention exactly | Ritter B, Meskers AJH, Miles O, et al. A space weather information service based upon remote and in-situ measurements of coronal mass ejections heading for Earth. Submitted to *Journal of Space Weather and Space Climate*; arXiv:1502.01846 (2015). Table 5: "Structure (20% of dry mass) 66.01 kg." | **V** — full PDF read directly in this session; figure quoted verbatim from the table |
| Demonstrated deep-space laser-terminal pointing stability | 2.5 µrad achieved (4.2 µrad RMS requirement) | NASA/MIT Lincoln Laboratory Lunar Laser Communication Demonstration (LLCD) mission documentation, located via web search | **P** — located by search summary in this session, not independently read in full; specific authorship not confirmed |
| Demonstrated CubeSat-scale laser-terminal pointing | ~500 µrad (0.5 mrad) in orbit | NASA/Aerospace Corporation OCSD (Optical Communications and Sensor Demonstration) 1.5U CubeSat program, located via web search | **P** — same caveat as above |
| Railgun shot-to-shot muzzle-velocity deviation | 0.92%–5.6%, depending on armature mass and burst | He Y, Song S, Guan Y, Cheng C, Dai W, Qiu X, Li Y. An investigation into muzzle velocity repeatability of a railgun. *IEEE Transactions on Plasma Science* 2015;43:1647–1651. | **P** — bibliographic details and figures confirmed via search-engine summary of the paper, not independently read in full text |

**Honest result, stated plainly.** The pointing assumption (200 µrad) is now bracketed by real demonstrated hardware, sitting inside the achievable range but on the side that requires LLCD-class actively-stabilized optics rather than CubeSat-cheap hardware — a real strengthening, with a real caveat (which class of hardware, and its cost, is not separately budgeted). The muzzle-speed repeatability assumption (0.1%) is **not** strengthened by this search — the closest available literature data (railguns) is 9–56× worse than assumed — though a stated, technically grounded reason exists for why that specific data may not transfer (railguns have a sliding-contact/friction noise source a contactless coilgun does not share). No coilgun-specific repeatability measurement was found. This is reported as an open question, not resolved in either direction. See `LIMITATIONS.md` L3 and L9 for the full discussion.

## v9 addition — real-world platform build-cost benchmark (checked against, not confirmed)

| Quantity | Value | Source | Status |
|---|---|---|---|
| Modern (post-1990) small-satellite bus cost | ~$100,000/kg dry mass (FY97 dollars) | Bearden DA. Small-satellite costs. *Crosslink* (The Aerospace Corporation) 2000/2001;Winter:33–44. Fig. "Dollars-per-kilogram comparison" — 100 (modern small sats), 150 (traditional small sats), 500 (DOD large sats), thousand FY97 dollars/kg, bus cost only | **V** — full article read in this session |
| Traditional NASA mission cost | ~$900,000/kg launch mass (FY98 dollars); faster-better-cheaper NASA missions ~$120,000/kg | same, Fig. "dollars-per-kilogram ratio" (total mission cost, not bus-only) | **V** |
| General spacecraft cost anchor | "often cost more than $200 thousand per kilogram and could reach $1 million per kilogram with delivery-to-space costs included" | same, quoted directly | **V** |

**This is an uncomfortable finding, reported without softening.** This study's headline $21.2M costed platform — the figure behind the primary $0.9–4.2M/object result, not the acknowledged $500k lower-bound sensitivity case — rests on $50k/kg build cost. Every real benchmark located in this search, for the cheapest demonstrated category (modern, post-1990, faster-better-cheaper small satellites), sits at $100,000–120,000/kg or above, in late-1990s dollars not adjusted for inflation to the present. This study's primary platform-cost assumption is therefore below the cheapest real-world benchmark found, not merely unvalidated in the neutral sense the earlier entry implied. This does not change this study's qualitative conclusion — the concept is already found not economically competitive at $21.2M, and a more realistic (higher) platform cost only widens that gap — but it does mean the "$21.2M costed design" should not be read as a conservative or even neutral estimate relative to demonstrated small-satellite cost history; it is optimistic. The bus-cost-only comparison (Bearden's modern-smallsat figure) is the fairer read-across, since FG01's $50k/kg is also described as covering the bus-equivalent subsystems (structure, avionics, power, propulsion, thermal) rather than a full mission cost; even on that narrower comparison, the gap is real. This is not a re-derivation of Section 3.11's figures — it is a flag that they should be revisited against this benchmark before being called a costed, realistic baseline.

## Summary for a reviewer

- **Verified in-session (V):** the 40 J/g threshold and its regime definitions;
  the Ganguli/NRL tungsten prior art in full; Starshine 1's geometry, orbit,
  reentry date, **and now mass** (confirmed 39 kg against eoPortal, resolving
  the prior weakest-link status — see v9 addition); the NRLMSISE-00
  implementation version and the table's fidelity to it; the 20%-of-dry-mass
  structure fraction (Ritter et al. 2015) and the real small-satellite $/kg
  cost benchmarks (Bearden 2000) added in this revision.
- **Highest-priority unverified/unfavorable entries:** the **2026 rhenium spot
  price** and the **SBM non-catastrophic effective-mass convention** remain
  unverified (P). More consequential than either: this study's **$50k/kg
  platform build-cost figure**, used as the primary (not sensitivity-only)
  cost basis, is now checked against real small-satellite cost history and
  found to sit *below* every real benchmark located — not merely unsourced,
  but optimistic. This does not change the study's verdict (a higher real cost
  only strengthens "not competitive") but means the headline $21.2M / $0.9–4.2M
  figures should be read as a lower bound, not a conservative center. Muzzle-
  speed repeatability (0.1%) also remains unresolved: the closest analog data
  (railguns, not coilguns) is 9–56× worse than assumed, with a stated but
  unconfirmed reason it may not transfer.
- **No number in this study is presented as sourced when it is an estimate.**
  All **E** entries are tagged UNVALIDATED in the outputs and are listed in
  `LIMITATIONS.md`.
