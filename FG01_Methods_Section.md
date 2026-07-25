# 2. Materials and methods

## 2.1 Analytical framework, falsification gates, and traceability

The analysis is implemented as a linked pipeline of thirty-four simulation modules organised into five families: the core physics and system sequence (SIM-0 – SIM-14), the cloud–debris interaction model (INT-0 – INT-6), the engagement-economics model (ECO-1 – ECO-5), the comparative laser benchmark (LAS-1 – LAS-4), and launcher-specific ballistics (GUN-1, COIL-1). Modules share a common library providing orbital mechanics, atmospheric density, relative motion, breakup modelling, interaction physics, and debris-population routines, so that every downstream module reads upstream results at runtime rather than through re-entered constants.

The architecture is deliberately falsification-first. Rather than sizing a design and checking constraints afterwards, hard physical boundaries are evaluated *before* any design work: the catastrophic-disruption boundary (Section 2.4) and the orbital-decay requirement (Section 2.7) are each formulated as pass/fail gates that a candidate operating point must clear before guidance, launcher, platform, or cost analysis proceeds. Where no design point satisfies a gate, that outcome is reported as a primary result rather than resolved by parameter adjustment.

Every quantitative input carries one of three traceability tags: **DERIVED** (obtained algebraically or numerically from stated, non-controversial physics), **SOURCED** (traceable to a specific citation), or **UNVALIDATED** (an engineering assumption for which no direct experimental or literature source exists at the relevant length- and velocity-scale). Inputs carrying uncertainty are propagated as distributions rather than point values (Section 2.12), and results dependent on UNVALIDATED inputs are reported as envelopes. The complete pipeline is reproducible from a single command with pinned dependency versions and a fixed pseudorandom seed.

## 2.2 Scope: target classes and deorbit modes

Two target classes are analysed, because the governing physics is mass-independent while the delivery economics are not.

**Small debris** comprises solid aluminium-6061 spheres of 0.5–10 cm diameter, representing the untracked fragment population that lies below the reliable ground-surveillance threshold in low Earth orbit yet retains sufficient mass and closing velocity to disable an operational spacecraft [1,2]. **Mid-mass objects** comprise intact tracked bodies of approximately 10²–10⁴ kg, representing spent upper stages and defunct satellites of the class targeted by conventional active debris removal [3]. The specific impact energy governing catastrophic disruption (Section 2.4) is independent of target mass, so a single interaction model spans both classes; what differs is the areal efficiency of cloud delivery (Section 2.5.2), which improves as target cross-section grows, and the projectile mass required, which scales linearly with target mass.

Two deorbit modes are evaluated against both classes.

**Direct reentry** lowers the target perigee to 100 km in a single engagement, producing reentry within minutes. **Accelerated reentry** lowers the orbit by a finite decrement Δh (evaluated at 150, 200 and 250 km), after which natural atmospheric drag completes the deorbit over months to years. The two modes impose substantially different impulse requirements at the same altitude, and therefore different specific impact energies; both are evaluated across the full altitude range so that the trade between impulse magnitude and post-engagement orbital lifetime is resolved rather than assumed. Success is defined differently for each: perigee below 100 km for direct reentry, and post-engagement orbital lifetime below the applicable regulatory threshold for accelerated reentry (Section 2.7).

## 2.3 Orbital mechanics and the deorbit impulse requirement

For a target in a circular orbit at altitude *h*, the velocity decrement required at apogee to lower perigee to a target altitude *h*<sub>p</sub> follows from the vis-viva relation [4]:

Δ*v*(*h*, *h*<sub>p</sub>) = √[μ/(*R*<sub>E</sub>+*h*)] − √{μ[2/(*R*<sub>E</sub>+*h*) − 2/(2*R*<sub>E</sub>+*h*+*h*<sub>p</sub>)]}  (1)

where μ is Earth's gravitational parameter and *R*<sub>E</sub> the mean equatorial radius [DERIVED]. Direct reentry sets *h*<sub>p</sub> = 100 km; accelerated reentry sets *h*<sub>p</sub> = *h* − Δ*h*. Equation (1) is evaluated across *h* = 400–1200 km. Because the required impulse for a fixed altitude decrement decreases with altitude while the impulse for direct reentry increases, the two modes diverge substantially over the tested range, and this divergence propagates directly into the disruption gate through Equation (3).

Momentum conservation for a projectile of mass *m*<sub>p</sub> striking a target of mass *m*<sub>t</sub> at relative velocity *v*<sub>rel</sub>, with momentum-enhancement factor β accounting for ejecta-augmented coupling, gives

*m*<sub>t</sub> Δ*v* = β *m*<sub>p</sub> *v*<sub>rel</sub>  (2)

from which the required projectile mass and the design mass ratio *f*(*h*) = *m*<sub>p,design</sub>/*m*<sub>t</sub> follow, with an explicit engineering margin γ applied to projectile mass. The mass ratio is independent of target mass and of engagement range; it depends only on altitude (through Δ*v*), on *v*<sub>rel</sub>, and on β.

## 2.4 Catastrophic disruption gate

Substituting the minimum projectile mass from Equation (2) into the specific-energy definition *E*<sub>s</sub> = ½*m*<sub>p</sub>*v*<sub>rel</sub>²/*m*<sub>t</sub> — the governing quantity of the NASA Standard Breakup Model [5] — yields

*E*<sub>s</sub>(*h*, *v*<sub>rel</sub>, β) = γ Δ*v*(*h*) *v*<sub>rel</sub> / (2β)  (3)

Equation (3) is notable for being independent of target mass, which is why a single gate serves both target classes. A collision is treated as catastrophic when *E*<sub>s</sub> equals or exceeds 40 kJ kg⁻¹, evaluated across the 30–50 kJ kg⁻¹ range to span the threshold's reported uncertainty [5,6]. Inverting Equation (3) gives the maximum relative velocity admissible at a given altitude without crossing into the catastrophic regime. This gate is evaluated first, for both deorbit modes, before any subsequent module executes; the substitution in Equation (3) is verified algebraically and Equation (1) is cross-checked against an independent Hohmann-transfer formulation.

## 2.5 Cloud–debris interaction model

Treating a dispersed projectile swarm as a single macroscopic impactor of equal total mass is a modelling assumption, not a validated equivalence. The interaction model therefore resolves the swarm explicitly, so that bulk momentum delivery and per-impact specific energy are computed as separate quantities rather than conflated.

### 2.5.1 Momentum coupling from crater-ejecta scaling

The momentum-enhancement factor β is derived rather than assumed, using crater-ejecta scaling theory [7,8]. Momentum enhancement above unity requires target material to be ejected at speeds exceeding the impactor's own momentum contribution; the mass of such ejecta is governed by the ratio of impact velocity to the target's strength velocity. For aluminium-6061 the strength velocity is of order 3×10² m s⁻¹, so an impact at the engagement velocities admitted by Equation (3) occurs at only a small multiple of that value, and the ejecta contribution is correspondingly modest. The projectile predominantly embeds, penetrating several times its own diameter, which represents the physical floor β = 1. The same scaling evaluated at conventional hypervelocity conditions (≈6 km s⁻¹) reproduces the momentum enhancements of β ≈ 2–3 reported in the impact and asteroid-deflection literature [8,9], which serves as the model's consistency check. Because no direct measurement exists for rhenium- or tungsten-on-aluminium impacts at centimetre scale and sub-kilometre-per-second velocity, β retains an UNVALIDATED tag and is propagated as a distribution.

### 2.5.2 Cloud areal density and mass efficiency

Following release, the projectile swarm's spatial evolution is modelled using linearised Clohessy–Wiltshire relative motion in the Hill frame [10], propagated to the intercept plane to yield a two-dimensional areal number-density field. Only projectile mass incident on the target silhouette transfers useful momentum, so the governing quantity is the mass efficiency ε, the ratio of on-target mass to launched mass. For a Gaussian transverse profile and a target small relative to the cloud, ε scales as the square of the ratio of target diameter to cloud diameter, with a coefficient approximately three times more favourable than the uniform-disc idealisation. This quantity is where the two target classes diverge most sharply: the areal penalty that dominates engagement of a centimetre-scale fragment becomes negligible for a metre-scale body.

### 2.5.3 Aggregate delivery and aim tolerance

Total delivered impulse is obtained by summing grain-level coupling over the incident population, giving Δ*v*<sub>imparted</sub> = ε *m*<sub>L</sub> β̄ *v*<sub>rel</sub>/*m*<sub>t</sub> for launched mass *m*<sub>L</sub>. Aim error is treated not as a binary hit-or-miss outcome but as a continuous reduction in ε: for a Gaussian cloud displaced from the target centre, on-target mass falls off exponentially with the square of the miss distance, so holding delivered impulse constant requires launched mass to increase by the reciprocal factor. Engagement reliability is therefore purchased with launched mass along a smooth, finite penalty curve rather than encountering a discontinuous failure threshold. Launched mass scales as the square of the miss distance, making launcher pointing accuracy and muzzle-velocity repeatability the dominant engineering sensitivities; both are treated as UNVALIDATED engineering estimates and swept.

### 2.5.4 Secondary ejecta and net debris

Because each engagement comprises many individual impacts, secondary-ejecta production is evaluated at the aggregate level. Per-impact crater volume from the scaling relations of Section 2.5.1 is summed over the incident population, and the resulting ejecta mass distribution is classified against the 1 mm and 1 cm trackability thresholds using the non-catastrophic (cratering) branch of the Standard Breakup Model [5], rather than the catastrophic-fragmentation branch applicable above the disruption threshold. Net debris flux is computed as objects removed minus trackable fragments generated minus projectile mass remaining in orbit.

## 2.6 Engagement geometry and terminal guidance

Engagement geometry is resolved in three dimensions rather than as a scalar closing speed, because the two differ materially in their consequences. Decomposing the relative velocity into components shows that only the component anti-parallel to the target's orbital velocity produces the along-track impulse that lowers the orbit; a crossing engagement expends the remainder on a cross-track impulse that rotates the orbital plane without reducing altitude, while incurring the full specific energy of the encounter in Equation (3). Delivering a given orbit-lowering impulse from a crossing geometry therefore requires substantially more specific energy than from a co-orbital geometry, and for realistic crossing angles this places the engagement well inside the catastrophic regime. Crossing engagements are consequently excluded on kinematic grounds before guidance performance is considered, and the admissible engagement geometry is restricted to a narrow retrograde, in-plane cone.

Within that cone, terminal accuracy is computed from an explicit error budget combining launcher pointing uncertainty, muzzle-velocity repeatability, target state uncertainty over the tracking arc, and control-loop latency. The transverse miss contribution is evaluated using the true cross-line-of-sight rate obtained from the relative-motion solution, not the scalar closing velocity; because the admissible geometry is co-orbital, this rate is small and the latency contribution to miss distance is correspondingly modest across the tested range. Miss distance is then mapped to launched-mass penalty through Section 2.5.3.

## 2.7 Orbital decay and reentry

### 2.7.1 Atmosphere model

Atmospheric density is obtained from the NRLMSISE-00 empirical model [11], evaluated at solar-minimum, nominal and solar-maximum conditions (F10.7 = 70, 150, 250). To make repeated propagation tractable, density is pre-tabulated on an altitude grid and interpolated; interpolation error against direct model evaluation is verified to remain below 0.02%.

### 2.7.2 Propagation

Post-engagement orbital lifetime is computed by two independent methods. Near-circular orbits are integrated using the orbit-averaged decay relation d*a*/d*t* = −*C*<sub>D</sub>(*A*/*m*)ρ(*h*)√(μ*a*) [12], while eccentric post-impulse orbits — which traverse several orders of magnitude of density variation per revolution — are propagated by full three-dimensional Cartesian integration including J2 and atmospheric co-rotation. The two methods are cross-validated against one another; comparison is performed on mean rather than osculating elements, since J2 induces kilometre-scale periodic oscillations in osculating semi-major axis that would otherwise mask genuine agreement. Model output is additionally compared against the documented decay history of an object of known geometry and area-to-mass ratio.

Orbital lifetime is evaluated for target debris in both deorbit modes, and separately for projectile mass that misses its intended target. For missed mass, the post-launch trajectory is computed from launch geometry rather than assumed: a retrograde launch at relative velocity exceeding the direct-reentry impulse requirement places the missed projectile on a trajectory whose perigee lies below the sensible atmosphere, so its orbital persistence is determined by that geometry rather than by its ballistic coefficient alone.

### 2.7.3 Reentry ablation

Reentry survival is evaluated using rate-limited stagnation-point convective heating after Fay–Riddell theory [13], rather than a total-energy ablation budget, since reentry heating is rate- rather than energy-limited. The ablation threshold is set at the melting transition with aerodynamic melt-stripping, consistent with the definition of the effective ablation enthalpies used for metallic materials, rather than at the vaporisation point. Ablated mass fraction, surviving ballistic coefficient and ground-casualty expectation are computed for both target and projectile materials across the tested size range, and compared against the accepted 10⁻⁴ per-event casualty threshold [14].

## 2.8 Electromagnetic launcher

Launcher requirements are computed from standard electromagnetic-launch energy relations, with kinetic energy ½*m*<sub>p</sub>*v*² converted to required electrical energy through a stage efficiency swept over 0.3–0.5, and stage count, bore length, conductor gauge and capacitor mass derived from established coilgun scaling relations [15]. Required performance is compared against demonstrated laboratory and test-article benchmarks rather than projected capability, so that the achievable mass–velocity envelope constrains the design space rather than being fitted to it.

## 2.9 Debris population and engagement cadence

Engagement opportunity is modelled from orbital mechanics rather than assumed. Fragments originating in a common breakup event share an inclination but disperse in right ascension of the ascending node, because differential nodal precession under J2 varies with semi-major axis and inclination [4,16]; the resulting spread is computed as a function of time since the parent event. A platform can engage only those objects whose geometry falls within the admissible cone of Section 2.6, so opportunity counting combines the population's spatial density with the relative drift rate and the cone geometry.

Reaching successive targets, however, requires the platform to match each target's orbit. The propulsive cost per engagement is derived in closed form from the nodal-precession rate, orbital velocity and mission duration, and validated against numerical tour simulation. This quantity is central to the economics because opportunity count and per-engagement propulsive cost both scale linearly with the width of the altitude band serviced, so the two effects cancel and the per-engagement cost is invariant with respect to that design choice. Achievable engagement count over a mission is therefore set by the platform's total velocity budget divided by this per-engagement cost. Population inputs are taken from established debris-environment models [2,17], with the caveat that the majority of the ≥1 cm population is statistically modelled rather than individually catalogued.

## 2.10 Cost model

The cost model is formulated parametrically, as a structural relation between inputs rather than an evaluation at a single design point. Cost per object removed is

*C*<sub>object</sub> = (*C*<sub>platform</sub> + *C*<sub>ops</sub>)/*N*<sub>eng</sub> + *m*<sub>L</sub>*P*<sub>material</sub> + *E*<sub>shot</sub>*P*<sub>elec</sub>  (4)

where *C*<sub>platform</sub> is derived from a subsystem mass budget (structure, avionics, power, launcher, propulsion, thermal control) combined with a build cost per unit dry mass and a launch price per unit wet mass, with propellant sized by the Tsiolkovsky relation for the velocity budget of Section 2.9; *N*<sub>eng</sub> is the achievable engagement count; *m*<sub>L</sub> is launched projectile mass per engagement from Section 2.5.3; and *E*<sub>shot</sub> is launcher electrical energy from Section 2.8.

Each input is swept across a stated range rather than fixed: dry mass, build cost per kilogram, launch price per kilogram, fleet size with an associated learning-curve exponent, velocity budget, and material unit price. Results are therefore reported as surfaces over this parameter space, with the sensitivity of *C*<sub>object</sub> to each input reported explicitly (Section 2.12). Because Equation (4) is dominated by its first term wherever *N*<sub>eng</sub> is small, the model's structure implies that platform mass and achievable engagement count, rather than consumable cost, govern the outcome; revised subsystem mass estimates or launch prices substitute directly into the same relation without altering its form. Material prices are dated at the point of use, since the unit price of rhenium moved by a factor of approximately three over the period spanned by the analysis [18].

## 2.11 Comparative benchmarking

Candidate remediation methods are compared on a like-for-like basis using published performance and cost assessments, with each competing approach evaluated for applicability to the target class in question: laser ablation and other contactless methods for the small-debris class [19,20], and rendezvous-and-capture removal vehicles for the mid-mass class [3,21]. Because the concept under study and several of its comparators are at differing technology-readiness levels, comparisons state the maturity of each method alongside its cost estimate, and terminal-accuracy requirements are computed on a common basis rather than taken from source claims.

## 2.12 Uncertainty quantification

All parameters carrying a stated uncertainty range are propagated jointly through the full pipeline by Monte Carlo sampling with 10⁶ draws and a fixed seed, rather than by one-at-a-time sweeps, which do not correctly characterise variance in a nonlinear, multi-constraint system [22]. Viability is defined as the joint satisfaction of the constituent constraints — sub-catastrophic specific energy, achievable delivery accuracy, compliant post-engagement orbital lifetime, non-negative net debris flux, and a stated cost criterion — evaluated simultaneously on each sample. First-order and total-effect Sobol sensitivity indices [23,24] are computed to attribute output variance among inputs, with convergence diagnostics reported. Where a cost criterion expressed per unit mass is not meaningful for very small targets, a per-object criterion is reported alongside it rather than substituted silently.

## 2.13 Verification and validation

Each module carries verification appropriate to its physics: Equation (1) is cross-checked against an independent Hohmann formulation and the Equation (3) substitution verified algebraically; the relative-motion solution is checked against analytical solutions; the two orbital-decay propagators are cross-validated on mean elements and compared against an observed decay history; the crater-ejecta scaling is verified to reproduce published momentum enhancement at hypervelocity conditions; the closed-form per-engagement velocity cost is compared against numerical tour simulation; atmospheric interpolation is verified against direct model evaluation; and Monte Carlo convergence is confirmed by sample-size scaling. Verification residuals are reported with each result rather than aggregated, and all external numerical inputs are recorded in a citation log stating, for each value, whether it was independently verified against a primary source, inherited from prior work without independent verification, or estimated.

## 2.14 Mapping of methods to research questions

| Research question | Governing method |
|---|---|
| Viability of momentum-transfer removal | §2.4, §2.5, §2.12 |
| Dispersed cloud versus discrete grains; grain sizing | §2.5.1–2.5.2 |
| Sufficiency of imparted impulse | §2.3, §2.7.2 |
| Maximum serviceable altitude | §2.7.2 |
| Launch geometry and angle | §2.6 |
| Projectile mass per unit target mass, *f*(*h*) | §2.3, §2.5.3 |
| Delivery accuracy and its consequence | §2.5.3, §2.6 |
| Secondary debris generation | §2.5.4 |
| Energy per engagement; launcher sizing | §2.8 |
| Cost per object and per unit mass | §2.9, §2.10 |
| Comparison against alternative methods | §2.11 |
| Reentry survival of target and projectile | §2.7.3 |
| Regulatory and environmental compliance | §2.7.2, §2.7.3 |

---

## References

1. Liou, J.-C., Johnson, N. L. (2006). Risks in space from orbiting debris. *Science*, 311(5759), 340–341.
2. Klinkrad, H. (2006). *Space Debris: Models and Risk Analysis*. Springer, Berlin.
3. Shan, M., Guo, J., Gill, E. (2016). Review and comparison of active space debris capturing and removal methods. *Progress in Aerospace Sciences*, 80, 18–32.
4. Vallado, D. A. (2013). *Fundamentals of Astrodynamics and Applications*, 4th ed. Microcosm Press, Hawthorne, CA.
5. Johnson, N. L., Krisko, P. H., Liou, J.-C., Anz-Meador, P. D. (2001). NASA's new breakup model of EVOLVE 4.0. *Advances in Space Research*, 28(9), 1377–1384.
6. Krisko, P. H. (2011). Proper implementation of the 1998 NASA breakup model. *Orbital Debris Quarterly News*, 15(4), 4–5.
7. Holsapple, K. A. (1993). The scaling of impact processes in planetary sciences. *Annual Review of Earth and Planetary Sciences*, 21, 333–373.
8. Housen, K. R., Holsapple, K. A. (2011). Ejecta from impact craters. *Icarus*, 211(1), 856–875.
9. Cheng, A. F., et al. (2023). Momentum transfer from the DART mission kinetic impact on asteroid Dimorphos. *Nature*, 616, 457–460.
10. Clohessy, W. H., Wiltshire, R. S. (1960). Terminal guidance system for satellite rendezvous. *Journal of the Aerospace Sciences*, 27(9), 653–658.
11. Picone, J. M., Hedin, A. E., Drob, D. P., Aikin, A. C. (2002). NRLMSISE-00 empirical model of the atmosphere: Statistical comparisons and scientific issues. *Journal of Geophysical Research: Space Physics*, 107(A12), 1468.
12. King-Hele, D. (1987). *Satellite Orbits in an Atmosphere: Theory and Applications*. Blackie, Glasgow.
13. Fay, J. A., Riddell, F. R. (1958). Theory of stagnation point heat transfer in dissociated air. *Journal of the Aeronautical Sciences*, 25(2), 73–85.
14. NASA (2019). *NASA Procedural Requirements for Limiting Orbital Debris and Evaluating the Meteoroid and Orbital Debris Environments*, NPR 8715.6B. National Aeronautics and Space Administration, Washington, DC.
15. McNab, I. R. (2003). Launch to space with an electromagnetic railgun. *IEEE Transactions on Magnetics*, 39(1), 295–304.
16. Vallado, D. A., Finkleman, D. (2014). A critical assessment of satellite drag and atmospheric density modeling. *Acta Astronautica*, 95, 141–165.
17. Flegel, S., et al. (2011). The MASTER-2009 space debris environment model. *Proceedings of the 5th European Conference on Space Debris*, ESA SP-672.
18. U.S. Geological Survey (2025). *Mineral Commodity Summaries 2025*. U.S. Geological Survey, Reston, VA.
19. Phipps, C. R., et al. (2012). Removing orbital debris with lasers. *Advances in Space Research*, 49(9), 1283–1300.
20. Ganguli, G., Crabtree, C., Rudakov, L., Chappie, S. (2012). A concept for elimination of small orbital debris. *Transactions of the Japan Society for Aeronautical and Space Sciences, Aerospace Technology Japan*, 10(ists28), Pr\_2\_7–Pr\_2\_13.
21. Forshaw, J. L., et al. (2016). RemoveDEBRIS: An in-orbit active debris removal demonstration mission. *Acta Astronautica*, 127, 448–463.
22. Saltelli, A., et al. (2008). *Global Sensitivity Analysis: The Primer*. John Wiley & Sons, Chichester.
23. Sobol', I. M. (2001). Global sensitivity indices for nonlinear mathematical models and their Monte Carlo estimates. *Mathematics and Computers in Simulation*, 55(1–3), 271–280.
24. Saltelli, A., et al. (2010). Variance based sensitivity analysis of model output. Design and estimator for the total sensitivity index. *Computer Physics Communications*, 181(2), 259–270.

> **Reference verification note.** Consistent with the project's provenance requirements, every reference above should be confirmed against the primary source and recorded in the citation log before submission, with volume, page and DOI verified individually. References [1]–[5], [7]–[13], [15], [19], [20], [22]–[24] correspond to standard, widely cited works in their respective fields; references [6], [16], [17], [21] should be checked for the most current edition or superseding publication, and [14] and [18] for the applicable revision and year of issue at the time of submission.
