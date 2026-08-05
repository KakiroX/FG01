# FG01/KIDR: Technical Summary

## 1. Problem statement

Low Earth orbit (LEO) contains an estimated 1.2 million objects with characteristic length ≥1 cm, of which only ~3% are individually tracked by ground-based space surveillance networks. At typical LEO relative encounter velocities (7–15 km/s), the kinetic energy of even a centimeter-scale fragment is sufficient to catastrophically disable an operational spacecraft on impact, since kinetic energy scales as $v^2$ and orbital velocities are an order of magnitude above terrestrial projectile speeds. Because this population lies below the size threshold at which individual cataloging is feasible, existing active debris removal (ADR) architectures — robotic-arm capture, harpoon/net rendezvous, electrodynamic tether deorbit — cannot be applied, since all require a tracked, individually approachable target.

## 2. Concept

FG01 proposes targeted momentum transfer without rendezvous: a dedicated orbital platform carrying an electromagnetic launcher (coilgun) fires a briefly-dispersing cloud of rhenium grains at a terminally-tracked (but not catalogued) debris object. The impulse is sized not to force immediate atmospheric reentry, but to lower the object's orbit by a fixed decrement (~200 km), after which enhanced atmospheric drag completes deorbit over a shortened but still finite timescale — an "accelerated reentry" architecture, distinct from a direct-impact deorbit design.

### 2.1 Energy constraint (catastrophic-disruption threshold)

Impact specific energy $E_s$ (kinetic energy per unit target mass) must remain below the catastrophic-fragmentation threshold ($E_{s,c} \approx 40$ J/g, per the NASA Standard Breakup Model), above which the target fragments into a population of new trackable-and-untrackable debris rather than being displaced as a single body. The design point operates at $E_s = 26.7$ kJ/kg, i.e., 0.67× the threshold — a derived margin, not an assumed one.

### 2.2 Momentum coupling

Momentum delivered to the target exceeds the incident projectile momentum alone, due to recoil from crater ejecta expelled opposite the impact direction (a momentum-enhancement factor β, by Newton's third law applied to the ejecta stream). β was not assumed a priori; it was derived from Housen–Holsapple crater-ejecta scaling, giving a central value of β ≈ 1.201 at the 619 m/s design-point velocity — an embedding-dominated rather than ejecta-dominated coupling regime, since the impact velocity is only marginally above the target material's strength velocity.

### 2.3 Engagement geometry

Terminal aiming precision is governed by relative velocity at intercept, not absolute velocity. A crossing-orbit engagement geometry (platform and target on intersecting trajectories) implies a relative velocity of several km/s, which is kinematically incompatible with the pointing/timing precision achievable by any modeled terminal guidance system. A co-orbital engagement geometry — platform matched to the target's orbital plane and altitude prior to intercept, reducing relative velocity to ~10 m/s — is therefore not a design preference but a kinematic necessity, derived directly from the terminal-guidance error budget. This geometry also converts the projectile cloud's aim-error response from a binary hit/miss criterion (as in a focused-beam system) to a continuous, non-catastrophic falloff in delivered momentum with miss distance.

### 2.4 Reentry and casualty risk

Projectile grains (50–1000 μm) and target debris up to 5 cm undergo complete thermal ablation on reentry under stagnation-point heating (Sutton–Graves correlation); expected ground-casualty probability is negligible against the 10⁻⁴-per-event regulatory threshold.

## 3. Economic constraint (the binding result)

Each of the four technical gates above (energy, momentum coupling, geometry, reentry safety) is cleared with margin — the concept is **technically feasible**. The binding limitation is orbital-mechanical, not physical.

Engagement rate — the frequency at which a single platform can re-orient to intercept a new target within its serviced debris cluster — is governed by a closed-form invariant:

$$\Delta v_{\text{per engagement}} = \frac{\pi v_{\text{orbital}}}{3.5\,\dot\Omega\, T_{\text{mission}}}$$

where $\dot\Omega$ is the target cluster's nodal regression rate and $T_{\text{mission}}$ is mission duration. This relation was validated to 1.14% against a full numerical orbital tour simulation. Evaluated at the most favorable debris cluster examined (Fengyun-1C fragmentation debris, sun-synchronous, ~865 km), it yields on the order of 16–55 engagements per platform-*decade*, independent of platform mass, material, or launcher design.

Since platform amortization dominates cost per object at these engagement rates (≈99% of total cost under the costed subsystem-mass-fraction platform design, $21.2M), cost per object of 1 cm debris is $0.9–4.2 million. A sensitivity analysis using an unsupported, optimistic platform-cost floor 42× below the costed design still yields $9,100–$31,300 per object — 1.6–82× more expensive than a matched-physics-derived cost for ground-based laser ablation ($381–$5,702 per object), the only competing non-rendezvous method addressing this size class.

Critically, the same $\Delta v$-invariant constrains a space-based laser-ablation implementation identically, since ablative recoil is directed along the line of sight and therefore requires the same co-orbital geometry. The economic finding is thus a structural property of the *architecture class* — targeted, co-orbital, single-object-per-engagement, non-rendezvous concepts — not a deficiency specific to FG01's material or hardware implementation.

## 4. Scope exception

For intermediate-mass targets (~100 kg – 3.7 t, e.g., derelict upper stages), the architecture is cost-*competitive* (not efficient) against robotic-tug rendezvous, because the projectile-based approach carries no dead mass beyond the projectile itself, whereas a tug must accelerate its own bus mass. This niche has not been engineered to platform level (magazine and propellant mass have not been re-budgeted into a vehicle design) and represents a distinct research question from the original sub-10 cm target class.

## 5. Principal conclusion

Technical feasibility and economic viability are governed by independent constraint sets for this class of concept. Under present technology, the decisive variable is achievable engagement cadence (an orbital-mechanics quantity), not projectile material, launcher type, or platform cost — each of which was tested as a sensitivity parameter and found non-determinative of the economic outcome.
