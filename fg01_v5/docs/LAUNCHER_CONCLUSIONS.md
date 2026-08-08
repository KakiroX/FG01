# Launcher Trade Study — Conclusions and Recommendation

Document 3 of 3 · [Methodology](LAUNCHER_METHODOLOGY.md) · [Results](LAUNCHER_RESULTS.md)

---

## 1. Recommendation

> **Adopt the chemical (gunpowder) launcher: 30 mm smoothbore, 1.0 m barrel,
> 27.2 g Retumbo-class propellant, 100 g projectile at 619 m/s, 89 MPa peak
> chamber pressure, Ti-6Al-4V tapered wall, combustible cases.**
>
> **Exclude the railgun** on scoping grounds (§4). **Retain the coilgun as the
> backup architecture**, selected only if shot demand exceeds ~500 per platform
> or if consumables handling proves unacceptable.

Confidence: **moderate-to-high** on the mass argument, **moderate** on the
absolute chemical numbers (validated to ±3.5% on velocity, but only against
small-arms cartridges — see §6).

---

## 2. The decision in one number

**Mass crossover: 545 shots.**

Below it the chemical launcher is lighter; above it the coilgun is. ECO-1 puts
realised demand at **~48 shots per platform-decade** (24 engagements × 2 shots) —
**an order of magnitude below crossover**. At that demand:

| | Chemical | Coilgun | Ratio |
|---|---|---|---|
| System mass, 48 shots | **15.3 kg** | 98.0 kg | **6.4×** |
| Length | **1.0 m** | 2.95 m | 2.95× |
| Prime cost driver | 0.72 kg barrel | 63.6 kg capacitors | — |

The coilgun's advantage is that it carries **no consumable**, so its mass is flat
in shot count. That is decisive for a magazine-limited weapon fired thousands of
times. FG01 is not that: it is a low-cadence, few-dozen-shot platform, and it
pays 83 kg up front for a benefit it never collects.

---

## 3. Why the earlier verdict is reversed

**[COILGUN_DESIGN.md](COILGUN_DESIGN.md) recommended the coilgun. This study
reverses that.** The reversal is recorded here rather than by editing the
earlier document, so the reasoning chain stays auditable.

Three of the four arguments that favoured the coilgun have since failed:

| Original argument | Status now | Why |
|---|---|---|
| **Per-shot velocity tunability** | **withdrawn** | CHK-01 shows 400 m/s costs 34–37% of engageable targets, so there is no useful lower setting; CART-01 shows dosing is unnecessary and dispersion control replaces metering. There is nothing left to tune. |
| **Thermal advantage** | **reversed — was stated backwards** | The chemical gun dumps ~23.7 kJ into the barrel and vents the remainder overboard with the propellant gas. The coilgun must **dissipate ~76 kJ inside the vehicle** (80% of 95 kJ stored) with no working fluid to carry it away. Thermal management favours the chemical gun, not the coilgun. |
| **No erosion, no moving parts** | **stands** | Genuine and unchanged. This is the coilgun's real advantage. |
| **Mass competitiveness** | **failed at the design velocity** | The coilgun was originally sized at 400 m/s (42 kJ). At the 619 m/s CHK-01 requires it needs 95 kJ, 67 stages, 2.95 m and 63.6 kg of capacitors. |

The single largest driver of the reversal is that **CHK-01 moved the design
velocity from 400 to 619 m/s.** Coilgun energy scales as v², so this is a 2.3×
penalty on the dominant mass item; the chemical gun absorbs the same change with
6 g more powder and a larger bore.

---

## 4. Railgun — excluded, not modelled

The brief asked for a railgun. **No railgun simulation was written**, and this
document does not claim one. The exclusion is a scoping decision made from
closed-form estimates:

Using `F = ½L′I²` with a representative rail inductance gradient
L′ ≈ 0.4–0.5 µH/m, delivering 19.2 kJ of muzzle energy to a 100 g projectile
requires **451–637 kA** through a sliding contact.

| | FG01 requirement | US Navy EMRG (reference) |
|---|---|---|
| Muzzle energy | 19.2 kJ | 32 MJ |
| Rail current | 451–637 kA | 3–5 MA |
| Bore life | needs hundreds of shots | tens of shots before rail refurbishment |
| Efficiency | — | ~20–25% at the muzzle, lower at the wall plug |

The disqualifying issue is not the current magnitude but that a railgun draws it
through a **sliding physical contact**. Rail erosion, armature transition to arc
mode, and pulsed-power switching at that current are all unsolved for a
long-life, unattended, on-orbit system, and none of it buys anything a coilgun
does not — the coilgun reaches the same velocity with no contact at all. **A
railgun is strictly dominated by the coilgun for this application**, so it was
excluded before modelling rather than after.

This is a scoping judgement, tagged `UNVALIDATED`. If the railgun must be ruled
out formally, that requires its own study.

---

## 5. What the chemical launcher must accept

The recommendation is not free. Four honest costs:

**1. Peak acceleration 62,090 g.** Higher than the coilgun's 49,689 g, and it is
now the binding structural constraint on the rhenium canister. Loading density
0.60 g/cm³ was chosen partly to hold this down — going to 0.90 g/cm³ (a nearly
full case) would drive it to 112,676 g at 159 MPa, which the canister does not
survive at reasonable safety factor.

**2. Consumables handling.** 27.2 g of propellant per round plus a case, an
autoloader, and a live-propellant magazine on an unattended spacecraft. The
coilgun has none of this. Combustible cases (6.8 g vs 139 g brass) are strongly
recommended — they are 5.1× lighter than brass and remove case ejection
entirely — but they are the **least mature** item in the design.

**3. Barrel erosion.** Present and cumulative. A 1 mm Inconel 718 liner over the
first ~15 cm is specified; life is not characterised. `UNVALIDATED`

**4. Cadence.** ~150 s between shots, thermally limited, against the coilgun's
21 s charge time. Irrelevant at FG01's engagement rate, but it forecloses any
future rapid-fire concept of operations.

---

## 6. Confidence and what would change the answer

| Claim | Confidence | Basis |
|---|---|---|
| Chemical is lighter below ~500 shots | **High** | Robust to ±50% on every input; the gap is 6.4× |
| 619 m/s is the right design velocity | **High** | CHK-01, four independent breakup clouds |
| 30 mm bore is the minimum at ≤100 MPa | **High** | `DERIVED`, converged, validated model |
| 27.2 g charge, 89.1 MPa | **Moderate** | Model validated only against 7.62 mm rifle cartridges; 30 mm is a 4× bore extrapolation |
| 62,090 g peak acceleration | **Moderate** | Same extrapolation; also the binding canister constraint |
| Coilgun 95.3 kJ at 619 m/s | **Moderate** | Internally consistent (0.4%) but no external validation; DC resistance assumed, so 20.4% efficiency is optimistic |
| Railgun exclusion | **Low–moderate** | Scoping estimate only, `UNVALIDATED` |

**The recommendation would flip if:**

- shot demand rose above ~500 per platform (crossover), or
- the rhenium canister could not be qualified to 62,000 g while surviving 50,000 g
  (removing the chemical gun's structural headroom), or
- live propellant on an unattended long-duration platform were ruled
  inadmissible on safety grounds — a programmatic constraint this study does not
  evaluate, and the most likely single reason to choose the coilgun anyway.

**The recommendation would NOT flip if** the coilgun's efficiency, capacitor
energy density or wall thicknesses improved substantially: capacitors are 82.7%
of its mass, and even a 2× energy-density improvement leaves it at ~66 kg
against 15.3 kg.

---

## 7. Open items

| # | Item | Owner discipline |
|---|---|---|
| 1 | Re-validate the ballistics model against a **large-bore** (≥30 mm) reference gun | interior ballistics |
| 2 | Qualify the rhenium canister to 62,000 g | structures |
| 3 | Combustible-case maturity and long-duration propellant storage | propulsion / materials |
| 4 | Barrel erosion life with the Inconel liner | materials |
| 5 | Recoil management, 67 N·s per shot, and its attitude-control cost | GNC |
| 6 | Formal railgun exclusion study, if required for the record | EM launch |
| 7 | Re-run COIL-01 with AC resistance before quoting 20.4% | electrical |

---

## 8. Statement of limitations

This trade study is a **0-D lumped-parameter comparison**. It does not model:
two-phase propellant flow, heat transfer to the bore, projectile balloting or
in-bore yaw, sabot discard, coil AC losses or proximity effect, switch losses,
armature heating, or any failure mode. It compares two architectures on mass,
envelope, energy and peak load at a fixed velocity requirement, and nothing
else. The conclusions are sound at that level of fidelity and should not be used
to size flight hardware without the open items in §7 closed.

All simulations here are **additive**: they import the v5/v6 library read-only,
write only their own `gun0*`/`coil0*`/`chk01`/`cart01` outputs, and are not
registered in `run_all.py`. No existing simulation or result was modified.
