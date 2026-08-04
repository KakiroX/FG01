# Charge Dosing — Decision Memo

**Question:** does the rhenium charge need to be adjustable in flight, or is a
fixed charge sufficient?

**Answer: fix the charge at 81 g. Make the *dispersion* adjustable instead.**

Source: [`sims/cart01_charge_law.py`](../sims/cart01_charge_law.py) ·
detail in [CARTRIDGE_CHARGE_LAW.md](CARTRIDGE_CHARGE_LAW.md)

---

## The options

| scheme | powder | total accel¹ | fits 100 g? | capacitor | shots / 30 kg | vs 48 needed² |
|---|---|---|---|---|---|---|
| constant 100 g | 100 g | 119 g | **no** | 48.6 kJ / 32.4 kg | 300 | 6.2× |
| **constant 81 g** | **81 g** | **100 g** | **yes** | **40.8 kJ / 27.2 kg** | **370** | **7.7×** |
| dosed 50–100 g | 75 g avg | 94 g | yes | 38.4 kJ / 25.6 kg | 400 | 8.3× |

¹ powder + 19 g aluminium canister/armature · ² ECO-1: ~24 engagements per
platform-decade, 2 shots each

---

## Why metering isn't worth a mechanism

**1. It buys shots that aren't needed.** A 30 kg magazine yields 300–600 shots
against a demand of ~48 per decade — 6–12× margin either way. The magazine has
never been the binding constraint; the Δv-per-engagement invariant is.

**2. Rhenium economy is trivial.** 100 g versus a 75 g dosed average saves
25 g/shot = **$0.18/shot ≈ $9 per decade**, against a platform costing $50M+.

**3. Dosing cannot solve the problem it appears to solve.** Engaging 1 cm debris
safely requires a **0.295 g** charge — 170× below the 50 g floor. No setting in
the 50–100 g range makes small debris safe. Only defocusing does.

---

## Why 81 g and not 100 g

100 g of powder plus the 19 g canister is **119 g accelerated — 19% over** the
100 g the coilgun is sized for, requiring **48.6 kJ / 32.4 kg** of capacitors
instead of 40.8 kJ / 27.2 kg.

> A constant 100 g charge pays **~5 kg of capacitor mass permanently** to avoid a
> mechanism it does not need anyway.

**81 g + 19 g = exactly 100 g**, matching the existing COIL-01 design point
(42 kJ, 28 kg), giving 370 shots and single-shot coverage to **581 g** of debris.
581 g – 1.16 kg takes two shots.

**Note on scale:** this 388–581 g single-shot envelope is a narrow niche in its
own right (6–8 cm objects) and should not be conflated with ECO-5's separately
identified 100 kg–3.7 t intermediate-mass window — the two are two to three
orders of magnitude apart and do not overlap.

---

## What is genuinely required: dispersion control

Delivered momentum is `charge × ε(cloud width)`, so **either** knob can size a
shot — but only dispersion reaches the sub-gram incident masses small debris
needs. Firing a fixed charge at a small target with a focused cloud is
catastrophic:

| debris | E_s if focused (81 g) | vs threshold | defocus cloud to |
|---|---|---|---|
| 1.4 g (1 cm) | 313 kJ/kg | **7.8×** | **27.6 cm** |
| 20 g | 129 kJ/kg | **3.2×** | **17.7 cm** |
| 100 g | 75 kJ/kg | **1.9×** | **13.6 cm** |
| 388 g | 40 kJ/kg | 1.00× | — (safe focused) |

$$d_{\text{cloud}} = d_t\sqrt{2.995\,\frac{m_{\text{powder}}}{m_{\text{incident,req}}}}
\qquad m_{\text{incident,req}} = \frac{\gamma\,m_t\,\Delta v(h)}{\beta\,v_{rel}}$$

After defocusing, **E_s = 26.7 kJ/kg (0.67× threshold) for every target mass** —
this invariance is the check that the rule is correct. Surplus powder
self-disposes (SIM-2: retrograde launch above Δv_direct reenters within one
orbit), so the only cost is wasted rhenium at ~$7.28/kg.

**If dispersion is not controllable, do not engage debris below ~388 g** — the
shot would be net debris-generating regardless of charge setting.

---

## Recommendation

| Parameter | Setting |
|---|---|
| Rhenium charge | **fixed 81 g** |
| Total accelerated mass | 100 g |
| Capacitor bank | 42 kJ (28 kg) — unchanged from COIL-01 |
| Adjustable parameter | **cloud dispersion, ~5–28 cm at 10 m range** |
| Single-shot envelope | debris 388 g – 581 g focused; any smaller mass via defocus |
| Metering mechanism | **not required** |

**Caveats.** The defocus branch needs ~2.5 m/s of transverse dispersion at 10 m
range to reach 28 cm — not demonstrated hardware. β = 1.2007 is derived, not
measured, and enters the charge inversely. σ_miss = 5.775 mm is assumed; charge
scales as σ² in the small-target regime.
