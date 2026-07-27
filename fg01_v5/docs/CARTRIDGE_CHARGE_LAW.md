# CART-01 — Cartridge Charge Law

**How many grams of rhenium powder to load, given debris mass and orbit.**

Code: [`sims/cart01_charge_law.py`](../sims/cart01_charge_law.py) · Outputs: `outputs/cart01_*`

> **Additive only.** This module imports from `fg01/` read-only, writes only
> `cart01_*` outputs, and is **not** registered in `run_all.py`. No existing
> simulation, output or document is modified.

---

## 1. The formula

$$m_{\text{powder}} = \text{clamp}\!\left[\; C(h)\cdot m_t \cdot P(m_t),\;\; 50\ \text{g},\;100\ \text{g}\right]$$

with **m_t in kg** and

$$C(h)=\frac{\gamma\,\Delta v(h)}{\beta\, v_{rel}}\qquad
P(m_t)=\max\!\left(0.274\, m_t^{-2/3},\; 1\right)$$

| Symbol | Meaning | Value used |
|---|---|---|
| γ | engineering margin | 2.0 |
| β | momentum enhancement (INT-2) | 1.2007 at 619 m/s |
| v_rel | intercept velocity | 619 m/s |
| Δv(h) | orbit-lowering Δv for 200 km | 48.7–57.7 m/s (table below) |
| P | areal-capture penalty (INT-4) | `max(K/A_t, 1)`, K = 1.706×10⁻³ m² |

**C(h) — the only altitude dependence:**

| h (km) | 400 | 500 | 600 | 700 | 800 | 900 | 1000 | 1100 | 1200 |
|---|---|---|---|---|---|---|---|---|---|
| Δv (m/s) | 57.72 | 56.45 | 55.23 | 54.05 | 52.91 | 51.81 | 50.74 | 49.72 | 48.72 |
| **C(h)** | **0.1553** | 0.1519 | 0.1486 | 0.1454 | 0.1424 | **0.1394** | 0.1365 | 0.1338 | **0.1311** |

Altitude is a **weak** lever — only 18% across the whole LEO band. Debris mass
dominates completely.

### Two regimes

$$m_{\text{powder}} \;\propto\; m_t^{1/3}\ \ (\text{small debris, penalty dominates})
\qquad m_{\text{powder}} \;\propto\; m_t\ \ (\text{large debris, penalty} = 1)$$

The crossover is at **A_t = K**, i.e. **d = 4.7 cm, m_t = 143 g**. Below it the
target is smaller than the shot pattern and most powder misses; above it
essentially all of it lands.

---

## 2. Lookup table — grams of powder to load

| debris mass | 400 km | 600 km | 800 km | 900 km | 1000 km | 1200 km |
|---|---|---|---|---|---|---|
| 1.4 g (1 cm) | 50 ⚠ | 50 ⚠ | 50 ⚠ | 50 ⚠ | 50 ⚠ | 50 ⚠ |
| 10 g | 50 ⚠ | 50 ⚠ | 50 ⚠ | 50 ⚠ | 50 ⚠ | 50 ⚠ |
| 50 g | 50 ⚠ | 50 ⚠ | 50 ⚠ | 50 ⚠ | 50 ⚠ | 50 ⚠ |
| 100 g | 50 ⚠ | 50 ⚠ | 50 ⚠ | 50 ⚠ | 50 ⚠ | 50 ⚠ |
| 200 g | 50 ⚠ | 50 ⚠ | 50 ⚠ | 50 ⚠ | 50 ⚠ | 50 ⚠ |
| 300 g | 50 ⚠ | 50 ⚠ | 50 ⚠ | 50 ⚠ | 50 ⚠ | 50 ⚠ |
| **400 g** | **62.1** | **59.4** | **56.9** | **55.8** | **54.6** | **52.4** |
| **500 g** | **77.7** | **74.3** | **71.2** | **69.7** | **68.3** | **65.6** |
| **600 g** | **93.2** | **89.2** | **85.4** | **83.6** | **81.9** | **78.7** |
| 700 g | 100 (×2) | 100 (×2) | 99.7 | 97.6 | 95.6 | 91.8 |
| 800 g | 100 (×2) | 100 (×2) | 100 (×2) | 100 (×2) | 100 (×2) | 100 (×2) |
| 1000 g | 100 (×2) | 100 (×2) | 100 (×2) | 100 (×2) | 100 (×2) | 100 (×2) |

⚠ = charge floored at 50 g; **the cloud must be defocused — see §4**.
(×2) = two shots needed.

---

## 3. Single-cartridge envelope

The mass range one cartridge services **without** flooring or capping:

| h (km) | min debris | max debris | diameter range |
|---|---|---|---|
| 400 | 322 g | 644 g | 6.11 – 7.69 cm |
| 600 | 336 g | 673 g | 6.20 – 7.81 cm |
| **900** | **359 g** | **717 g** | **6.33 – 7.98 cm** |
| 1200 | 381 g | 763 g | 6.46 – 8.14 cm |

**The 50–100 g cartridge is a 2:1 adjustment, so it maps onto a 2:1 debris-mass
band: roughly 320–760 g, i.e. 6–8 cm objects.** That is a narrow niche, and it
sits neatly inside the 100 kg–3.7 t window ECO-5 identified as the only regime
where the momentum kick beats a tug — at the small end of it.

---

## 4. ⚠ The 50 g floor is a hazard for small debris

Delivering exactly the required momentum always lands at the same place:

$$\frac{E_s}{E_{s,c}}=\frac{\gamma\,\Delta v\, v_{rel}}{2\beta E_{s,c}}=0.67
\quad\text{— independent of target mass}$$

So a *correctly sized* shot is always sub-catastrophic with the same 33% margin.
**The danger is being unable to size it down.** If you fire the 50 g minimum at a
small target with a normally-focused cloud, you deliver far more momentum than
the target can absorb:

| debris | if fired focused → E_s | vs threshold | catastrophic? | required incident | **defocus cloud to** |
|---|---|---|---|---|---|
| 1.4 g (1 cm) | 313 kJ/kg | **7.8×** | **YES** | 0.20 g | **27.6 cm** |
| 5 g | 205 kJ/kg | **5.1×** | **YES** | 0.70 g | **22.3 cm** |
| 20 g | 129 kJ/kg | **3.2×** | **YES** | 2.79 g | **17.7 cm** |
| 50 g | 95 kJ/kg | **2.4×** | **YES** | 6.97 g | **15.2 cm** |
| 100 g | 75 kJ/kg | **1.9×** | **YES** | 13.9 g | **13.6 cm** |
| 150 g | 64 kJ/kg | **1.6×** | **YES** | 20.9 g | **12.7 cm** |
| 240 g | 39.9 kJ/kg | **1.00×** | borderline | 33.5 g | 11.7 cm |
| 360 g | 26.6 kJ/kg | 0.67× | no | 50.2 g | 4.9 cm |

**Rule: below ~240 g debris, a focused 50 g shot shatters the target.**

### The fix — defocus, don't refuse

Widen the dispersion so only the required mass lands and the surplus misses.
After defocusing, E_s returns to **26.7 kJ/kg (0.67× threshold) for every target
mass** — the last column of the table above is constant, which is the check that
the rule is right.

Required cloud diameter to make the 50 g floor safe:

$$d_{\text{cloud}} = d_t\sqrt{\frac{2.995\, m_{\text{powder}}}{m_{\text{incident,req}}}}$$

This is a genuinely useful capability rather than a workaround: the surplus
powder still self-disposes (SIM-2 — retrograde launch above Δv_direct puts it on
an immediate-reentry trajectory), so the only cost is wasted rhenium at ~$0.36/g.

**Operationally:** the tuning knob you need is *dispersion*, not just charge
mass. If the hardware can only vary charge and not cloud width, then **do not
engage debris below ~240 g** — the shot would be net debris-generating.

---

## 4b. Do you need to dose at all? — **No. Use a constant 81 g.**

Metering is a mechanism, and mechanisms fail. Here is what it actually buys:

| scheme | charge | total accel | fits 100 g? | cap energy | cap mass | shots / 30 kg | vs 48 needed | 1-shot target max | safe w/o defocus above |
|---|---|---|---|---|---|---|---|---|---|
| constant 100 g | 100 g | 119 g | **no** | 48.6 kJ | 32.4 kg | 300 | 6.2× | 717 g | 479 g |
| **constant 81 g** | **81 g** | **100 g** | **yes** | **40.8 kJ** | **27.2 kg** | **370** | **7.7×** | **581 g** | **388 g** |
| constant 50 g | 50 g | 69 g | yes | 28.2 kJ | 18.8 kg | 600 | 12.5× | 359 g | 240 g |
| dosed 50–100 g | 75 g avg | 94 g | yes | 38.4 kJ | 25.6 kg | 400 | 8.3× | 538 g | 359 g |

**Three reasons dosing is not worth a mechanism:**

1. **It buys shots you don't need.** A 30 kg magazine gives 300–600 shots. ECO-1's
   realised demand is ~24 engagements per decade, ~48 shots at two per
   engagement — a **6–12× margin either way**. The magazine has never been the
   binding constraint.
2. **Rhenium economy is negligible.** 100 g vs a 75 g dosed average is 25 g/shot
   = **$0.18/shot**, about **$9 per decade**. Against a platform costing $50M+.
3. **Dosing cannot solve the problem it looks like it solves.** To engage 1 cm
   debris safely you would need a **0.295 g** charge — 170× below your 50 g
   floor. *No* setting in the 50–100 g range makes small debris safe. Only
   defocusing does.

**So the real trade is charge metering versus dispersion control — and you need
dispersion control regardless.** Once you have it, metering is redundant: cloud
width alone can size any shot from any fixed charge.

### Why 81 g and not 100 g

100 g of powder plus the 19 g canister is **119 g accelerated**, 19% over the
100 g the coilgun is sized for. That costs **48.6 kJ / 32.4 kg of capacitors**
instead of 40.8 kJ / 27.2 kg — **you pay ~5 kg of capacitor mass permanently, to
avoid a mechanism you don't need anyway.**

**81 g of powder + 19 g canister = exactly 100 g**, which fits the existing
COIL-01 design point (42 kJ, 28 kg), gives 370 shots, and still covers
single-shot targets up to 581 g. Targets from 581 g to 1.16 kg take two shots.

**Recommendation: fix the charge at 81 g. Make the dispersion adjustable
instead.**

---

## 5. Note on total accelerated mass

The 50–100 g is *powder*. Adding the 19 g aluminium canister/armature (COIL-01):

| | powder | + canister | total accelerated |
|---|---|---|---|
| minimum charge | 50 g | 19 g | **69 g** |
| maximum charge | 100 g | 19 g | **119 g** |

**COIL-01 was sized for 100 g *total*.** A full 100 g powder charge is therefore
**19% over** that design point and needs **≈50 kJ rather than 42 kJ** — muzzle
energy scales linearly with mass at fixed velocity. Either:

- size the capacitor bank for **50 kJ** (34 kg of capacitors) to cover the full
  charge range, or
- cap the powder at **81 g** to stay inside the 100 g total, or
- accept ~365 m/s instead of 400 m/s at full charge on a 42 kJ bank.

---

## 6. Validation

The charge law is derived from the same expressions as the existing modules and
reproduces them **exactly**, which is why it can be trusted without re-running
anything:

| Quantity, 1 cm target at 900 km | CART-01 | existing | error |
|---|---|---|---|
| areal penalty | 21.727 | 21.727 (INT-4) | **0.00%** |
| mass on target | 0.1971 g | 0.19709 g (INT-3) | **0.00%** |
| ideal launched powder | 4.282 g | 4.282 g (INT-4) | **0.00%** |

## 7. Limitations

- **Debris modelled as a solid Al-6061 sphere** (AM-2). Real fragments are
  irregular with larger silhouette per unit mass, which *lowers* the areal
  penalty and the required charge — so the law is conservative for real debris.
  If a measured diameter is available, use A_t directly instead of inferring it
  from mass.
- **β = 1.2007 assumed** at 619 m/s (INT-2, derived not measured). β enters
  inversely: a 10% error in β is a 10% error in charge.
- **σ_miss = 5.775 mm assumed** (the SIM-4 budget). Charge scales as σ², so this
  is the most sensitive input in the small-target regime — a 2× worse aim
  quadruples the required powder.
- **Fixed v_rel = 619 m/s.** CHK-01 showed 400 m/s passes the gates but tightens
  the firing cone by 35%; if the operating velocity changes, both β and C(h)
  change and the law must be re-evaluated (the code takes `v_rel` as an argument).
- **The defocus branch assumes cloud width is controllable** to ~40 cm at 10 m
  range, i.e. a transverse dispersion of ~2.5 m/s. Not demonstrated hardware.
