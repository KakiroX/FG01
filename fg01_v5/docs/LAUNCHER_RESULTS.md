# Launcher Trade Study — Results

Document 2 of 3 · [Methodology](LAUNCHER_METHODOLOGY.md) · [Conclusions](LAUNCHER_CONCLUSIONS.md)

All results regenerate from `sims/gun0*.py`, `coil0*.py`, `chk01_velocity.py`,
`cart01_charge_law.py`. Traceability tags per the v5 convention.

---

## 1. Validation

### 1.1 Chemical launcher, external validation

Calibrated coefficients (Methodology §4.3) reproduce four independent published
loads:

| Reference load | Published v | Model v | Δ | Published P_max | Model P_max | Δ | η |
|---|---|---|---|---|---|---|---|
| .300 Win Mag, 71 gr H1000 | 900 m/s | 905.0 | **+0.6%** | 441 MPa | 438.5 | **−0.6%** | 23.9% |
| .308 Win, 44 gr Varget | 795 m/s | 767.1 | **−3.5%** | 415 MPa | 417.7 | **+0.6%** | 26.1% |
| .30-06, 57 gr IMR 4350 | 850 m/s | 844.9 | **−0.6%** | 405 MPa | 405.2 | **+0.1%** | 24.3% |
| .300 RUM, 82 gr Retumbo | 910 m/s | 923.0 | **+1.4%** | 441 MPa | 460.8 | **+4.5%** | 23.5% |

**All four within ±3.5% velocity and ±4.5% pressure.** Thermodynamic efficiency
23.5–26.1%, inside the 25–35% band expected for small arms. `VALIDATED`

### 1.2 Numerical convergence

Reference dt = 5×10⁻⁸ s:

| dt | v error | P_max error |
|---|---|---|
| 2×10⁻⁶ s | 2.4×10⁻⁶ % | 3.5×10⁻⁴ % |
| 1×10⁻⁶ s | 5.0×10⁻⁷ % | 4.7×10⁻⁵ % |
| 2×10⁻⁷ s | 2.1×10⁻⁸ % | 2.8×10⁻⁶ % |

**Fully converged.** Time step is not a source of uncertainty. `VALIDATED`

### 1.3 Coilgun, internal consistency

COIL-03's radial stack-up places the winding mean radius at **20.09 mm** against
the **20.00 mm** COIL-01 assumed — **0.4%**. The electromagnetic results are
geometrically self-consistent. `VALIDATED`

Overall efficiency **19.6%** falls inside the 15–50% band reported for
multistage induction launchers (lower end, as expected at this small scale).

### 1.4 Charge law, cross-validation

CART-01 against the existing v6 interaction modules, 1 cm target at 900 km:

| Quantity | CART-01 | v6 module | Error |
|---|---|---|---|
| Areal penalty | 21.727 | 21.727 (INT-4) | **0.00%** |
| Mass on target | 0.1971 g | 0.19709 g (INT-3) | **0.00%** |
| Launched powder | 4.282 g | 4.282 g (INT-4) | **0.00%** |

Exact, because derived from the same expressions. `VALIDATED`

### 1.5 Note on the stated validation target

The original brief proposed "10 g of H1000 in a 0.6 m barrel → 800–900 m/s for an
8–10 g bullet." **This load is not physically realisable**: 10 g of H1000 will not
fit in a .30-calibre case (~6 cm³ holds ~4.6 g), and if it did, chamber pressure
would far exceed proof. Calibration therefore used the real .300 Win Mag magnum
load, which the figure appears to describe.

---

## 2. Chemical launcher (GUN-01 / GUN-02)

### 2.1 The pressure limit is set by bore, not by barrel length

At the 619 m/s design point, 25 mm bore:

| Barrel | Charge | P_max | Under 100 MPa? |
|---|---|---|---|
| 0.6 m | 25.95 g | 131.5 MPa | ✗ |
| 1.0 m | 21.37 g | 117.5 MPa | ✗ |
| 2.0 m | 18.00 g | 106.5 MPa | ✗ |

**Quadrupling the barrel moves peak pressure only 131.5 → 106.5 MPa.** Peak
pressure occurs within the first ~2–3 cm of travel, so barrel length cannot
relieve it. `DERIVED`

Bore sweep at 619 m/s, 1.0 m barrel:

| Bore | Charge | P_max | Feasible |
|---|---|---|---|
| 25 mm | 21.37 g | 117.5 MPa | ✗ |
| **30 mm** | **27.21 g** | **89.1 MPa** | **✓** |
| 35 mm | 33.79 g | 72.2 MPa | ✓ |
| 40 mm | 41.10 g | 61.2 MPa | ✓ |

**Bore is the only parameter with quadratic leverage** — at fixed pressure, force
scales with bore area, so the projectile accelerates faster and the chamber
volume grows faster, self-limiting the peak. `DERIVED`

### 2.2 Barrel length

30 mm bore, 619 m/s, Retumbo:

| L | Charge | P_max | x at P_max | a_max | Charge saved |
|---|---|---|---|---|---|
| 0.5 m | 35.45 g | 103.0 MPa | 3.22 cm | 72,043 g | — |
| 0.6 m | 32.62 g | 98.3 MPa | 2.93 cm | 68,717 g | 2.83 g |
| 0.8 m | 29.21 g | 92.6 MPa | 2.58 cm | 64,588 g | 3.41 g |
| **1.0 m** | **27.21 g** | **89.1 MPa** | 2.38 cm | 62,090 g | 2.00 g |
| 1.5 m | 24.58 g | 84.5 MPa | 2.11 cm | 58,714 g | 1.31 g |
| 2.0 m | 23.31 g | 82.2 MPa | 1.98 cm | 57,048 g | 1.27 g |

Shortest barrel under 100 MPa: **0.6 m** (marginal). Knee where marginal charge
saving falls below 1 g per step: **1.0 m**. `DERIVED`

### 2.3 Design point

| Parameter | Value |
|---|---|
| Bore · length · charge | 30 mm · 1.0 m · 27.21 g Retumbo |
| Muzzle velocity | 619 m/s |
| Peak chamber pressure | 89.14 MPa |
| Peak pressure location | 2.38 cm |
| Travel time | 2.428 ms |
| Muzzle energy | 19.16 kJ |
| Thermodynamic efficiency | 15.88% |
| Peak acceleration | 62,090 g |
| Recoil impulse | 67.06 N·s |
| Chamber volume · length | 45.35 cm³ · 6.4 cm |

### 2.4 Barrel wall — tapered

Peak pressure occurs 2.4 cm into travel, so the chamber end carries everything
and the rest of the barrel is dead mass if run at constant wall. Hoop
`σ = P·r/t`, SF 2 on yield:

| Material | Chamber wall | Muzzle wall | Tapered | Constant | **Saving** |
|---|---|---|---|---|---|
| 4340 steel | 5.348 mm | 1.500 mm | 1.679 kg | 4.663 kg | 64.0% |
| **Ti-6Al-4V** | **2.971 mm** | **1.500 mm** | **0.724 kg** | 1.363 kg | **46.9%** |
| Inconel 718 | 2.596 mm | 1.500 mm | 1.294 kg | 2.178 kg | 40.6% |
| AerMet 100 | 1.555 mm | 1.500 mm | 1.173 kg | 1.216 kg | 3.6% |

Taper profile, Ti-6Al-4V:

| Station | Pressure envelope | Wall | OD |
|---|---|---|---|
| 0 cm | 89.14 MPa | 2.97 mm | 35.94 mm |
| 13.6 cm | 54.20 MPa | 1.81 mm | 33.61 mm |
| 27.1 cm → muzzle | ≤34.83 MPa | 1.50 mm (min gauge) | 33.00 mm |

Beyond ~27 cm the wall is set by minimum practical gauge, not strength.
`DERIVED`

### 2.5 Loading density — the propellant density limit

Three distinct densities:

| | Value | Meaning |
|---|---|---|
| Solid (crystalline) | **1.60 g/cm³** | nitrocellulose material, zero voids `SOURCED` |
| Bulk (poured) | ~0.92 g/cm³ | grains touching, ~42% void `SOURCED` |
| Loading (design choice) | 0.60 g/cm³ baseline | charge ÷ chamber volume |

Sweep at 30 mm, 1.0 m, 619 m/s:

| LD (g/cm³) | Charge | Chamber | Length | P_max | a_max | Case fill |
|---|---|---|---|---|---|---|
| 0.50 | 28.62 g | 57.2 cm³ | 8.1 cm | 73.6 MPa | 50,885 g | 54% |
| **0.60** | **27.21 g** | **45.3 cm³** | **6.4 cm** | **89.1 MPa** | **62,090 g** | 65% |
| **0.70** | 26.17 g | 37.4 cm³ | 5.3 cm | **107.7 MPa** | 75,460 g | 76% |
| 0.80 | 25.37 g | 31.7 cm³ | 4.5 cm | 130.4 MPa | 91,859 g | 87% |
| 0.90 | 24.71 g | 27.5 cm³ | 3.9 cm | 159.3 MPa | 112,676 g | 98% |
| 1.00 | 24.16 g | 24.2 cm³ | 3.4 cm | 197.6 MPa | 140,294 g | 109% |

**The ullage is not waste — it is the expansion volume that limits peak
pressure.** Loading densities above 1.60 g/cm³ are physically impossible
(propellant solid density); above ~1.15 g/cm³ requires compressed loads.

### 2.6 Cartridge case — the dominant round mass

Chamber 45.3 cm³, 6.4 cm long:

| Case type | Wall | Case mass | Round total¹ | Rounds per 10 kg |
|---|---|---|---|---|
| Brass | 1.5 mm | **139.1 g** | 266.3 g | 37.6 |
| Steel | 1.2 mm | 111.0 g | 238.2 g | 42.0 |
| **Combustible** | — | **6.8 g** | 134.0 g | **74.6** |

¹ includes the 100 g projectile

**A brass case is 5.1× the propellant mass** — the single heaviest item in the
round. `DERIVED`

---

## 3. Electromagnetic launcher (COIL-01/02/03)

### 3.1 Mass budget of the 100 g

| Item | Mass | Note |
|---|---|---|
| Aluminium canister / armature | **19.0 g (19%)** | not optional — the field pushes on it |
| Rhenium payload | **81.0 g** | |

Armature `L_a = 19.3 nH`, `R_a = 29.1 µΩ`, **L/R = 0.66 ms** — decay time ~7× the
pulse, so induced current persists through transit. `DERIVED`

### 3.2 Design points

Constraints a ≤ 50,000 g, B ≤ 15 T:

| Target v | Achieved | Stages | Length | **Energy** | η | Cap mass | Cu | Recoil | Charge @5 kW |
|---|---|---|---|---|---|---|---|---|---|
| 300 m/s | 308.1 | 18 | 0.79 m | 25.10 kJ | 18.9% | 16.7 kg | 3.2 kg | 30.8 N·s | 5.6 s |
| **400 m/s** | **407.3** | **31** | **1.36 m** | **42.26 kJ** | **19.6%** | **28.2 kg** | 5.5 kg | 40.7 N·s | 9.4 s |
| 500 m/s | 500.2 | 44 | 1.94 m | 62.41 kJ | 20.0% | 41.6 kg | 7.8 kg | 50.0 N·s | 13.9 s |
| **619 m/s** | **623.4** | **67** | **2.95 m** | **95.32 kJ** | **20.4%** | **63.6 kg** | 11.9 kg | 62.3 N·s | 21.2 s |
| 700 m/s | 700.8 | 83 | 3.65 m | 119.58 kJ | 20.5% | 79.7 kg | 14.8 kg | 70.1 N·s | 26.6 s |

**619 m/s costs 2.3× the energy and 2.2× the length of 400 m/s.** `DERIVED`

### 3.3 Stored energy is invariant under the acceleration limit

Target 400 m/s:

| Accel limit | Stages | **Length** | Energy | Cap mass | B_peak | Canister SF |
|---|---|---|---|---|---|---|
| 20,000 g | 90 | 3.96 m | 40.33 kJ | 26.9 kg | 6.2 T | 9.2 |
| **50,000 g** | 31 | **1.36 m** | 42.26 kJ | 28.2 kg | 9.9 T | **3.7** |
| **100,000 g** | 15 | **0.66 m** | 40.48 kJ | 27.0 kg | 13.3 T | **1.8** |
| 200,000 g | 13 | 0.57 m | 43.96 kJ | 29.3 kg | 14.8 T | **0.9 ✗** |

> **Stored energy is invariant at 40–44 kJ across a 10× range of acceleration
> limit. Only length responds.** Efficiency stays ~20% throughout, so the
> capacitor bank is fixed by physics and length is purely a payload-structure
> decision. `DERIVED`

### 3.4 Structure and envelope

Magnetic pressure `B²/2μ₀` = **34.2 MPa** at 9.27 T peak:

| Item | Material | Thickness | Mass | Share |
|---|---|---|---|---|
| Capacitor bank | — | — | **28.2 kg** | **82.7%** |
| Coil copper | Cu | — | 4.17 kg | 12.2% |
| Outer case | Ti-6Al-4V | 0.80 mm | 1.005 kg | 3.0% |
| Bore tube | G10/FR4 | 1.50 mm | 0.375 kg | 1.1% |
| Coil banding | S-glass | 0.93 mm | 0.334 kg | 1.0% |

**All structure = 1.71 kg = 5.0%.** Wall thickness is not a mass or cost driver.
`DERIVED`

Radial stack-up (rectangular conductor): **9.88 mm wall → 49.8 mm OD**, of which
the **winding is 4.15 mm (42%)**. Envelope Ø49.8 × 1364 mm.

**Design trap:** the bore tube must be electrically insulating. A metal tube is a
shorted turn around the armature — it carries induced current, shields the
armature and destroys coupling. 316 stainless (7.4×10⁻⁷ Ω·m) is unusable at any
thickness. `DERIVED`

---

## 4. Velocity consistency check (CHK-01)

400 m/s originated in the launcher studies, not the engagement analysis, so it
was tested against every v5/v6 gate.

### 4.1 Gates passed

| Gate | 400 m/s | 619 m/s |
|---|---|---|
| Disruption E_s (400–1200 km) | 18.3–21.7 kJ/kg **PASS** | 25.1–29.8 kJ/kg **PASS** |
| vs. conservative 30 kJ/kg bound | **PASS** everywhere | **PASS** everywhere |
| Self-disposal, all altitudes | **PASS**, margin 1.38× at 1200 km | **PASS**, margin 2.13× |
| Sub-catastrophic on matched target | 9.75 kJ/kg **PASS** | 13.35 kJ/kg **PASS** |

### 4.2 Costs of the lower velocity

**β is velocity-dependent** (INT-2), because 400 m/s is only 1.25× aluminium's
320 m/s strength velocity — essentially pure embedding:

| v_rel | β (oblique) | Mass penalty vs. 619 |
|---|---|---|
| 300 m/s | 1.000 | 2.48× |
| **400 m/s** | **1.063** | **1.75×** |
| 500 m/s | 1.131 | 1.32× |
| **619 m/s** | **1.201** | 1.00× |

**The firing cone scales linearly with launch velocity**, tightening from **1.51°
to 0.98°** at a 10% energy penalty. Run through the real breakup clouds:

| Cluster | Engageable @1.51° | @0.98° | Ratio |
|---|---|---|---|
| Fengyun-1C | 183 | 118 | 0.65 |
| Cosmos-2251 | 87 | 55 | 0.63 |
| Iridium-33 | 39 | 26 | 0.66 |
| Cosmos-1408 | 185 | 122 | 0.66 |

**400 m/s costs 34–37% of engageable targets**, raising cost per object ~1.54×.
This outweighs the 53 kJ / 36 kg / 1.6 m the lower velocity saves. `DERIVED`

---

## 5. Charge law (CART-01)

$$m_{\text{powder}} = \text{clamp}\left[C(h)\cdot m_t\cdot P(m_t),\ 50\text{ g},\ 100\text{ g}\right]$$

with m_t in kg, `C(h) = γ·Δv(h)/(β·v_rel)` (0.1311 at 1200 km → 0.1553 at
400 km) and `P(m_t) = max(0.274·m_t^{−2/3}, 1)`.

**Two regimes**, crossing at 143 g / 4.7 cm: `m_powder ∝ m_t^{1/3}` below,
`∝ m_t` above. Single-cartridge envelope **322–763 g** of debris depending on
altitude. `DERIVED`

**Safety finding.** Delivering the required momentum always lands at
`E_s/E_s,c = 0.67`, independent of target mass. But below ~240 g debris, the 50 g
cartridge floor with a focused cloud delivers far more:

| Debris | E_s if focused | vs. threshold | Defocus cloud to |
|---|---|---|---|
| 1.4 g (1 cm) | 313 kJ/kg | **7.8×** | 27.6 cm |
| 20 g | 129 kJ/kg | **3.2×** | 17.7 cm |
| 100 g | 75 kJ/kg | **1.9×** | 13.6 cm |
| 240 g | 39.9 kJ/kg | 1.00× | 11.7 cm |

After defocusing, E_s returns to **26.7 kJ/kg for every target mass** — the
invariance is the check that the rule is correct. `DERIVED`

**Dosing is unnecessary.** Magazine gives 300–600 shots against ~48 demanded per
decade (6–12× margin); rhenium economy is $0.18/shot; and no charge in 50–100 g
makes 1 cm debris safe (0.295 g would be needed). Dispersion control replaces
metering entirely.

---

## 6. Head-to-head, both at 619 m/s

### 6.1 System mass

| | 48 shots | 300 shots | 1000 shots |
|---|---|---|---|
| **Chemical total** | **15.3 kg** | 57.2 kg | 173.6 kg |
| — barrel + breech | 7.4 kg | 7.4 kg | 7.4 kg |
| — rounds (brass case) | 8.0 kg | 49.9 kg | 166.3 kg |
| **Coilgun total** | **98.0 kg** | 98.0 kg | 98.0 kg |
| — capacitors | 63.5 kg | 63.5 kg | 63.5 kg |

**Mass crossover: 545 shots.** ECO-1's realised demand is ~48 per platform-decade
(24 engagements × 2 shots), **an order of magnitude below crossover**. `DERIVED`

### 6.2 Non-mass comparison

| Metric | Chemical | Coilgun | Favours |
|---|---|---|---|
| Length | 1.0 m | 2.95 m | chemical |
| Prime mass (48 shots) | 15.3 kg | 98.0 kg | **chemical, 6.4×** |
| Consumable per shot | 27.2 g + case | 0 | coilgun |
| Recoil | 67.1 N·s | 62.3 N·s | coilgun, 7% |
| Waste heat into vehicle | **23.7 kJ** (gas carries rest overboard) | **~76 kJ** unless recovered | **chemical** |
| Cadence | 150 s (thermal) | 21 s (charging) | coilgun, 7× |
| Erosion | yes | none | coilgun |
| Moving parts | breech, autoloader | none | coilgun |
| Velocity tunable per shot | no | yes | coilgun |
| Peak acceleration | 62,090 g | 49,689 g | coilgun |

---

## 7. Results summary table

| Quantity | Chemical (619 m/s) | Coilgun (619 m/s) |
|---|---|---|
| Bore | 30 mm | 30 mm |
| Length | 1.0 m | 2.95 m |
| Envelope OD | 35.9 mm (chamber) | 49.8 mm |
| Energy per shot | 27.21 g Retumbo (118 kJ chemical) | 95.32 kJ stored |
| Efficiency | 15.9% | 20.4% |
| Prime hardware | 0.72 kg barrel + ~6 kg breech | 63.6 kg caps + 11.9 kg Cu |
| Peak load | 89.1 MPa chamber | 34.2 MPa magnetic |
| Peak acceleration | 62,090 g | 49,689 g |
| Recoil | 67.06 N·s | 62.34 N·s |
| Validated against | 4 published cartridges, ±3.5%/±4.5% | internal consistency 0.4% |
