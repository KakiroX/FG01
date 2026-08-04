# Coilgun Wall Thickness — Cost-Model Reference

Design point: COIL-01, 30 mm bore, 100 g at 400 m/s, 31 stages, 1.36 m.
Sources: [`coil02_structure.py`](../sims/coil02_structure.py) (walls) ·
[`coil03_stackup.py`](../sims/coil03_stackup.py) (full radial build)

---

## Total tube dimensions — the number for stowage

| | |
|---|---|
| Bore | **30 mm** |
| **Total radial wall** | **10.9 mm** |
| **Outer diameter** | **51.8 mm** |
| Length | 1364 mm |
| Tube mass (all layers, 31 stages) | 8.9 kg |
| Swept volume | 2.9 litres |

### Radial stack-up, bore to outside

| Layer | Thickness | r_outer | Mass/stage |
|---|---|---|---|
| running clearance | 0.50 mm | 15.50 mm | — |
| bore tube (G10) | 1.50 mm | 17.00 mm | 11.3 g |
| former / ground insulation | 0.50 mm | 17.50 mm | 3.0 g |
| **WINDING** (2 layers × 15 turns, AWG 10) | **5.18 mm** | 22.68 mm | **234.1 g** |
| ground insulation | 0.50 mm | 23.18 mm | 4.0 g |
| banding (S-glass) | 0.93 mm | 24.11 mm | 11.1 g |
| thermal gap | 1.00 mm | 25.11 mm | — |
| outer case (Ti) | 0.80 mm | 25.91 mm | 22.7 g |

**The winding is 47% of the wall and 82% of the tube mass** — 22 turns of AWG 10
will not fit in one layer of a 40 mm coil (only 15 fit), so it needs two layers,
each a full wire diameter of radial build. The *structural* walls are only
5.7 mm of the 10.9 mm.

> **Consistency check:** this stack puts the winding mean radius at **20.09 mm**
> against the **20.00 mm** COIL-01 assumed — 0.4%. The electromagnetics stand.

### If you need a smaller tube

| Conductor | Area | Layers | Wall | **OD** | Cu (31 stages) |
|---|---|---|---|---|---|
| AWG 8 | 8.37 mm² | 2 | 12.3 mm | 54.5 mm | 6.7 kg |
| **AWG 10** (baseline) | 5.26 mm² | 2 | 10.9 mm | **51.8 mm** | 4.1 kg |
| **3×2 mm rectangular** | **6.00 mm²** | 2 | 9.7 mm | **49.5 mm** | 4.5 kg |
| AWG 12 | 3.31 mm² | 2 | 9.8 mm | 49.7 mm | 2.5 kg |
| AWG 14 | 2.08 mm² | **1** | 7.4 mm | **44.7 mm** | 1.5 kg |

**Rectangular conductor is the free win:** 3×2 mm gives a *smaller* tube (49.5 mm)
at *larger* conductor area (6.0 mm²) than AWG 10 — so lower resistance too.

Finer round wire shrinks the tube more (AWG 14 fits a single layer, 44.7 mm OD)
but **conductor area falls 2.5×, so coil resistance rises 2.5×**. COIL-01's
19.6% efficiency was computed at AWG 10; going finer requires re-running it.

---

## The rectangular-conductor build

**Conductor: two parallel strips of 3 mm (axial) × 1 mm (radial), wound
bifilar** — 6.0 mm² total, same as one 3×2 mm bar but split for skin depth
(see below).

```
                    ┌─ outer case, Ti 0.8 mm ──────────── r 24.88
                    │  ┌─ thermal gap 1.0 ─────────────── r 24.08
                    │  │  ┌─ S-glass banding 0.93 ─────── r 23.08
                    │  │  │  ┌─ ground insulation 0.5 ─── r 22.15
   ╔════════════════╪══╪══╪══╪═════════════════════════╗
   ║ ▓▓▓▓▓▓▓▓▓▓▓▓▓  LAYER 2 — 9 turns  ▓▓▓▓▓▓▓▓▓▓▓▓▓▓ ║  r 21.65
   ║ ─────────────  interlayer Kapton 0.15  ────────── ║  r 19.65
   ║ ▓▓▓▓▓▓▓▓▓▓▓▓▓  LAYER 1 — 13 turns ▓▓▓▓▓▓▓▓▓▓▓▓▓▓ ║  r 19.50
   ╠════════════════════════════════════════════════════╣  r 17.50
   ║              former / ground insulation 0.5        ║
   ╠════════════════════════════════════════════════════╣  r 17.00
   ║              BORE TUBE — G10, 1.5 mm               ║
   ╚════════════════════════════════════════════════════╝  r 15.50
                    running clearance 0.5
   ────────────────────────────────────────────────────    r 15.00  BORE
   ←──────────────── 40 mm coil length ────────────────→
```

| | Round AWG 10 | **Rectangular 3×2** |
|---|---|---|
| Turns per layer | 15 | **13** |
| Layers | 2 | 2 |
| Winding radial build | 5.18 mm | **4.15 mm** (incl. interlayer) |
| **Fill factor** | 55.9% | **82.5%** |
| Conductor area | 5.26 mm² | **6.00 mm²** |
| DC resistance | 1.00× | **0.88×** |
| Total wall | 10.91 mm | **9.88 mm** |
| **Outer diameter** | 51.8 mm | **49.8 mm** |

**Fill factor is the whole story:** round wire wastes 44% of the winding window
on the gaps between circles. Rectangular packs to 83%, so you get *more* copper
in *less* radial space — smaller tube and lower resistance at once.

### ⚠ Why two 1 mm strips, not one 2 mm bar

The pulse is not DC. Ring frequency `f = 1/(2π√(LC))` rises through the gun as
the capacitors get smaller, and skin depth falls with it:

| Stage | C | f | Skin depth δ | Max useful conductor (2δ) |
|---|---|---|---|---|
| 1 | 1451 µF | 1.15 kHz | 1.92 mm | 3.85 mm |
| 9 | 47 µF | 6.4 kHz | 0.82 mm | 1.63 mm |
| **31** | **13.5 µF** | **11.9 kHz** | **0.60 mm** | **1.19 mm** |

A solid 2 mm conductor is **3.3× the skin depth** at the last stages, so current
crowds to its surface and AC resistance rises ~30–50%. **Splitting into two
1 mm strips keeps every conductor inside ~1.7δ**, cutting that penalty to a few
percent, at zero cost in area or build.

This also partly addresses the limitation flagged in
[COILGUN_DESIGN.md](COILGUN_DESIGN.md) §8 — COIL-01 used DC resistance, so its
19.6% efficiency is optimistic. A thin-strip winding is the cheapest way to make
the DC assumption closer to true. *(It does not remove the caveat: efficiency
should still be budgeted at 15–20% until COIL-01 is re-run with AC resistance.)*

### Practical notes

- **Wound flatwise**, 3 mm dimension along the axis, 1 mm radially — standard
  magnet-wire practice, not edge-winding.
- **Layer 1 takes 13 turns, layer 2 takes 9** (22 total); the 4 spare positions
  in layer 2 give room for the crossover and lead-out.
- **Interlayer insulation 0.15 mm** Kapton or Nomex; each strip
  polyimide-enamelled.
- **Transpose the two strips between layers** so both see the same flux linkage,
  otherwise circulating current flows between them.
- Rectangular magnet wire is a stock item — no special manufacturing.

---

## Bottom line for a cost model

> **Don't model walls. They are 5% of launcher mass and negligible in cost.**
> **Model: capacitors + copper, × 1.3 integration factor.**

| Item | Mass | Share | Driver |
|---|---|---|---|
| **Capacitor bank** | **28.2 kg** | **83%** | stored energy |
| Coil copper | 4.17 kg | 12% | turns × stages |
| Outer case | 1.01 kg | 3.0% | minimum gauge |
| Bore tube | 0.38 kg | 1.1% | coupling gap |
| Coil banding | 0.33 kg | 1.0% | magnetic pressure |
| **Subtotal** | **34.1 kg** | | |
| **+30% switches, bus, mounts** | **44.3 kg** | | |

```
m_launcher ≈ 1.3 × ( E_stored/1500 [kg] + 4.2 [kg] + 1.7 [kg] )
```

---

## Recommended thicknesses

| Component | Material | **Wall** | Mass | Sized by |
|---|---|---|---|---|
| Bore tube | **G10/FR4** or PEEK | **1.5 mm** | 0.38 kg | coupling gap, *not* strength |
| Coil banding | **S-glass/epoxy** | **0.93 mm** (SF 2) | 0.33 kg | magnetic pressure |
| Outer case | Ti-6Al-4V or Al | **0.8 mm** | 1.01 kg | minimum gauge |

---

## Why these numbers

**A coilgun barrel is not a gun barrel — there is no chamber pressure.** Three
unrelated loads:

**1. Bore tube — carries no pressure at all.** Thickness is set by the
coil-to-armature gap, because coupling falls as the gap grows (1 mm → 1.5 mm
costs ~3.5% of k). Use the thinnest wall that survives handling: **1.5 mm**.

**2. Coil banding — the only real pressure vessel.** Magnetic pressure
B²/2μ₀ = **34.2 MPa** at the 9.3 T peak field pushes windings outward, reacted
by hoop tension `t = P·r/σ`:

| Material | SF 1.5 | SF 2.0 | Mass (31 stages, SF 2) |
|---|---|---|---|
| 4340 steel | 2.37 mm | 3.16 mm | 4.46 kg |
| Ti-6Al-4V | 1.32 mm | 1.76 mm | 1.40 kg |
| **S-glass/epoxy** | 0.70 mm | **0.93 mm** | **0.33 kg** |
| Carbon fibre T1000 | 0.40 mm | 0.53 mm | 0.15 kg |

S-glass is the pick: nearly as light as carbon fibre and **electrically
insulating**, so it can sit directly on the winding. (For context, 34 MPa is
well under the chemical gun's 85 MPa chamber.)

**3. Outer case — recoil and alignment, not pressure.** 48.7 kN peak needs only
0.94 mm even in steel; titanium and composites are **minimum-gauge limited** at
0.8 mm. Strength is not the driver.

---

## ⚠ The one structural choice that can break the launcher

**The bore tube must be electrically insulating.** A metal tube is a **shorted
turn** around the armature — it carries its own induced current, shields the
armature and destroys the coupling the launcher depends on.

| Material | Resistivity | Usable? |
|---|---|---|
| G10/FR4 | 10¹³ Ω·m | ✅ |
| PEEK | 10¹⁴ Ω·m | ✅ |
| Alumina | 10¹² Ω·m | ✅ (brittle) |
| **316 stainless** | **7.4×10⁻⁷ Ω·m** | ❌ **never, at any thickness** |

Carbon fibre is also conductive — fine as banding *outside* an insulating
layer, never as the bore.

---

## Caveats

- 9.3 T is back-computed from COIL-01's peak force treating it as pressure over
  the coil cross-section — an equivalent field, ±20% on the real peak. Banding
  mass scales with B², but at 0.33 kg even a 2× error is immaterial.
- Thin-wall hoop formula (t/r ≈ 0.04, well within validity).
- Fatigue not assessed: SF 2 is a static margin. At ~300 shots per magazine this
  is low-cycle and unlikely to bind, but a real design should check it.
- Coil banding must also react **axial** recoil per stage; not sized here (the
  outer case carries it).
