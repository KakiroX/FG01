# FG01 / KIDR v5 — Accelerated-Reentry Simulation Suite

Executes the simulation plan in `../FG01_Agent_Context_and_Simulation_Plan_v5.md`:
momentum-transfer remediation of 1–10 cm LEO debris by a targeted kick that lowers
the orbit ~200 km, after which atmospheric drag completes the deorbit.

## Quick start

```bash
pip install -r requirements.txt
python run_all.py
```

~12 minutes on a laptop. Add `--fast` to reuse cached density tables.

## Read this first

| Document | Contents |
|---|---|
| **[docs/V6_RESULTS.md](docs/V6_RESULTS.md)** | **Start here.** Interaction model (INT-0…5) + engagement economics (ECO-1…5), and the final verdict |
| **[docs/LASER_COMPARISON.md](docs/LASER_COMPARISON.md)** | LAS-1…3: the laser benchmark rebuilt from physics — accuracy, technology readiness, corrected cost |
| [docs/MASTER_RESULTS.md](docs/MASTER_RESULTS.md) | v5 results, SIM-0 … SIM-13, with the design card and the failure envelope |
| [docs/RQ_RESULTS.md](docs/RQ_RESULTS.md) | Results ordered by the *FG01 techs.md* research questions |
| [docs/CITATION_LOG.md](docs/CITATION_LOG.md) | Every external number, with verified / inherited / estimated status |
| [docs/LIMITATIONS.md](docs/LIMITATIONS.md) | All unvalidated assumptions in one place, ranked by consequence |

## Headline

**The kick physics works. Visiting the targets is what does not.**

v5 lifted v4's primary falsification: design-point specific energy drops to
**12.1–35.7 kJ/kg**, below the 40 kJ/kg catastrophic threshold everywhere. Two
corrections to the plan's own model mattered as much: missed mass **self-disposes**
(retrograde launch above Δv_direct puts it on an immediate-reentry trajectory), and
the guidance wall is a **sweep-rate artefact** — the corrected co-orbital budget
gives a 5.8 mm miss, latency-independent from 0.1 µs to 100 ms.

v6 resolved the swarm explicitly (β = **1.20** derived from crater-ejecta scaling,
not the assumed 1.5; ε **3× better** than a uniform disc; net debris **+1** with the
largest possible ejecta fragment 7× below the trackable threshold) — and then found
the economics fail on orbital mechanics:

> **Δv per engagement = π·v_orbital / (3.5·Ω̇·T_mission) ≈ 107 m/s**

Widening the altitude band raises opportunities and hop cost *linearly and
identically*, so this **cannot be tuned**. 8,918 fragments per decade sweep through
the firing cone at Fengyun-1C; charging the Δv to reach them collapses that to
**~24 engagements**. Cost per 1 cm object: **$0.9M–4.2M** against a $5,000 target.

**The constructive result, and its boundary:** the same kick applied to larger
objects improves the economics, because the areal penalty vanishes on a
square-metre target while the Δv to visit it is unchanged. But a momentum kick is
an *intrinsically inefficient* way to move mass — effective **I_sp = 75.8 s**,
worse than cold gas — so the projectile is a fixed **14% of the target's mass**.
That is invisible at grams and decisive at tonnes. Against an electric tug of
equal platform cost, FG01 is cheaper below ~56 kg, break-even near 100 kg, and
**+$4.9M/object at 8.3 t**. The defensible window is **~100 kg to ~3.7 t**, and
only if avoiding a docking operation is worth a few million per object.

## Layout

```
fg01/            core physics
  constants.py     tagged parameter table (DERIVED/SOURCED/UNVALIDATED/DESIGN)
  orbital.py       vis-viva, Δv_shift, elements, area-to-mass
  atmosphere.py    NRLMSISE-00 tabulation (AM-6: real model, not a fitted curve)
  decay.py         orbit-averaged Gauss propagator + 3-D Cartesian validator
  delivery.py      cloud → target momentum delivery (the areal-density coupling)
  gnc.py           three side-by-side terminal miss models
  sbm.py           NASA Standard Breakup Model, both regime branches
  relmotion.py     [v6] engagement geometry, retrograde efficiency, nodal rates
  interaction.py   [v6] resolved swarm: areal field, β from ejecta scaling, ejecta
  population.py    [v6] real breakup clouds, SBM synthesis, J2 ageing
  io_utils.py      tagged table/figure output
sims/            SIM-0 … SIM-13, INT-0 … INT-5, ECO-1 … ECO-4
                 one module each, all runnable standalone
outputs/         CSV + Markdown tables, JSON conclusions and validation records
figures/         SVG + PDF + PNG, colourblind-safe
docs/            the five deliverable documents
```

Each SIM module runs on its own, e.g. `python sims/sim02_decay.py`.
SIM-9 onward read upstream results from `outputs/` at runtime rather than from
re-typed constants (AM-11).

## Validation performed

| Check | Result |
|---|---|
| Δv_shift vs. two-burn Hohmann first burn | exact (0.0 relative error) |
| E_s substitution, 2,000 random draws | 4.9×10⁻¹⁶ max relative error |
| Density table vs. direct NRLMSISE-00 calls | 1.9×10⁻⁴ max relative error |
| Orbit-averaged vs. 3-D Cartesian, 400 revolutions | 0.01% |
| Cartesian step-size convergence | 0.157% |
| Starshine 1 (258 d observed) | 233 d nominal, −9.8% (v4 reported −34%) |
| Closed-form vs. numerical delivery sizing | 0.1–0.2% |
| v4 P_hit collapse reproduced | 2.4×10⁻⁵⁸ at 10 ms (v4: ~6×10⁻⁵⁸) |
| Lifetime ∝ debris diameter | exact, as A/m ∝ 1/d requires |
| Monte-Carlo convergence, 10³→10⁶ | 0.5760 → 0.5749 |
| [v6] Closed forms vs. explicit 3-D states | \|w\| 2.1×10⁻¹⁶, cos ψ 2.8×10⁻¹¹ |
| [v6] Mass efficiency vs. grain-by-grain MC | 2.16 MC standard errors |
| [v6] β model at 6 km/s vs. literature | 1.96–8.27 against reported 2–3 |
| [v6] Δv invariant vs. tour simulation | 1.14% |
| [v6] INT-4 sizing vs. Monte-Carlo reliability | 94.7% against a 95% target |

## Deviations from the plan, and why

Three, all reported side by side with the plan's own formulation rather than
replacing it:

1. **SIM-4 miss model.** The plan's `v_rel·τ` lag term is retained as model M1 and
   reproduces v4 exactly. Models M2/M3 use the platform–target relative speed as
   the sweep rate, which follows from the engagement geometry the energy gate
   forces. All three are tabulated together.
2. **SIM-4 success variable.** `P_hit` is retained and reported; the design
   criterion is the delivered-momentum distribution, because a partial miss
   delivers partial momentum rather than nothing.
3. **SIM-12/13 cost criterion.** The plan's $100,000/kg is computed and reported;
   a per-object criterion is carried alongside because the per-kg bar is $141 for
   a 1 cm object and is unreachable by any orbital system.

A fourth addition: SIM-13 applies a **remediation-value test** — an object that
would decay unaided inside 25 years has not been remediated. Without it the
optimiser selects operating points that pass every gate while achieving nothing.

Reproducibility: fixed seed `20260723`, pinned `requirements.txt`,
NRLMSISE-00 implementation and version stated in every decay output.
