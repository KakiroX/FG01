**sim05_swarm_per_impact** — AM-16: per-impact specific energy for the swarm. Reported ALWAYS alongside the bulk figure, never instead of it. The bulk (single-impactor) energy at the v5 design point is already sub-catastrophic, so this branch is a robustness reserve, not a load-bearing assumption.

Traceability: `swarm equivalence`=UNVALIDATED (AM-15/AM-16)

| m_grain_mg | v_rel_m_s | Es_per_impact_kJ_kg | per_impact_catastrophic | frac_of_threshold | N_gt_1cm_per_impact_SBM |
|------------|-----------|---------------------|-------------------------|-------------------|-------------------------|
| 0.01       | 619       | 0.001355            | no                      | 3.388e-05         | 0.0001836               |
| 0.01       | 1000      | 0.003537            | no                      | 8.842e-05         | 0.000263                |
| 0.1        | 619       | 0.01355             | no                      | 0.0003388         | 0.001032                |
| 0.1        | 1000      | 0.03537             | no                      | 0.0008842         | 0.001479                |
| 1          | 619       | 0.1355              | no                      | 0.003388          | 0.005805                |
| 1          | 1000      | 0.3537              | no                      | 0.008842          | 0.008318                |
| 10         | 619       | 1.355               | no                      | 0.03388           | 0.03264                 |
| 10         | 1000      | 3.537               | no                      | 0.08842           | 0.04677                 |
| 100        | 619       | 13.55               | no                      | 0.3388            | 0.1836                  |
| 100        | 1000      | 35.37               | no                      | 0.8842            | 0.263                   |