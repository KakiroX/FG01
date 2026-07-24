**sim14_laser_pulse_compounding** — P(every one of N pulses hits) = p_hit^N for a laser engagement needing sustained tracking across N pulses to accumulate useful Delta-v (per the brief: the system 'must stay focused... on a small, fast-moving target'). FG01's design point needs exactly 1 kick per object (SIM-13).

Traceability: `p_hit^N`=DERIVED (independent-pulse approximation, conservative if tracking errors are correlated shot-to-shot)

| scenario              | d_target_m | n_pulses | p_hit_single | p_all_n_hit |
|-----------------------|------------|----------|--------------|-------------|
| optimistic            | 0.005      | 1        | 5.309e-05    | 5.309e-05   |
| optimistic            | 0.005      | 2        | 5.309e-05    | 2.819e-09   |
| optimistic            | 0.005      | 3        | 5.309e-05    | 1.496e-13   |
| optimistic            | 0.005      | 5        | 5.309e-05    | 4.217e-22   |
| optimistic            | 0.005      | 10       | 5.309e-05    | 1.779e-43   |
| optimistic            | 0.01       | 1        | 0.0002123    | 0.0002123   |
| optimistic            | 0.01       | 2        | 0.0002123    | 4.509e-08   |
| optimistic            | 0.01       | 3        | 0.0002123    | 9.574e-12   |
| optimistic            | 0.01       | 5        | 0.0002123    | 4.317e-19   |
| optimistic            | 0.01       | 10       | 0.0002123    | 1.864e-37   |
| optimistic            | 0.02       | 1        | 0.0008491    | 0.0008491   |
| optimistic            | 0.02       | 2        | 0.0008491    | 7.21e-07    |
| optimistic            | 0.02       | 3        | 0.0008491    | 6.122e-10   |
| optimistic            | 0.02       | 5        | 0.0008491    | 4.414e-16   |
| optimistic            | 0.02       | 10       | 0.0008491    | 1.948e-31   |
| optimistic            | 0.05       | 1        | 0.005295     | 0.005295    |
| optimistic            | 0.05       | 2        | 0.005295     | 2.804e-05   |
| optimistic            | 0.05       | 3        | 0.005295     | 1.485e-07   |
| optimistic            | 0.05       | 5        | 0.005295     | 4.162e-12   |
| optimistic            | 0.05       | 10       | 0.005295     | 1.733e-23   |
| optimistic            | 0.1        | 1        | 0.02101      | 0.02101     |
| optimistic            | 0.1        | 2        | 0.02101      | 0.0004415   |
| optimistic            | 0.1        | 3        | 0.02101      | 9.278e-06   |
| optimistic            | 0.1        | 5        | 0.02101      | 4.096e-09   |
| optimistic            | 0.1        | 10       | 0.02101      | 1.678e-17   |
| as_stated_requirement | 0.005      | 1        | 1.266e-06    | 1.266e-06   |
| as_stated_requirement | 0.005      | 2        | 1.266e-06    | 1.602e-12   |
| as_stated_requirement | 0.005      | 3        | 1.266e-06    | 2.028e-18   |
| as_stated_requirement | 0.005      | 5        | 1.266e-06    | 3.248e-30   |
| as_stated_requirement | 0.005      | 10       | 1.266e-06    | 1.055e-59   |
| as_stated_requirement | 0.01       | 1        | 5.063e-06    | 5.063e-06   |
| as_stated_requirement | 0.01       | 2        | 5.063e-06    | 2.563e-11   |
| as_stated_requirement | 0.01       | 3        | 5.063e-06    | 1.298e-16   |
| as_stated_requirement | 0.01       | 5        | 5.063e-06    | 3.326e-27   |
| as_stated_requirement | 0.01       | 10       | 5.063e-06    | 1.106e-53   |
| as_stated_requirement | 0.02       | 1        | 2.025e-05    | 2.025e-05   |
| as_stated_requirement | 0.02       | 2        | 2.025e-05    | 4.101e-10   |
| as_stated_requirement | 0.02       | 3        | 2.025e-05    | 8.305e-15   |
| as_stated_requirement | 0.02       | 5        | 2.025e-05    | 3.406e-24   |
| as_stated_requirement | 0.02       | 10       | 2.025e-05    | 1.16e-47    |

_(75 rows total; full data in sim14_laser_pulse_compounding.csv)_
