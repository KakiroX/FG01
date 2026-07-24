**eco01_dv_invariant** — Closed form: Delta-v per engagement = pi*v/(3.5*Omega_dot*T). Independent of the altitude band, so it cannot be tuned away. Validated against the structural table (107 m/s at Fengyun-1C, 10 yr) and against the tour simulation.

| case                         | h_km | inc_deg | mission_yr | nodal_rate_deg_day | dv_per_engagement_m_s | engagements_per_2000_m_s | dv_for_10000_engagements_km_s |
|------------------------------|------|---------|------------|--------------------|-----------------------|--------------------------|-------------------------------|
| Fengyun-1C (sun-sync)        | 865  | 98.8    | 10         | 0.9779             | 106.9                 | 18.72                    | 1069                          |
| Fengyun-1C (sun-sync)        | 865  | 98.8    | 20         | 0.9779             | 53.43                 | 37.43                    | 534.3                         |
| Cosmos-2251                  | 789  | 74      | 10         | 1.828              | 57.46                 | 34.8                     | 574.6                         |
| Cosmos-2251                  | 789  | 74      | 20         | 1.828              | 28.73                 | 69.61                    | 287.3                         |
| Iridium-33                   | 789  | 86.4    | 10         | 0.4165             | 252.3                 | 7.929                    | 2523                          |
| Iridium-33                   | 789  | 86.4    | 20         | 0.4165             | 126.1                 | 15.86                    | 1261                          |
| Cosmos-1408                  | 480  | 82.5    | 10         | 1.01               | 106.3                 | 18.81                    | 1063                          |
| Cosmos-1408                  | 480  | 82.5    | 20         | 1.01               | 53.15                 | 37.63                    | 531.5                         |
| hypothetical low-inclination | 900  | 51.6    | 10         | 3.904              | 26.7                  | 74.89                    | 267                           |
| hypothetical low-inclination | 900  | 51.6    | 20         | 3.904              | 13.35                 | 149.8                    | 133.5                         |
| hypothetical i=30 deg        | 900  | 30      | 10         | 5.443              | 19.15                 | 104.4                    | 191.5                         |
| hypothetical i=30 deg        | 900  | 30      | 20         | 5.443              | 9.577                 | 208.8                    | 95.77                         |