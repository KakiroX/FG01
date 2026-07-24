**las02_ablation_feedback** — The coupled ablation-and-estimation problem: each pulse changes the target's velocity, so an estimator that does not model its own effect walks off the target.

| pulses_without_feedback | elapsed_s | unmodelled_displacement_m | fraction_of_spot_radius | target_lost |
|-------------------------|-----------|---------------------------|-------------------------|-------------|
| 1                       | 0.1       | 0.01075                   | 0.04257                 | no          |
| 3                       | 0.3       | 0.09673                   | 0.3831                  | no          |
| 10                      | 1         | 1.075                     | 4.257                   | yes         |
| 30                      | 3         | 9.673                     | 38.31                   | yes         |
| 100                     | 10        | 107.5                     | 425.7                   | yes         |