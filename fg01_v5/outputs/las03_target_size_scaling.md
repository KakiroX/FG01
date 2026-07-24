**las03_target_size_scaling** — How the comparison scales with target size. FG01's mass penalty falls as 1/A_target and reaches 1; the laser's engagement time grows because a larger target needs proportionally more total impulse while the fluence per unit area is unchanged.

| target_diameter_cm | target_mass_g | laser_pulses_per_object | laser_engagement_time_s | laser_beam_fraction_on_target | fg01_areal_penalty | fg01_launched_mass_g |
|--------------------|---------------|-------------------------|-------------------------|-------------------------------|--------------------|----------------------|
| 1                  | 1.414         | 294.2                   | 29.42                   | 0.0005135                     | 3.626              | 0.3573               |
| 3                  | 38.17         | 882.7                   | 88.27                   | 0.004622                      | 1                  | 2.661                |
| 10                 | 1414          | 2942                    | 294.2                   | 0.05135                       | 1                  | 98.54                |
| 30                 | 3.817e+04     | 8827                    | 882.7                   | 0.4622                        | 1                  | 2661                 |
| 100                | 1.414e+06     | 1.511e+05               | 1.511e+04               | 1                             | 1                  | 9.854e+04            |