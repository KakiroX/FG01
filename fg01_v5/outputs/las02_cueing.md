**las02_cueing** — External orbit knowledge each system needs before an engagement can start. FG01 uses external cueing only to get within range of its own sensor, which then closes the loop; the laser must be handed a solution accurate to its spot radius at 1,186 km because it has no way to improve on it.

| system                                   | required_position_knowledge_m | cue_used_at_range_m | angular_equivalent_urad | loop_closed_onboard |
|------------------------------------------|-------------------------------|---------------------|-------------------------|---------------------|
| FG01: close to optical acquisition range | 1000                          | 5e+04               | 2e+04                   | yes                 |
| Ground laser: place spot on target       | 0.2525                        | 1.186e+06           | 0.2129                  | no                  |