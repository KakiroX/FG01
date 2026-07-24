**las04_recoil_geometry** — Useful retrograde fraction of ablation recoil as a function of where the illuminator sits relative to the target. Recoil is along the line of sight, so only an illuminator AHEAD of the target along its velocity vector produces a retrograde push.

| illuminator_position                      | useful_retrograde_fraction | energy_penalty_for_same_kick | direction                 | usable |
|-------------------------------------------|----------------------------|------------------------------|---------------------------|--------|
| directly ahead (co-orbital, leading)      | 1                          | 1                            | retrograde (lowers orbit) | yes    |
| ahead and 10 deg above                    | 0.9848                     | 1.015                        | retrograde (lowers orbit) | yes    |
| ahead and 45 deg above                    | 0.7071                     | 1.414                        | retrograde (lowers orbit) | yes    |
| directly above (higher orbit, same plane) | 0                          | ∞                            | no along-track component  | no     |
| directly below                            | 0                          | ∞                            | no along-track component  | no     |
| abeam (cross-track, different plane)      | 0                          | ∞                            | no along-track component  | no     |
| directly behind (trailing)                | -1                         | ∞                            | prograde (RAISES orbit)   | no     |