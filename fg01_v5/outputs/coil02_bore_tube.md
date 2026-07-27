**coil02_bore_tube** — Bore tube options. 'usable' is False for conductors: a metal tube is a shorted turn that shields the armature.

| material            | wall_mm | mass_kg | resistivity_ohm_m | coupling_factor_rel | usable | note                          |
|---------------------|---------|---------|-------------------|---------------------|--------|-------------------------------|
| G10/FR4 glass-epoxy | 1       | 0.2458  | 1e+13             | 0.6037              | yes    | standard, cheap, machinable   |
| G10/FR4 glass-epoxy | 1.5     | 0.3746  | 1e+13             | 0.5827              | yes    | standard, cheap, machinable   |
| G10/FR4 glass-epoxy | 2       | 0.5074  | 1e+13             | 0.563               | yes    | standard, cheap, machinable   |
| G10/FR4 glass-epoxy | 3       | 0.7848  | 1e+13             | 0.5267              | yes    | standard, cheap, machinable   |
| PEEK                | 1       | 0.1753  | 1e+14             | 0.6037              | yes    | tougher, lower friction       |
| PEEK                | 1.5     | 0.2673  | 1e+14             | 0.5827              | yes    | tougher, lower friction       |
| PEEK                | 2       | 0.362   | 1e+14             | 0.563               | yes    | tougher, lower friction       |
| PEEK                | 3       | 0.56    | 1e+14             | 0.5267              | yes    | tougher, lower friction       |
| Alumina Al2O3       | 1       | 0.5181  | 1e+12             | 0.6037              | yes    | hard, wear-resistant, brittle |
| Alumina Al2O3       | 1.5     | 0.7896  | 1e+12             | 0.5827              | yes    | hard, wear-resistant, brittle |
| Alumina Al2O3       | 2       | 1.07    | 1e+12             | 0.563               | yes    | hard, wear-resistant, brittle |
| Alumina Al2O3       | 3       | 1.654   | 1e+12             | 0.5267              | yes    | hard, wear-resistant, brittle |
| 316 stainless       | 1       | 1.063   | 7.4e-07           | 0.6037              | no     | DO NOT USE - shorted turn     |
| 316 stainless       | 1.5     | 1.62    | 7.4e-07           | 0.5827              | no     | DO NOT USE - shorted turn     |
| 316 stainless       | 2       | 2.194   | 7.4e-07           | 0.563               | no     | DO NOT USE - shorted turn     |
| 316 stainless       | 3       | 3.394   | 7.4e-07           | 0.5267              | no     | DO NOT USE - shorted turn     |