**lnch02_mass_budget** — FG01 platform mass budget with the chemical launcher, against the coilgun baseline SIM-11 assumed. 48 shots, 30 kg magazine. The 100 g projectile is payload in both cases and cancels.

Traceability: `launcher masses`=DERIVED (GUN-02, COIL-01), `breech/autoloader mass`=UNVALIDATED (parametric), `subsystem fractions`=UNVALIDATED (SIM-11 convention)

| architecture                | round_g | launcher_kg | power_kg | thermal_kg | consumable_kg | avionics_kg | structure_kg | dry_kg | wet_kg | launch_cost_low_MUSD | launch_cost_high_MUSD |
|-----------------------------|---------|-------------|----------|------------|---------------|-------------|--------------|--------|--------|----------------------|-----------------------|
| chemical (brass case)       | 166.3   | 7.424       | 16.67    | 9.368      | 7.983         | 50          | 24.29        | 115.7  | 145.7  | 0.4372               | 1.02                  |
| chemical (steel case)       | 138.2   | 7.424       | 16.67    | 9.368      | 6.634         | 50          | 24.02        | 114.1  | 144.1  | 0.4323               | 1.009                 |
| chemical (combustible case) | 34.01   | 7.424       | 16.67    | 9.368      | 1.632         | 50          | 23.02        | 108.1  | 138.1  | 0.4143               | 0.9668                |
| coilgun                     | 0       | 98.09       | 166.7    | 30         | 0             | 50          | 74.95        | 419.7  | 449.7  | 1.349                | 3.148                 |