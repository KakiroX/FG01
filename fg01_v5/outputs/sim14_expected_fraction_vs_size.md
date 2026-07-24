**sim14_expected_fraction_vs_size** — Expected fraction of the intended per-shot effect actually delivered, vs target diameter. FG01: continuous areal-capture model (delivery.expected_fraction_random_miss, already established in SIM-3/4). Laser: binary per the brief's own comparison table (hit=full effect, miss=zero), so expected fraction = P(miss < target radius) under a Rayleigh miss.

Traceability: `metric`=DERIVED (unified single-shot expected-effect fraction)

| system | scenario              | d_target_m | sigma_total_m | expected_fraction |
|--------|-----------------------|------------|---------------|-------------------|
| FG01   | optimistic            | 0.005      | 0.0005033     | 0.2796            |
| FG01   | as_designed           | 0.005      | 0.002007      | 0.2093            |
| FG01   | degraded              | 0.005      | 0.019         | 0.008404          |
| Laser  | optimistic            | 0.005      | 0.2426        | 5.309e-05         |
| Laser  | as_stated_requirement | 0.005      | 1.571         | 1.266e-06         |
| Laser  | degraded_uncorrected  | 0.005      | 6.965         | 6.443e-08         |
| FG01   | optimistic            | 0.01       | 0.0005033     | 0.7306            |
| FG01   | as_designed           | 0.01       | 0.002007      | 0.6091            |
| FG01   | degraded              | 0.01       | 0.019         | 0.0332            |
| Laser  | optimistic            | 0.01       | 0.2426        | 0.0002123         |
| Laser  | as_stated_requirement | 0.01       | 1.571         | 5.063e-06         |
| Laser  | degraded_uncorrected  | 0.01       | 6.965         | 2.577e-07         |
| FG01   | optimistic            | 0.02       | 0.0005033     | 0.9947            |
| FG01   | as_designed           | 0.02       | 0.002007      | 0.9767            |
| FG01   | degraded              | 0.02       | 0.019         | 0.1263            |
| Laser  | optimistic            | 0.02       | 0.2426        | 0.0008491         |
| Laser  | as_stated_requirement | 0.02       | 1.571         | 2.025e-05         |
| Laser  | degraded_uncorrected  | 0.02       | 6.965         | 1.031e-06         |
| FG01   | optimistic            | 0.05       | 0.0005033     | 1                 |
| FG01   | as_designed           | 0.05       | 0.002007      | 1                 |
| FG01   | degraded              | 0.05       | 0.019         | 0.57              |
| Laser  | optimistic            | 0.05       | 0.2426        | 0.005295          |
| Laser  | as_stated_requirement | 0.05       | 1.571         | 0.0001266         |
| Laser  | degraded_uncorrected  | 0.05       | 6.965         | 6.443e-06         |
| FG01   | optimistic            | 0.1        | 0.0005033     | 1                 |
| FG01   | as_designed           | 0.1        | 0.002007      | 1                 |
| FG01   | degraded              | 0.1        | 0.019         | 0.9658            |
| Laser  | optimistic            | 0.1        | 0.2426        | 0.02101           |
| Laser  | as_stated_requirement | 0.1        | 1.571         | 0.0005061         |
| Laser  | degraded_uncorrected  | 0.1        | 6.965         | 2.577e-05         |