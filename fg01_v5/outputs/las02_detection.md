**las02_detection** — Detection geometry and photon budget for a 1 cm, albedo-0.1 sphere. 'resolved_ratio' is target angular size divided by the diffraction limit: above 1 the target is an extended object whose centroid is directly measurable; below 1 it is an unresolved point source.

Traceability: `photometry`=DERIVED, `albedo 0.1`=UNVALIDATED

| sensor                 | range_m   | aperture_m | target_angular_size_urad | diffraction_limit_urad | resolved_ratio | resolved | apparent_magnitude | photoelectrons_per_s | photoelectrons_per_10ms | snr_10ms  |
|------------------------|-----------|------------|--------------------------|------------------------|----------------|----------|--------------------|----------------------|-------------------------|-----------|
| FG01 terminal sensor   | 10        | 0.1        | 1000                     | 6.71                   | 149            | yes      | -7.735             | 2.368e+11            | 2.368e+09               | 4.866e+04 |
| Ground telescope       | 1.186e+06 | 5          | 0.008432                 | 0.1342                 | 0.06283        | no       | 17.64              | 4.208e+04            | 420.8                   | 18.44     |
| Ground telescope, 10 m | 1.186e+06 | 10         | 0.008432                 | 0.0671                 | 0.1257         | no       | 17.64              | 1.683e+05            | 1683                    | 39.86     |