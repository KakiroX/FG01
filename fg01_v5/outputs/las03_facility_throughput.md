**las03_facility_throughput** — Ground-laser facility cost and throughput. Throughput is limited by the terminator window (the target must be sunlit while the site is dark for the tracking camera to work), by weather, and by acquisition overhead.

Traceability: `facility cost`=UNVALIDATED (parametric), `twilight/weather`=SOURCED (site statistics)

| aperture_m | facility_cost_musd | telescope_musd | laser_musd | engagement_time_s | twilight_hr | weather_frac | objects_per_night | objects_per_decade | cost_per_object_usd | electricity_per_object_usd |
|------------|--------------------|----------------|------------|-------------------|-------------|--------------|-------------------|--------------------|---------------------|----------------------------|
| 5          | 156.8              | 150            | 6.773      | 29.42             | 3           | 0.6          | 72.46             | 2.647e+05          | 592.3               | 0.5536                     |
| 5          | 156.8              | 150            | 6.773      | 29.42             | 4           | 0.7          | 112.7             | 4.117e+05          | 380.8               | 0.5536                     |
| 10         | 604.8              | 600            | 4.75       | 29.42             | 3           | 0.6          | 72.46             | 2.647e+05          | 2285                | 0.3883                     |
| 10         | 604.8              | 600            | 4.75       | 29.42             | 4           | 0.7          | 112.7             | 4.117e+05          | 1469                | 0.3883                     |
| 15         | 1509               | 1500           | 9.036      | 29.42             | 3           | 0.6          | 72.46             | 2.647e+05          | 5702                | 0.7385                     |
| 15         | 1509               | 1500           | 9.036      | 29.42             | 4           | 0.7          | 112.7             | 4.117e+05          | 3665                | 0.7385                     |