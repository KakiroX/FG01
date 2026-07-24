"""Pre-build and verify the NRLMSISE-00 density tables (one per solar case)."""
import time
from fg01.atmosphere import build_table, verify_table
from fg01.io_utils import save_json

res = {}
for case in ("solar_min", "nominal", "solar_max"):
    t0 = time.time()
    alt, rho = build_table(case, force=True)
    dt = time.time() - t0
    mx, mn = verify_table(case)
    res[case] = {"build_seconds": round(dt, 1),
                 "max_interp_rel_err": mx, "mean_interp_rel_err": mn,
                 "rho_400km": float(rho[alt == 400][0]),
                 "rho_800km": float(rho[alt == 800][0]),
                 "rho_1000km": float(rho[alt == 1000][0])}
    print(case, res[case], flush=True)
save_json(res, "sim02_atmosphere_table_verification")
print("done")
