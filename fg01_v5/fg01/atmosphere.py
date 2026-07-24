"""
Atmosphere: a *tabulation of the real NRLMSISE-00 model* (AM-6).

No fitted density curve is used anywhere in this study.  What is built here is
a lookup table whose entries are direct evaluations of the NRLMSISE-00 C model
(via the `nrlmsise00` Python wrapper), spatially and diurnally averaged, then
log-linearly interpolated in altitude on a 1 km grid.  The interpolation error
against direct model calls is measured and reported by `verify_table()`.

Why average: an object decaying over months-to-decades samples all longitudes,
all local solar times and (for inclined orbits) a wide latitude band, so the
orbit-averaged drag is driven by the globally averaged density, not by a single
lat/lon/LST sample.  Latitude weighting is cos(lat) (equal-area).

Solar/geomagnetic drivers are held fixed per case (F10.7 = 70/150/250 with
Ap = 4/15/45).  This is the plan's specification; it brackets, but does not
reproduce, a real 11-year cycle.  Reported as a three-case envelope, never as a
single number.
"""
import datetime as dt
import numpy as np

from .constants import OUT, NRLMSISE_IMPL, F107_CASES, AP_CASES

try:
    from nrlmsise00 import msise_flat
except ImportError as exc:  # pragma: no cover
    raise ImportError("nrlmsise00 is required (pip install nrlmsise00)") from exc

# 1 km grid, 0 km -> 1600 km (extends below 80 km for the SIM-8 reentry
# integration; the decay propagator never uses altitudes below 100 km)
ALT_GRID_KM = np.arange(0.0, 1601.0, 1.0)

_LATS = np.arange(-82.5, 82.6, 15.0)          # equal-area weighted
_LONS = np.arange(0.0, 360.0, 45.0)
_HOURS = np.arange(0.0, 24.0, 6.0)            # local-time coverage
_DAYS = [dt.datetime(2020, 3, 20), dt.datetime(2020, 6, 21),
         dt.datetime(2020, 9, 22), dt.datetime(2020, 12, 21)]   # seasonal coverage

_CACHE = {}


def _table_path(case):
    return OUT / f"atmos_table_{case}.npz"


def build_table(case="nominal", force=False):
    """Globally/diurnally/seasonally averaged NRLMSISE-00 density [kg/m^3]."""
    if case in _CACHE and not force:
        return _CACHE[case]
    path = _table_path(case)
    if path.exists() and not force:
        d = np.load(path)
        _CACHE[case] = (d["alt_km"], d["rho"])
        return _CACHE[case]

    f107 = F107_CASES[case]
    ap = AP_CASES[case]
    w_lat = np.cos(np.deg2rad(_LATS))
    w_lat = w_lat / w_lat.sum()

    rho = np.zeros_like(ALT_GRID_KM)
    for i, alt in enumerate(ALT_GRID_KM):
        acc = 0.0
        for day in _DAYS:
            for hour in _HOURS:
                t = day + dt.timedelta(hours=float(hour))
                for j, lat in enumerate(_LATS):
                    out = msise_flat(t, float(alt), float(lat), _LONS,
                                     f107, f107, ap)
                    # out[..., 5] = total mass density in g/cm^3 -> kg/m^3
                    acc += w_lat[j] * np.mean(np.asarray(out)[..., 5]) * 1000.0
        rho[i] = acc / (len(_DAYS) * len(_HOURS))

    np.savez_compressed(path, alt_km=ALT_GRID_KM, rho=rho,
                        impl=NRLMSISE_IMPL, f107=f107, ap=ap)
    _CACHE[case] = (ALT_GRID_KM, rho)
    return _CACHE[case]


def density(alt_m, case="nominal"):
    """Interpolated averaged density [kg/m^3] at geometric altitude alt_m [m]."""
    alt_km_grid, rho = build_table(case)
    alt_km = np.asarray(alt_m, dtype=float) / 1000.0
    alt_km = np.clip(alt_km, alt_km_grid[0], alt_km_grid[-1])
    lr = np.interp(alt_km, alt_km_grid, np.log(rho))
    return np.exp(lr)


def verify_table(case="nominal"):
    """
    Interpolation-error check: compare the table against fresh direct
    NRLMSISE-00 calls at off-grid altitudes, using the same averaging.
    Returns (max_rel_err, mean_rel_err).
    """
    f107, ap = F107_CASES[case], AP_CASES[case]
    w_lat = np.cos(np.deg2rad(_LATS)); w_lat /= w_lat.sum()
    test_alts = np.array([137.3, 251.7, 403.4, 555.5, 662.2, 800.5,
                          911.9, 1000.4, 1150.8, 1300.3])
    errs = []
    for alt in test_alts:
        acc = 0.0
        for day in _DAYS:
            for hour in _HOURS:
                t = day + dt.timedelta(hours=float(hour))
                for j, lat in enumerate(_LATS):
                    out = msise_flat(t, float(alt), float(lat), _LONS,
                                     f107, f107, ap)
                    acc += w_lat[j] * np.mean(np.asarray(out)[..., 5]) * 1000.0
        direct = acc / (len(_DAYS) * len(_HOURS))
        interp = float(density(alt * 1000.0, case))
        errs.append(abs(interp - direct) / direct)
    errs = np.array(errs)
    return float(errs.max()), float(errs.mean())


def density_pointwise(alt_m, lat_deg, lon_deg, when, f107, ap):
    """Direct (non-averaged) NRLMSISE-00 call -- used by the 3-D validator."""
    out = msise_flat(when, float(np.asarray(alt_m) / 1000.0),
                     float(lat_deg), float(lon_deg), f107, f107, ap)
    return float(np.asarray(out)[..., 5]) * 1000.0


def scale_height_km(alt_m, case="nominal"):
    """Local density scale height, from the table (diagnostic)."""
    a = np.asarray(alt_m, dtype=float)
    r1 = density(a - 2500.0, case)
    r2 = density(a + 2500.0, case)
    return 5.0 / np.log(r1 / r2)
