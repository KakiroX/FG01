"""
FG01 / KIDR v5 -- Parameter table (Part 3 of the v5 plan) with traceability tags.

Every value carries a tag:
    DERIVED     -- follows from stated physics / definitions
    SOURCED     -- specific external citation (see docs/CITATION_LOG.md, key given)
    UNVALIDATED -- assumption with no direct source at the relevant scale
    DESIGN      -- design variable swept by the study

SI units internally (m, kg, s) unless a name says otherwise.
"""
from dataclasses import dataclass

# --------------------------------------------------------------------------
# Traceability helper
# --------------------------------------------------------------------------
DERIVED, SOURCED, UNVALIDATED, DESIGN = "DERIVED", "SOURCED", "UNVALIDATED", "DESIGN"


@dataclass(frozen=True)
class P:
    """A tagged parameter."""
    value: float
    units: str
    tag: str
    cite: str = ""          # key into docs/CITATION_LOG.md
    note: str = ""

    def __float__(self):
        return float(self.value)


# --------------------------------------------------------------------------
# Gravitational / geodetic
# --------------------------------------------------------------------------
MU = 3.986004418e14                 # m^3/s^2   DERIVED (WGS-84/EGM96 GM)
R_E = 6.371e6                       # m         DERIVED (mean radius)
OMEGA_E = 7.292115e-5               # rad/s     DERIVED (sidereal rotation)
J2 = 1.08262668e-3                  # -         SOURCED [egm96]
G0 = 9.80665                        # m/s^2     DERIVED

MU_KM = 398600.4418                 # km^3/s^2  (for cross-checks vs the plan)
R_E_KM = 6371.0                     # km

# --------------------------------------------------------------------------
# Materials
# --------------------------------------------------------------------------
RHO_AL = 2700.0                     # kg/m^3    SOURCED [matweb_al6061]  (Al-6061, AM-2)
RHO_RE = 21010.0                    # kg/m^3    SOURCED [usgs_re]
RHO_W = 19250.0                     # kg/m^3    SOURCED [usgs_w]
RHO_TA = 16650.0                    # kg/m^3    SOURCED [crc]
RHO_STEEL = 7850.0                  # kg/m^3    SOURCED [crc]
RHO_CU = 8960.0                     # kg/m^3    SOURCED [crc]

# Effective heat of ablation (MJ/kg -> J/kg), range in *_RANGE
H_ABL_RE = 4.5e6                    # J/kg      SOURCED [ablation]  range 3.5-5.5
H_ABL_RE_RANGE = (3.5e6, 5.5e6)
H_ABL_AL = 12.0e6                   # J/kg      SOURCED [ablation]  range 10-14
H_ABL_AL_RANGE = (10.0e6, 14.0e6)
H_ABL_W = 8.0e6                     # J/kg      SOURCED [ablation]  range 6.0-10.0
H_ABL_W_RANGE = (6.0e6, 10.0e6)

# --------------------------------------------------------------------------
# Impact / breakup physics
# --------------------------------------------------------------------------
ES_C = 40.0e3                       # J/kg      SOURCED [johnson2001]  catastrophic threshold
ES_C_RANGE = (30.0e3, 50.0e3)       # J/kg      SOURCED (spread in the literature)
BETA_RANGE = (1.0, 2.5)             # -         UNVALIDATED (AM-4a) flat
BETA_CONTEXT_MAX = 12.0             # -         UNVALIDATED, context only (strengthless bound)
GAMMA = 2.0                         # -         DERIVED (engineering mass margin)

# --------------------------------------------------------------------------
# Drag
# --------------------------------------------------------------------------
CD_SPHERE = 2.2                     # -         SOURCED [kinghele] free-molecular sphere
CD_RANGE = (2.0, 2.4)

# --------------------------------------------------------------------------
# Design variables
# --------------------------------------------------------------------------
DH_NOMINAL = 200.0e3                # m         DESIGN
DH_RANGE = (150.0e3, 250.0e3)       # m         DESIGN
ALT_SWEEP_KM = [400, 500, 600, 700, 800, 900, 1000, 1100, 1200]     # DESIGN
DEBRIS_DIAM_M = [0.005, 0.01, 0.02, 0.05, 0.10]                     # DESIGN
GRAIN_DIAM_M = [0.05e-3, 0.1e-3, 0.5e-3, 1.0e-3, 2.0e-3]            # DESIGN
TAU_MS = [0.1, 1.0, 5.0, 10.0, 50.0, 100.0]                         # DESIGN (control lag)
RANGE_M = [1.0, 5.0, 10.0, 50.0, 100.0]                             # DESIGN (engagement range)

# --------------------------------------------------------------------------
# Space environment drivers
# --------------------------------------------------------------------------
F107_CASES = {"solar_min": 70.0, "nominal": 150.0, "solar_max": 250.0}   # sfu SOURCED [f107]
AP_CASES = {"solar_min": 4.0, "nominal": 15.0, "solar_max": 45.0}        # -   UNVALIDATED pairing
NRLMSISE_IMPL = "nrlmsise00 (Python wrapper of Dominik Brodowski C port of NRLMSISE-00), v0.1.2"

# --------------------------------------------------------------------------
# Economics (all dated -- see CITATION_LOG)
# --------------------------------------------------------------------------
P_RE_2026 = 7283.0                  # USD/kg    SOURCED [re_price_2026] spot 2026-07-22
P_RE_2026Q1 = 6389.0                # USD/kg    SOURCED [re_price_2026]
P_RE_2025 = 2486.0                  # USD/kg    SOURCED [re_price_2025]
P_RE_RANGE = (2486.0, 7283.0)
RE_PRODUCTION_TPY = (60.0, 81.0)    # t/yr      SOURCED [usgs_re]
P_W_RANGE = (30.0, 100.0)           # USD/kg    SOURCED [usgs_w]
W_PRODUCTION_TPY = 84000.0          # t/yr      SOURCED [usgs_w]
P_TA = 200.0                        # USD/kg    SOURCED [usgs_ta]
P_STEEL = 1.0                       # USD/kg    SOURCED [usgs_fe]
C_LAUNCH_RANGE = (3000.0, 7000.0)   # USD/kg    SOURCED [launch_cost]
C_LAUNCH_STARSHIP = (100.0, 500.0)  # USD/kg    labelled sensitivity ONLY (AM-9)
P_ELEC = 0.10                       # USD/kWh   SOURCED [elec]
C_PLATFORM_RANGE = (50e6, 500e6)    # USD       UNVALIDATED (parametric)
PLATFORM_LIFE_YR = 10.0

# Competitor cost baselines, USD/kg removed  SOURCED [competitors]
COMPETITORS = {
    "Ground-based laser ablation": (100.0, 500.0, "theoretical"),
    "Electrodynamic tether":       (1000.0, 5000.0, "demonstrated (sub-scale)"),
    "Aerogel / foam capture":      (5000.0, 20000.0, "concept"),
    "ADR tug (robotic capture)":   (10000.0, 100000.0, "demonstrated (large targets)"),
}

# --------------------------------------------------------------------------
# Regulatory
# --------------------------------------------------------------------------
CASUALTY_THRESHOLD = 1e-4           # per event SOURCED [nasa_std_8719]
FCC_DISPOSAL_YR = 5.0               # yr        SOURCED [fcc_25_283]
ODMSP_YR = 25.0                     # yr        SOURCED [odmsp]

# --------------------------------------------------------------------------
# EM launcher benchmarks (demonstrated only -- AM-9)
# --------------------------------------------------------------------------
EML_BENCHMARKS = [
    # (label, projectile mass kg, muzzle velocity m/s, muzzle energy J, citation key)
    ("US Navy/ONR EMRG 2008 record shot", 3.2, 2520.0, 10.16e6, "onr2008"),
    ("US Navy/ONR EMRG 2010 33 MJ shot", 10.4, 2520.0, 33.0e6, "onr2010"),
    ("Coilgun, laboratory (Turman 1996)", 0.010, 2000.0, 2.0e4, "turman1996"),
    ("Light-gas/EM hybrid, 3 g class", 0.003, 5900.0, 5.2e4, "mcnab2003"),
]
COILGUN_ETA_RANGE = (0.30, 0.50)    # -         SOURCED [turman1996, mcnab2003]

# --------------------------------------------------------------------------
# Sensors
# --------------------------------------------------------------------------
OPTICAL_APERTURE_M = 0.10           # m         DESIGN
OPTICAL_LAMBDA_M = 550e-9           # m         DESIGN
LIDAR_BW_HZ = 1.0e6                 # Hz        DESIGN
LIDAR_SNR = 10.0                    # -         DESIGN

# --------------------------------------------------------------------------
# Output paths
# --------------------------------------------------------------------------
import pathlib
ROOT = pathlib.Path(__file__).resolve().parent.parent
OUT = ROOT / "outputs"
FIG = ROOT / "figures"
OUT.mkdir(exist_ok=True)
FIG.mkdir(exist_ok=True)

SEED = 20260723                     # fixed Monte-Carlo seed (reproducibility)
