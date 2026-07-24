"""
FG01 / KIDR v5 -- regenerate every result from scratch.

    python run_all.py            # full pipeline
    python run_all.py --fast     # skip the density-table rebuild if cached

Pipeline order follows the v5 plan:
    SIM-0 -> SIM-1 -> SIM-2 (gate) -> SIM-3 -> SIM-4 (gate)
          -> {SIM-5 ... SIM-11} -> SIM-12 (UQ) -> SIM-13 (operating point)
"""
import sys, time, importlib, pathlib

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent / "sims"))

MODULES = [
    ("sim00_disruption", "SIM-0  catastrophic-disruption boundary"),
    ("sim01_momentum", "SIM-1  momentum transfer and mass scaling"),
    ("sim02_decay", "SIM-2  orbital decay gate (NRLMSISE-00)"),
    ("sim03_cloud", "SIM-3  cloud dispersion and cross-section export"),
    ("sim04_gnc", "SIM-4  terminal GNC and shot precision"),
    ("sim05_fragmentation", "SIM-5  fragmentation and secondary debris"),
    ("sim06_emlaunch", "SIM-6  electromagnetic launch envelope"),
    ("sim07_material", "SIM-7  material trade study"),
    ("sim08_reentry", "SIM-8  reentry ablation and environment"),
    ("sim09_cost", "SIM-9  cost model and competitors"),
    ("sim10_regulatory", "SIM-10 regulatory and aviation compliance"),
    ("sim11_platform", "SIM-11 platform, mass budget, end of life"),
    ("sim12_uq", "SIM-12 global uncertainty quantification"),
    ("sim13_operating_point", "SIM-13 operating-point selection"),
    # --- v6: interaction model -------------------------------------------
    ("int00_geometry", "INT-0  engagement geometry / co-orbital premise"),
    ("int01_cloud", "INT-1  cloud areal-density field"),
    ("int02_beta", "INT-2  momentum coupling beta at sub-km/s"),
    ("int03_delivery", "INT-3/4 aggregate delivery and aim tolerance"),
    ("int05_ejecta", "INT-5  cloud ejecta and net debris"),
    ("int06_target_size", "INT-6  target size: interception / momentum / drag"),
    # --- v6: economic architecture ---------------------------------------
    ("eco02_clusters", "ECO-2  cluster availability"),
    ("eco01_cadence", "ECO-1  engagement cadence and Delta-v"),
    ("eco03_platform", "ECO-3  platform cost and consumables"),
    ("eco04_optimize", "ECO-4  economic optimisation and verdict"),
    ("eco05_large_objects", "ECO-5  impulse efficiency vs. target mass"),
]


def main():
    t0 = time.time()
    if "--fast" not in sys.argv:
        from fg01.atmosphere import build_table
        for case in ("solar_min", "nominal", "solar_max"):
            print(f"building NRLMSISE-00 table: {case}", flush=True)
            build_table(case, force=False)
    for mod, label in MODULES:
        t = time.time()
        print(f"\n>>> {label}", flush=True)
        importlib.import_module(mod).run()
        print(f"    ({time.time()-t:.1f} s)", flush=True)
    from fg01.io_utils import write_manifest
    write_manifest()
    print(f"\nAll simulations complete in {time.time()-t0:.0f} s.")
    print("Tables -> outputs/   Figures -> figures/   Docs -> docs/")


if __name__ == "__main__":
    main()
