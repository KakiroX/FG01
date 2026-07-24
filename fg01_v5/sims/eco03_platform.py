"""
ECO-3 (SIM-16) -- Platform cost and consumable logistics.

Drives C_platform toward the low end and checks that the consumable chain
supports the engagement count ECO-1 actually delivers (which is far lower than
the plan assumed, so the supply question changes character entirely).
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from fg01.constants import (C_LAUNCH_RANGE, P_RE_2026, P_RE_2025, P_W_RANGE,
                            RE_PRODUCTION_TPY, W_PRODUCTION_TPY, G0)
from fg01.io_utils import save_table, save_fig, save_json, banner, CB

M_SHOT_SMALL = 2.14e-3      # kg, INT-4 at 95% reliability, 1 cm target
M_SHOT_LARGE = 1.0          # kg per shot against a large tracked object


def platform_cost(m_dry_kg, build_usd_per_kg, c_launch, n_units=1):
    """Build + launch, with a learning curve if a fleet is produced."""
    learning = n_units ** (np.log(0.85) / np.log(2))     # 85% curve
    return m_dry_kg * (build_usd_per_kg * learning + c_launch)


def run():
    banner("ECO-3  Platform cost and consumable logistics")

    rows = []
    for m_dry in (200.0, 400.0, 600.0):
        for build_rate in (50e3, 150e3, 400e3):      # USD/kg dry, smallsat->exquisite
            for c_l in C_LAUNCH_RANGE:
                for n_units in (1, 10, 100):
                    c = platform_cost(m_dry, build_rate, c_l, n_units)
                    rows.append(dict(
                        m_dry_kg=m_dry, build_usd_per_kg=build_rate,
                        c_launch_usd_per_kg=c_l, fleet_units=n_units,
                        platform_cost_musd=c / 1e6,
                        cost_per_unit_musd=c / 1e6))
    df = pd.DataFrame(rows)
    save_table(df, "eco03_platform_cost",
               "Platform cost across dry mass, build cost per kg, launch price "
               "and fleet size (85% learning curve).",
               tags={"build_usd_per_kg": "UNVALIDATED (smallsat benchmark range)"})
    cheap = df[(df.m_dry_kg == 400.0) & (df.build_usd_per_kg == 50e3)
               & (df.c_launch_usd_per_kg == 3000.0)]
    print("  minimised platform (400 kg dry, $50k/kg build, $3k/kg launch):")
    print(cheap[["fleet_units", "platform_cost_musd"]].to_string(index=False))

    # ---- consumables at the engagement counts ECO-1 actually delivers ------
    cons = []
    for n_eng, lab in ((24, "ECO-1 realised (2 km/s tour)"),
                       (55, "ECO-1 realised (10 km/s tour)"),
                       (1e4, "plan's niche target"),
                       (1e5, "plan's laser-competitive target")):
        for m_shot, tgt in ((M_SHOT_SMALL, "1 cm debris"),
                            (M_SHOT_LARGE, "large tracked object, per shot")):
            for mat, price in (("Re", P_RE_2026), ("W", P_W_RANGE[1])):
                m_tot = n_eng * m_shot
                cons.append(dict(
                    case=lab, n_engagements=n_eng, target=tgt, material=mat,
                    m_shot_kg=m_shot, consumable_kg=m_tot,
                    material_cost_usd=m_tot * price,
                    launch_cost_usd=m_tot * C_LAUNCH_RANGE[0],
                    total_consumable_usd=m_tot * (price + C_LAUNCH_RANGE[0]),
                    frac_world_production=(m_tot / 1e3)
                    / (RE_PRODUCTION_TPY[1] if mat == "Re" else W_PRODUCTION_TPY)))
    dfc = pd.DataFrame(cons)
    save_table(dfc, "eco03_consumables",
               "Consumable mass and cost at the engagement counts ECO-1 "
               "actually delivers, against the counts the plan assumed.")
    print("\n  consumables for 1 cm debris:")
    print(dfc[dfc.target == "1 cm debris"][
        ["case", "n_engagements", "material", "consumable_kg",
         "total_consumable_usd"]].to_string(index=False))

    concl = {
        "minimised_platform_musd": float(cheap[cheap.fleet_units == 1]
                                         .platform_cost_musd.iloc[0]),
        "fleet_unit_cost_musd": float(cheap[cheap.fleet_units == 100]
                                      .platform_cost_musd.iloc[0]),
        "consumable_finding": (
            "The plan's concern that rhenium supply would force a switch to "
            "tungsten at N >= 1e5 is moot: ECO-1 shows a platform achieves tens "
            "of engagements, not 1e5, so consumable mass over a whole mission "
            "is grams to kilograms and supply is irrelevant. Tungsten remains "
            "the right choice on cost and physics (v5 SIM-7), but the supply "
            "argument for it evaporates along with the engagement count."),
        "cost_structure_inversion": (
            "v5 found cost was 99.4% platform amortisation at an assumed 1e4 "
            "engagements. With the realised count of ~24, amortisation is "
            "essentially 100% and the consumable is a rounding error measured "
            "in dollars per mission. Nothing about the projectile, the "
            "material, the launcher or the grain size has any material effect "
            "on cost per object."),
    }
    save_json(concl, "eco03_conclusions")
    return df, concl


if __name__ == "__main__":
    run()
