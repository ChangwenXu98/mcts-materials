"""Do the failing pre-relaxations plateau just above fmax, or are they stuck?

That distinction decides the fix. If they settle at ~0.06-0.10 eV/A, the
criterion is simply too strict for a *cell* relaxation at 150 GPa - and since
the DFT pressure match re-measures and corrects the cell anyway, the pre-relax
only has to deliver sensible internal coordinates. If they are sitting at 0.5,
no step budget will save them and they should fail fast instead.
"""

import sys
import time

import numpy as np
from ase.filters import FrechetCellFilter
from ase.io import read
from ase.optimize import FIRE
from mace.calculators import mace_mp

GPA_TO_EV_PER_ANG3 = 1.0 / 160.21766208
PRESSURE_GPA = 150.0

CASES = [
    ("BaMgH8", {"La": "Ba", "Be": "Mg"}),   # failed in the campaign
    ("MgYH8", {"La": "Y", "Be": "Mg"}),     # failed
    ("BaLiH8", {"La": "Ba", "Be": "Li"}),   # failed
    ("BeCeH8", {"La": "Ce"}),               # failed
    ("LaBeH8", {}),                         # control: converges in 36 steps
]


def main() -> int:
    calc = mace_mp(model="medium", device="cpu", default_dtype="float64", dispersion=False)
    base = read("campaign/LaBeH8_root.cif")

    with open(sys.argv[1], "w", buffering=1) as out:
        for label, swaps in CASES:
            atoms = base.copy()
            atoms.set_chemical_symbols(
                [swaps.get(s, s) for s in atoms.get_chemical_symbols()]
            )
            atoms.calc = calc
            target = FrechetCellFilter(
                atoms, scalar_pressure=PRESSURE_GPA * GPA_TO_EV_PER_ANG3
            )
            optimizer = FIRE(target, logfile=None)

            started = time.time()
            trace = []
            for _ in range(6):  # sample fmax every 100 steps, up to 600
                optimizer.run(fmax=0.05, steps=optimizer.get_number_of_steps() + 100)
                fmax = float(np.sqrt((target.get_forces() ** 2).sum(axis=1).max()))
                trace.append(f"{optimizer.get_number_of_steps()}:{fmax:.4f}")
                if fmax < 0.05:
                    break
            out.write(
                f"{label:8s} {' '.join(trace)}  "
                f"wall {time.time() - started:5.1f}s  V={atoms.get_volume():.2f}\n"
            )
    return 0


if __name__ == "__main__":
    sys.exit(main())
