"""Pick the cheapest k-mesh that gives a stable H_DOS and phi for the campaign.

The geometry is pre-relaxed ONCE with MACE and reused, so the only thing
varying between rows is the Brillouin-zone sampling.
"""
import json
import sys
import time
from pathlib import Path

from ase.io import read

from mcts_framework.superhydride import MacePrerelax, belli2025_tc
from mcts_framework.superhydride.qe import QERunner, QESettings, run_ground_state

ROOT = Path(sys.argv[1])          # working root
PSEUDO = sys.argv[2]
BIN = sys.argv[3]
RANKS = int(sys.argv[4])

# (label, kspacing_scf, kspacing_nscf) -> grids for this ~5.14 A fcc cell
CONFIGS = [
    ("gamma",   3.00, 3.00),      # 1x1x1 / 1x1x1
    ("2x2x2",   1.50, 0.60),      # 2x2x2 / 4x4x4
    ("4x4x4",   0.60, 0.28),      # 4x4x4 / 8x8x8
    ("6x6x6",   0.40, 0.19),      # 6x6x6 / 12x12x12
]

atoms = read(str(ROOT / "LaBeH8_root.cif"))
t0 = time.time()
relaxed = MacePrerelax(model="medium", device="cpu").relax(atoms, pressure_gpa=100.0)
prerelax_s = time.time() - t0
print(f"MACE pre-relax: {prerelax_s:.1f} s, V = {relaxed.get_volume():.2f} A^3", flush=True)

runner = QERunner(bin_dir=BIN, mpi_command="mpirun", ranks=RANKS,
                  environment_setup=None, timeout_s=1500)

rows = []
for label, ks_scf, ks_nscf in CONFIGS:
    settings = QESettings(ecutwfc=60.0, ecutrho=240.0, pseudo_dir=PSEUDO,
                          kspacing_scf=ks_scf, kspacing_nscf=ks_nscf)
    workdir = ROOT / "calib" / label
    t0 = time.time()
    try:
        r = run_ground_state(relaxed, settings, runner, str(workdir), relax=False)
        wall = time.time() - t0
        tc = belli2025_tc(r.phi, r.phi_star, r.h_f, r.h_dos)
        row = dict(config=label, wall_s=round(wall, 1), grid="x".join(map(str, r.grid_shape)),
                   fermi_ev=round(r.fermi_ev, 4), pressure_gpa=round(r.pressure_gpa, 1),
                   phi=round(r.phi, 4), phi_star=round(r.phi_star, 4),
                   h_f=round(r.h_f, 4), h_dos=round(r.h_dos, 4), tc_k=round(tc, 1))
    except Exception as exc:
        row = dict(config=label, wall_s=round(time.time() - t0, 1), error=str(exc)[:200])
    rows.append(row)
    print(json.dumps(row), flush=True)

(ROOT / "calib" / "results.json").write_text(json.dumps(
    {"prerelax_s": round(prerelax_s, 1), "ranks": RANKS, "rows": rows}, indent=2))
print("\ndone")
