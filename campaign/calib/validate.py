"""Validate the campaign protocol: MACE pre-relax -> pressure-matched SCF."""
import json
import sys
import time
from pathlib import Path

from ase.io import read

from mcts_framework.superhydride import MacePrerelax, belli2025_tc
from mcts_framework.superhydride.qe import PressureMatch, QERunner, QESettings, run_ground_state

ROOT, PSEUDO, BIN, RANKS = Path(sys.argv[1]), sys.argv[2], sys.argv[3], int(sys.argv[4])
TARGET = 150.0

base = read(str(ROOT / "LaBeH8_root.cif"))
pre = MacePrerelax(model="medium", device="cpu")
runner = QERunner(bin_dir=BIN, mpi_command="mpirun", ranks=RANKS, timeout_s=2400)
settings = QESettings(ecutwfc=90.0, ecutrho=360.0, pseudo_dir=PSEUDO,
                      kspacing_scf=0.60, kspacing_nscf=0.28)
match = PressureMatch(target_gpa=TARGET, tolerance_gpa=10.0, max_scf=3)

for label, swap in [("LaBeH8", None), ("YBeH8", ("La", "Y")), ("LaBH8", ("Be", "B"))]:
    atoms = base.copy()
    if swap:
        atoms.set_chemical_symbols([swap[1] if s == swap[0] else s
                                    for s in atoms.get_chemical_symbols()])
    t0 = time.time()
    try:
        relaxed = pre.relax(atoms, pressure_gpa=TARGET)
        r = run_ground_state(relaxed, settings, runner, str(ROOT / "val" / label),
                             pressure_gpa=TARGET, pressure_match=match, relax=False)
        row = dict(compound=label, wall_s=round(time.time() - t0, 1),
                   dft_pressure_gpa=round(r.pressure_gpa, 1),
                   grid="x".join(map(str, r.grid_shape)), fermi_ev=round(r.fermi_ev, 3),
                   phi=round(r.phi, 4), phi_star=round(r.phi_star, 4),
                   h_f=round(r.h_f, 3), h_dos=round(r.h_dos, 4),
                   tc_k=round(belli2025_tc(r.phi, r.phi_star, r.h_f, r.h_dos), 1))
    except Exception as exc:
        row = dict(compound=label, wall_s=round(time.time() - t0, 1), error=str(exc)[:200])
    print(json.dumps(row), flush=True)
print("done")
