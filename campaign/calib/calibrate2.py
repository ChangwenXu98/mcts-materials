"""Round 2: how far the MACE cell sits from its nominal pressure, and the cost
of the paper's 90 Ry cutoff.

MACE-MP is fitted on ambient-pressure data, so its 100 GPa cell is not the DFT
100 GPa cell. The DFT pressure printed by the SCF is the measurement; this
scans the MACE target to find which one lands the campaign in a sensible
window, and checks a substituted composition end to end.
"""
import json, sys, time
from pathlib import Path
from ase.io import read

from mcts_framework.superhydride import MacePrerelax, belli2025_tc
from mcts_framework.superhydride.qe import QERunner, QESettings, run_ground_state

ROOT, PSEUDO, BIN, RANKS = Path(sys.argv[1]), sys.argv[2], sys.argv[3], int(sys.argv[4])

# 4x4x4 SCF / 8x8x8 NSCF for this ~5.1 A fcc cell (round 1: converged).
SETTINGS = dict(ecutwfc=90.0, ecutrho=360.0, pseudo_dir=PSEUDO,
                kspacing_scf=0.60, kspacing_nscf=0.28)

base = read(str(ROOT / "LaBeH8_root.cif"))
pre = MacePrerelax(model="medium", device="cpu")
runner = QERunner(bin_dir=BIN, mpi_command="mpirun", ranks=RANKS, timeout_s=1500)

CASES = [("LaBeH8", "Be", 40.0), ("LaBeH8", "Be", 70.0), ("LaBeH8", "Be", 100.0),
         ("YBeH8", "Be", 70.0), ("LaBH8", "B", 70.0)]

rows = []
for label, _site, p_target in CASES:
    atoms = base.copy()
    if label == "YBeH8":
        atoms.set_chemical_symbols(["Y" if s == "La" else s for s in atoms.get_chemical_symbols()])
    elif label == "LaBH8":
        atoms.set_chemical_symbols(["B" if s == "Be" else s for s in atoms.get_chemical_symbols()])

    t0 = time.time()
    try:
        relaxed = pre.relax(atoms, pressure_gpa=p_target)
        t_pre = time.time() - t0
        t1 = time.time()
        r = run_ground_state(relaxed, QESettings(**SETTINGS), runner,
                             str(ROOT / "calib2" / f"{label}_{int(p_target)}"), relax=False)
        row = dict(compound=label, mace_target_gpa=p_target,
                   prerelax_s=round(t_pre, 1), qe_s=round(time.time() - t1, 1),
                   grid="x".join(map(str, r.grid_shape)),
                   dft_pressure_gpa=round(r.pressure_gpa, 1), fermi_ev=round(r.fermi_ev, 3),
                   phi=round(r.phi, 4), phi_star=round(r.phi_star, 4),
                   h_dos=round(r.h_dos, 4),
                   tc_k=round(belli2025_tc(r.phi, r.phi_star, r.h_f, r.h_dos), 1))
    except Exception as exc:
        row = dict(compound=label, mace_target_gpa=p_target, error=str(exc)[:220])
    rows.append(row)
    print(json.dumps(row), flush=True)

(ROOT / "calib2_results.json").write_text(json.dumps(rows, indent=2))
print("done")
