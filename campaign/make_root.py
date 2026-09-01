"""Build the LaBeH8 root template (LaBH8 structure type, Fm-3m).

Geometry is fixed by symmetry apart from two numbers: the cubic lattice
constant and the H 32f coordinate. Both are refined by the MACE pre-relaxation
at the target pressure, so the starting values below only have to be close
enough for the optimiser, not right.

    La   4a  (0, 0, 0)
    Be   4b  (1/2, 1/2, 1/2)
    H   32f  (x, x, x)          x ~ 0.372 -> an H8 cube around Be

The primitive fcc cell holds 10 atoms: 1 La + 1 Be + 8 H, H_f = 0.8.
"""
import sys
from ase.io import write
from ase.spacegroup import crystal

A_START = 5.10        # angstrom, cubic lattice constant before relaxation
X_32F = 0.372         # H 32f coordinate

atoms = crystal(
    symbols=["La", "Be", "H"],
    basis=[(0.0, 0.0, 0.0), (0.5, 0.5, 0.5), (X_32F, X_32F, X_32F)],
    spacegroup=225,
    cellpar=[A_START, A_START, A_START, 90, 90, 90],
    primitive_cell=True,
)

print("formula      ", atoms.get_chemical_formula(mode="metal"))
print("natoms       ", len(atoms))
print("volume       ", round(atoms.get_volume(), 3), "A^3")
h = [i for i, s in enumerate(atoms.get_chemical_symbols()) if s == "H"]
be = [i for i, s in enumerate(atoms.get_chemical_symbols()) if s == "Be"][0]
la = [i for i, s in enumerate(atoms.get_chemical_symbols()) if s == "La"][0]
print("n_H          ", len(h), " H_f =", round(len(h) / len(atoms), 4))
print("Be-H (min)   ", round(min(atoms.get_distance(be, i, mic=True) for i in h), 4), "A")
print("La-H (min)   ", round(min(atoms.get_distance(la, i, mic=True) for i in h), 4), "A")
print("H-H  (min)   ", round(min(atoms.get_distance(i, j, mic=True)
                                for k, i in enumerate(h) for j in h[k+1:]), 4), "A")

if len(sys.argv) > 1:
    write(sys.argv[1], atoms)
    print("wrote", sys.argv[1])
