"""
Cheap pre-relaxation with a machine-learning interatomic potential.

A screening campaign cannot afford a `vc-relax` per candidate: two passes of
variable-cell DFT dominate the cost of the descriptor funnel by an order of
magnitude. Relaxing with MACE first and then running a single self-consistent
DFT step on that cell buys most of the geometry for a fraction of the time.

What this costs, stated plainly
-------------------------------
**The cell is not variationally relaxed at the DFT level.** MACE-MP is fitted
on Materials Project data, which is overwhelmingly ambient-pressure, so at
100+ GPa it is extrapolating. The cell it returns is a reasonable geometry, not
the DFT equilibrium one at the target pressure, and the difference shows up as
a residual stress in the SCF that follows.

That residual is measured rather than assumed: the SCF prints its stress, the
funnel records it as ``pressure_gpa``, and a candidate whose recorded pressure
is far from the target had its descriptors computed on a cell that does not
describe the pressure claimed. Screen on that column before believing a Tc.

For ranking candidates within one campaign - all pre-relaxed the same way, all
evaluated the same way - this is a reasonable trade. For quoting a Tc against
the literature it is not, and those candidates want a real two-pass vc-relax.

(c) 2026. Triad National Security, LLC. All rights reserved.
"""

import logging
from dataclasses import dataclass
from typing import TYPE_CHECKING, Optional

if TYPE_CHECKING:
    from ase import Atoms

logger = logging.getLogger(__name__)

#: 1 GPa in eV/Angstrom^3, the unit ASE's cell filters want for a pressure.
GPA_TO_EV_PER_ANG3 = 1.0 / 160.21766208


class PrerelaxError(RuntimeError):
    """A pre-relaxation that could not be run or did not converge."""


@dataclass
class MacePrerelax:
    """
    Relax a structure at constant pressure with a MACE-MP foundation model.

    Args:
        model: MACE-MP model size ('small', 'medium', 'large'). The weights are
            fetched to ~/.cache/mace on first use, so pre-warm the cache before
            submitting to compute nodes without outbound network.
        device: 'cpu' or 'cuda'. CPU is the right choice inside a DFT campaign -
            the relaxation is seconds and the node's cores are for pw.x.
        fmax: force convergence, eV/Angstrom.
        max_steps: optimiser step cap. Hitting it is reported, not ignored.
        dtype: MACE default dtype. float64 matches the calculator's own default
            for energy/force evaluation used elsewhere in this package.
    """

    model: str = "medium"
    device: str = "cpu"
    fmax: float = 0.05
    max_steps: int = 300
    dtype: str = "float64"
    #: Raise when the optimiser hits max_steps instead of returning its best
    #: geometry. Default False, and the reason is a search bias rather than a
    #: convenience: a candidate that fails here scores 0.0, which the search
    #: reads as "bad material" when it means "the optimiser gave up". That
    #: teaches the search to avoid regions where MACE struggles rather than
    #: regions where the physics is poor. In the campaign that surfaced this,
    #: 7 of 8 early failures were this and nothing else.
    #:
    #: Returning the unconverged geometry is safe here because a DFT step
    #: follows and measures it: with a pressure match the cell is corrected
    #: outright, and the recorded pressure says how far it landed. The caller
    #: is told via `last_relax_converged`.
    strict: bool = False

    _calculator: object = None
    #: Whether the most recent relax() reached fmax. False after a run that
    #: exhausted max_steps and was returned anyway.
    last_relax_converged: bool = True

    def _get_calculator(self):
        """Build the MACE calculator once and keep it; construction dominates."""
        if self._calculator is None:
            try:
                from mace.calculators import mace_mp
            except ImportError as exc:  # pragma: no cover - optional dependency
                raise PrerelaxError(
                    "MACE pre-relaxation needs mace-torch. Install it, or turn "
                    "the pre-relaxation off."
                ) from exc
            logger.info("Loading MACE-MP '%s' on %s", self.model, self.device)
            self._calculator = mace_mp(
                model=self.model,
                device=self.device,
                default_dtype=self.dtype,
                dispersion=False,
            )
        return self._calculator

    def relax(self, atoms: "Atoms", pressure_gpa: Optional[float] = None) -> "Atoms":
        """
        Return a relaxed copy of ``atoms``.

        Args:
            atoms: the structure to relax.
            pressure_gpa: external pressure. None relaxes the ions only, at
                fixed cell; a value relaxes cell and ions against that
                pressure.

        Returns:
            A new Atoms with the relaxed cell and positions, and no calculator
            attached (so it stays cheap to copy and safe to hand to the QE
            input writer).

        Raises:
            PrerelaxError: if MACE is unavailable, or the optimiser hits
                ``max_steps`` and ``strict`` is set.
        """
        from ase.filters import FrechetCellFilter
        from ase.optimize import FIRE

        working = atoms.copy()
        working.calc = self._get_calculator()

        if pressure_gpa is None:
            target = working
        else:
            target = FrechetCellFilter(
                working, scalar_pressure=pressure_gpa * GPA_TO_EV_PER_ANG3
            )

        optimizer = FIRE(target, logfile=None)
        converged = bool(optimizer.run(fmax=self.fmax, steps=self.max_steps))
        self.last_relax_converged = converged
        if not converged:
            message = (
                f"MACE pre-relaxation did not reach fmax={self.fmax} eV/A in "
                f"{self.max_steps} steps for {working.get_chemical_formula()}"
            )
            if self.strict:
                raise PrerelaxError(message)
            logger.warning("%s; using the geometry it reached", message)

        relaxed = working.copy()
        relaxed.calc = None
        logger.info(
            "MACE pre-relax %s: %d steps, volume %.2f -> %.2f A^3",
            working.get_chemical_formula(),
            optimizer.get_number_of_steps(),
            atoms.get_volume(),
            relaxed.get_volume(),
        )
        return relaxed
