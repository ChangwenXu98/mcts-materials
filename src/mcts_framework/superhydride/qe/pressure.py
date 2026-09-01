"""
Matching the DFT pressure without a variable-cell relaxation.

Why this exists
---------------
A MACE pre-relaxation gives a cell quickly, but the potential is fitted on
ambient-pressure data and its error under compression is strongly
compound-dependent - measured here, at a nominal target of 70 GPa, LaBeH8
landed near 160 GPa and YBeH8 near 300 GPa. Ranking those two against each
other would be comparing compression, not chemistry: phi, phi* and H_DOS all
move steeply with pressure.

A full `vc-relax` fixes that and costs an order of magnitude more. This is the
cheap middle: rescale the cell isotropically, re-run the SCF, and use the
stress it prints to steer. Three SCFs land within a few GPa of the target, and
each one is a step the funnel needed anyway - the accepted SCF is the one whose
density the ELF and the projected DOS are then read from.

What it does not do
-------------------
**The internal coordinates are not re-optimised.** Isotropic scaling preserves
fractional positions, so the free Wyckoff parameters keep the values the
pre-relaxation gave them. For a high-symmetry template whose internal
coordinate moves slowly with pressure this is a small error; it is not a
substitute for relaxing the ions, and a candidate being taken seriously wants a
real two-pass vc-relax.

(c) 2026. Triad National Security, LLC. All rights reserved.
"""

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING, Optional

if TYPE_CHECKING:
    from ase import Atoms


@dataclass(frozen=True)
class PressureMatch:
    """
    Drive the DFT pressure onto a target by isotropic cell scaling.

    Args:
        target_gpa: the pressure the campaign is defined at.
        tolerance_gpa: stop once the SCF stress is this close.
        max_scf: hard cap on SCF evaluations, including the first. 3 is
            normally enough: one to measure, one to correct with an assumed
            bulk modulus, one to correct with the measured one.
        bulk_modulus_guess_gpa: starting stiffness, used only for the first
            correction. Superhydrides above 100 GPa sit in the few-hundred-GPa
            range; the value stops mattering once two points are in hand.
        max_volume_step: largest single-step fractional change in volume, so a
            bad stiffness estimate cannot throw the cell across the map.
    """

    target_gpa: float
    tolerance_gpa: float = 10.0
    max_scf: int = 3
    bulk_modulus_guess_gpa: float = 400.0
    max_volume_step: float = 0.12

    def within_tolerance(self, pressure_gpa: float) -> bool:
        return abs(pressure_gpa - self.target_gpa) <= self.tolerance_gpa

    def next_volume(
        self,
        volume: float,
        pressure_gpa: float,
        previous: Optional[tuple] = None,
    ) -> float:
        """
        Volume to try next, from the current (V, P) and optionally the last one.

        With two points the bulk modulus is measured rather than assumed:
        B = -dP / dlnV. A nonsensical estimate (negative, or outside the range
        a compressed solid can plausibly have) is discarded in favour of the
        guess rather than trusted.
        """
        bulk = self.bulk_modulus_guess_gpa
        if previous is not None:
            volume_prev, pressure_prev = previous
            log_ratio = math.log(volume / volume_prev)
            if abs(log_ratio) > 1e-9:
                measured = -(pressure_gpa - pressure_prev) / log_ratio
                if 50.0 <= measured <= 3000.0:
                    bulk = measured

        log_step = (pressure_gpa - self.target_gpa) / bulk
        cap = math.log(1.0 + self.max_volume_step)
        log_step = max(-cap, min(cap, log_step))
        return volume * math.exp(log_step)


def scale_to_volume(atoms: "Atoms", volume: float) -> "Atoms":
    """
    Return a copy of ``atoms`` scaled isotropically to ``volume``.

    Fractional coordinates are preserved, so this changes the compression and
    nothing else about the structure - in particular it cannot break the
    symmetry the template was built with.
    """
    scaled = atoms.copy()
    factor = (volume / atoms.get_volume()) ** (1.0 / 3.0)
    scaled.set_cell(atoms.get_cell() * factor, scale_atoms=True)
    return scaled
