"""
Unit tests for matching the DFT pressure by isotropic cell scaling.

(c) 2026. Triad National Security, LLC. All rights reserved.
"""

import math

import pytest

from mcts_framework.superhydride.qe.pressure import PressureMatch, scale_to_volume

TARGET = 150.0


@pytest.fixture
def match():
    return PressureMatch(target_gpa=TARGET, tolerance_gpa=10.0, max_scf=3)


# --- Acceptance -----------------------------------------------------------


def test_tolerance_is_two_sided(match):
    assert match.within_tolerance(TARGET)
    assert match.within_tolerance(TARGET + 9.9)
    assert match.within_tolerance(TARGET - 9.9)
    assert not match.within_tolerance(TARGET + 10.1)
    assert not match.within_tolerance(TARGET - 10.1)


# --- Stepping -------------------------------------------------------------


def test_overpressure_expands_and_underpressure_compresses(match):
    assert match.next_volume(34.0, 300.0) > 34.0   # too compressed -> grow
    assert match.next_volume(34.0, 50.0) < 34.0    # too loose -> shrink


def test_a_step_uses_the_guessed_stiffness_when_it_has_one_point(match):
    """dlnV = (P - P_target) / B, with B the guess until a second point exists."""
    # 180 GPa is 30 above target, an 7.8% step - inside the volume cap, so this
    # tests the stiffness rather than the clamp.
    expected = 34.0 * math.exp((180.0 - TARGET) / match.bulk_modulus_guess_gpa)
    assert expected < 34.0 * (1.0 + match.max_volume_step)
    assert match.next_volume(34.0, 180.0) == pytest.approx(expected)


def test_a_second_point_replaces_the_guess_with_a_measurement(match):
    """B = -dP/dlnV, measured from the two SCFs actually run."""
    v_prev, p_prev, v, p = 34.0, 300.0, 38.0, 210.0
    measured = -(p - p_prev) / math.log(v / v_prev)
    assert 50.0 <= measured <= 3000.0
    expected = v * math.exp((p - TARGET) / measured)
    assert match.next_volume(v, p, (v_prev, p_prev)) == pytest.approx(expected)


def test_an_implausible_stiffness_is_discarded(match):
    """
    A negative or absurd B - from SCF noise, or two points too close together -
    would throw the cell across the map. The guess is used instead.
    """
    # Pressure rising with volume gives a negative bulk modulus.
    with_bad = match.next_volume(38.0, 320.0, (34.0, 300.0))
    with_guess = match.next_volume(38.0, 320.0)
    assert with_bad == pytest.approx(with_guess)


def test_identical_volumes_do_not_divide_by_zero(match):
    assert match.next_volume(34.0, 200.0, (34.0, 199.0)) == pytest.approx(
        match.next_volume(34.0, 200.0)
    )


def test_a_single_step_is_capped(match):
    """A wild first pressure must not send the cell somewhere unphysical."""
    huge = match.next_volume(34.0, 5000.0)
    assert huge == pytest.approx(34.0 * (1.0 + match.max_volume_step))
    tiny = match.next_volume(34.0, -5000.0)
    assert tiny == pytest.approx(34.0 / (1.0 + match.max_volume_step))


def test_repeated_steps_converge_on_the_target():
    """
    Against a synthetic solid with a known bulk modulus, the secant update
    reaches the tolerance within the SCF budget.
    """
    match = PressureMatch(target_gpa=TARGET, tolerance_gpa=5.0, max_scf=3)
    bulk, v_ref, p_ref = 500.0, 34.0, 300.0

    def pressure(volume):  # P = P_ref - B ln(V/V_ref)
        return p_ref - bulk * math.log(volume / v_ref)

    volume, previous = v_ref, None
    for _ in range(match.max_scf):
        p = pressure(volume)
        if match.within_tolerance(p):
            break
        nxt = match.next_volume(volume, p, previous)
        previous, volume = (volume, p), nxt
    assert match.within_tolerance(pressure(volume))


# --- Scaling --------------------------------------------------------------


def test_scaling_hits_the_requested_volume(make_superhydride_structure):
    pytest.importorskip("ase")
    atoms = make_superhydride_structure().atoms
    scaled = scale_to_volume(atoms, 40.0)
    assert scaled.get_volume() == pytest.approx(40.0)


def test_scaling_preserves_fractional_coordinates(make_superhydride_structure):
    """
    Isotropic scaling changes the compression and nothing else - in particular
    it cannot break the symmetry the template was built with.
    """
    pytest.importorskip("ase")
    atoms = make_superhydride_structure().atoms
    scaled = scale_to_volume(atoms, 55.0)
    assert scaled.get_scaled_positions() == pytest.approx(atoms.get_scaled_positions())
    assert scaled.get_chemical_symbols() == atoms.get_chemical_symbols()


def test_scaling_does_not_mutate_the_input(make_superhydride_structure):
    pytest.importorskip("ase")
    atoms = make_superhydride_structure().atoms
    before = atoms.get_volume()
    scale_to_volume(atoms, 60.0)
    assert atoms.get_volume() == pytest.approx(before)


def test_scaling_preserves_the_structure_identity(make_superhydride_structure):
    """The space group and Wyckoff decoration must survive a pressure step."""
    pytest.importorskip("spglib")
    from mcts_framework.superhydride import SuperhydrideStructure

    material = make_superhydride_structure()
    scaled = SuperhydrideStructure(scale_to_volume(material.atoms, 48.0))
    assert scaled.get_identifier() == material.get_identifier()
