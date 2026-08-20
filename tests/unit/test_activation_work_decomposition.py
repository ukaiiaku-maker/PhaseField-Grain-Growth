import numpy as np

from grain_growth_pf.disconnections.mode import ModeDriving
from grain_growth_pf.disconnections.spectrum import isotropic_surrogate_library


def test_non_easy_modes_carry_signed_activation_vacancy_count():
    modes = isotropic_surrogate_library(
        b_shells=(0.25, 0.5), directions=4, step_heights=(0.25,),
        barrier_core_ev=0.2, b_coefficient_ev=0.01, h_coefficient_ev=0.01,
    )
    activated = [mode for mode in modes if mode.family != "easy"]
    assert activated
    assert all(mode.activation_vacancies == mode.point_defect_quota for mode in activated)
    assert any(mode.activation_vacancies > 0 for mode in activated)
    assert any(mode.activation_vacancies < 0 for mode in activated)


def test_vacancy_chemical_potential_biases_opposite_climb_modes_oppositely():
    modes = isotropic_surrogate_library(
        b_shells=(0.25, 0.5), directions=4, step_heights=(0.25,),
        barrier_core_ev=0.2, b_coefficient_ev=0.01, h_coefficient_ev=0.01,
        attempt_frequency=1e5,
    )
    plus = next(mode for mode in modes if mode.family != "easy" and mode.activation_vacancies > 0)
    minus = next(mode for mode in modes if mode.family != "easy" and mode.activation_vacancies < 0)
    drive = ModeDriving(vacancy_chemical_potential=0.05)
    assert plus.activation_work_ev(drive) > 0
    assert minus.activation_work_ev(drive) < 0
    assert plus.rate(900.0, drive) > plus.rate(900.0, ModeDriving())
    assert minus.rate(900.0, drive) < minus.rate(900.0, ModeDriving())
    assert np.isfinite(plus.rate(900.0, drive))
