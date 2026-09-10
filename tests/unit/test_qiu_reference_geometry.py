import ast
from pathlib import Path

import numpy as np
import pytest

from grain_growth_pf.mechanics.qiu_reference_geometry import (
    qiu_reference_elastic_projection,
    resolved_reference_shear,
    two_reference_coupling_factors,
    two_reference_decomposition,
)


ARCHIVED = Path(
    "/Users/sdillon/PF-graingrowth/.external/qiu/PF_Codes/Polycrystals/functions_2ref.py"
)


def _load_archived_beta():
    if not ARCHIVED.exists():
        pytest.skip("local archived Qiu source is unavailable")
    tree = ast.parse(ARCHIVED.read_text())
    function = next(
        node for node in tree.body
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == "beta"
    )
    function.decorator_list = []
    module = ast.Module(body=[function], type_ignores=[])
    namespace = {"np": np}
    exec(compile(ast.fix_missing_locations(module), str(ARCHIVED), "exec"), namespace)
    return namespace["beta"]


def test_two_reference_beta_matches_archived_function_directly():
    archived = _load_archived_beta()
    pairs = [
        (0.0, 0.2), (0.3, 1.4), (2.8, 0.1),
        (np.pi, 0.0), (0.0, np.pi), (1.1, 1.1),
        (np.deg2rad(70), 0.0), (np.deg2rad(150), 0.0),
    ]
    for left, right in pairs:
        assert np.allclose(
            two_reference_coupling_factors(left, right),
            archived(left, right), rtol=0.0, atol=2e-15,
        )


def test_two_reference_factor_order_reverses_with_grain_order():
    forward = two_reference_coupling_factors(0.2, 1.1)
    reverse = two_reference_coupling_factors(1.1, 0.2)
    assert np.allclose(reverse, (-forward[0], -forward[1]), atol=2e-15)


def test_reference_selection_switches_sectors_without_changing_physical_sum():
    stress = np.asarray([[0.6, -0.3], [-0.3, 1.2]])
    first = two_reference_decomposition(0.1, 0.7, np.asarray([1.0, 1.0]))
    second = two_reference_decomposition(0.1, 0.7, np.asarray([1.0, -1.0]))
    for decomposition in (first, second):
        direct = -sum((
            resolved_reference_shear(stress, decomposition.reference_angle_1)
            * decomposition.coupling_1,
            resolved_reference_shear(stress, decomposition.reference_angle_2)
            * decomposition.coupling_2,
        ))
        assert qiu_reference_elastic_projection(stress, decomposition) == direct
    assert np.isclose(abs(first.reference_angle_2 - first.reference_angle_1), np.pi / 2)
    assert np.isclose(abs(second.reference_angle_1 - second.reference_angle_2), np.pi / 2)


def test_reference_projection_sign_reverses_with_stress_and_grain_order():
    stress = np.asarray([[0.4, 0.2], [0.2, -0.1]])
    normal = np.asarray([0.3, 0.8])
    forward = two_reference_decomposition(0.2, 0.9, normal)
    reverse = two_reference_decomposition(0.9, 0.2, -normal)
    value = qiu_reference_elastic_projection(stress, forward)
    assert qiu_reference_elastic_projection(-stress, forward) == -value
    assert np.isclose(qiu_reference_elastic_projection(stress, reverse), -value)
