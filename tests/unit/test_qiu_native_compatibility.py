import numpy as np
import numba as nb
import pytest

from scripts.qiu_native_compatibility import (
    prepend_junction_rows_compatibility,
    prepend_junction_rows_reference,
)


@pytest.mark.parametrize("boundary,first,last", [
    (np.empty((0, 2), dtype=np.int64), np.array([0, 9]), np.array([7, 2])),
    (np.array([[4, 5]], dtype=np.int64), np.array([0, 1]), np.array([2, 3])),
    (np.array([[4, 5], [9, 1], [3, 8]], dtype=np.int32), np.array([0, 1]), np.array([2, 3])),
    (np.array([[4, 5], [4, 5]], dtype=np.int64), np.array([4, 5]), np.array([4, 5])),
    (np.array([[0, 499], [499, 0]], dtype=np.int64), np.array([499, 499]), np.array([0, 0])),
])
def test_compatibility_assembly_matches_independent_reference(boundary, first, last):
    expected = prepend_junction_rows_reference(boundary, first, last)
    actual = prepend_junction_rows_compatibility(boundary, first, last)
    assert actual.dtype == expected.dtype == boundary.dtype
    assert actual.shape == expected.shape
    assert np.array_equal(actual, expected)


def test_reference_rejects_mismatched_shape():
    with pytest.raises(ValueError):
        prepend_junction_rows_reference(np.empty((0, 2)), np.array([1]), np.array([2]))


def test_compatibility_assembly_compiles_in_nopython_mode():
    compiled = nb.njit(prepend_junction_rows_compatibility)
    boundary = np.array([[3, 4], [5, 6]], dtype=np.int64)
    actual = compiled(boundary, np.array([1, 2]), np.array([7, 8]))
    assert compiled.nopython_signatures
    assert np.array_equal(actual, np.array([[1, 2], [7, 8], [3, 4], [5, 6]]))
