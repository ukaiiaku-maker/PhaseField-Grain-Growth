import pandas as pd

from audit_qiu_sourcework_stream import reconcile_streams


def test_reconcile_streams_sums_boundary_work_and_sweep():
    scalar = pd.DataFrame({
        "step": [10, 11],
        "predicted_source_work": [3.0, -4.0],
        "integrated_boundary_sweep": [7.0, 2.0],
    })
    boundary = pd.DataFrame({
        "step": [10, 10, 11, 11],
        "predicted_elastic_work": [1.0, 2.0, -1.5, -2.5],
        "integrated_sweep": [3.0, 4.0, 1.0, 1.0],
    })
    merged, result = reconcile_streams(scalar, boundary, transition_step=11)
    assert merged["boundary_predicted_work"].tolist() == [3.0, -4.0]
    assert result["all_steps"]["maximum_abs_predicted_work_difference"] == 0.0
    assert result["all_steps"]["maximum_abs_integrated_sweep_difference"] == 0.0
