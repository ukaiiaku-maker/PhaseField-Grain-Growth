import pandas as pd

from grain_growth_pf.analysis.aggregate import aggregate_summaries
from grain_growth_pf.analysis.campaign import SUMMARY_COLUMNS


def _row(regime: str, temperature: float, coefficient: float) -> dict:
    row = {column: 0 for column in SUMMARY_COLUMNS}
    row.update({
        "regime": regime, "temperature": temperature, "K": coefficient,
        "Git_SHA": "simulation-sha", "number_of_realizations": 5,
    })
    return row


def test_aggregate_summaries_prefers_later_matched_temperature_series(tmp_path):
    mechanism = tmp_path / "mechanism.csv"
    temperature = tmp_path / "temperature.csv"
    pd.DataFrame([_row("G2", 900.0, 1.0), _row("C1", 900.0, 2.0)]).to_csv(
        mechanism, index=False
    )
    pd.DataFrame([_row("G2", 800.0, 0.5), _row("G2", 900.0, 1.5)]).to_csv(
        temperature, index=False
    )
    output = aggregate_summaries(
        [mechanism, temperature], tmp_path / "final.csv"
    )
    result = pd.read_csv(output)
    assert len(result) == 3
    assert result[(result.regime == "G2") & (result.temperature == 900.0)].K.iloc[0] == 1.5
