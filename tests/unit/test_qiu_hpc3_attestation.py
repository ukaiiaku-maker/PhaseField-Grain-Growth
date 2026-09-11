import pytest

from attest_qiu_hpc3_source import extract_source_commit


def test_extract_source_commit_from_seed_and_nonseed_entrypoints():
    commit = "a" * 40
    assert extract_source_commit(["bash", "run.sh", commit, "b" * 64]) == commit
    assert extract_source_commit(["bash", "run.sh", "5102", commit, "b" * 64]) == commit


def test_extract_source_commit_rejects_ambiguous_entrypoint():
    with pytest.raises(ValueError):
        extract_source_commit(["bash", "run.sh", "a" * 40, "b" * 40])
