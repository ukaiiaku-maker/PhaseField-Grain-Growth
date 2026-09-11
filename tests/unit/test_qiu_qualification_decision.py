from validate_qiu_qualification_decision import REQUIRED_RUN_ROLES, validate


def _decision():
    return {
        "schema_version": 1,
        "created_at": "2026-09-11T00:00:00Z",
        "source_branch": "codex/qiu-full-field-qualification-v1",
        "source_commit": "a" * 40,
        "historical_source_commit": "b" * 40,
        "historical_run_path": "/immutable/history",
        "initial_state_hashes": {"5101": "c" * 64},
        "selected_backend": "FFT_EIGENSTRAIN_V2",
        "model_identity": "corrected FFT surrogate",
        "classification": "QIU_REFERENCE_MODEL_MISIDENTIFIED",
        "root_causes": [],
        "rejected_hypotheses": [],
        "mathematical_gates": {"equilibrium": True},
        "source_mapping_gates": {"work": True},
        "timestep_gates": {"matched_time": True},
        "production_runs": [
            {"role": role, "pooling_allowed": True, "terminal_status": "completed"}
            for role in sorted(REQUIRED_RUN_ROLES)
        ],
        "seed_results": [],
        "excluded_runs": [],
        "tests": {"total": 1, "failures": 0, "errors": 0, "skipped": 0},
        "commit_shas": [],
        "remaining_limitations": [],
    }


def test_complete_terminal_decision_is_valid():
    assert validate(_decision()) == []


def test_nonterminal_or_missing_seed_decision_is_rejected():
    decision = _decision()
    decision["production_runs"] = decision["production_runs"][:-1]
    decision["tests"]["failures"] = 1

    errors = validate(decision)

    assert any("missing production run roles" in error for error in errors)
    assert any("zero failures" in error for error in errors)
