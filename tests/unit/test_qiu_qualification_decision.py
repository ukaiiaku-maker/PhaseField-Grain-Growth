from validate_qiu_qualification_decision import REQUIRED_RUN_ROLES, validate


def _artifact():
    return {"path": "/absolute/artifact", "sha256": "d" * 64}


def _decision():
    return {
        "schema_version": 1,
        "created_at": "2026-09-11T00:00:00Z",
        "source_branch": "codex/qiu-full-field-qualification-v1",
        "source_commit": "a" * 40,
        "historical_source_commit": "b" * 40,
        "historical_run_path": "/immutable/history",
        "initial_state_hashes": {
            str(seed): "c" * 64 for seed in (5101, 5102, 5103)
        },
        "selected_backend": "FFT_EIGENSTRAIN_V2",
        "model_identity": "corrected FFT surrogate",
        "classification": "QIU_REFERENCE_MODEL_MISIDENTIFIED",
        "root_causes": ["legacy operator defect"],
        "rejected_hypotheses": ["extinction threshold"],
        "mathematical_gates": {"equilibrium": True},
        "source_mapping_gates": {"work": True},
        "timestep_gates": {"matched_time": True},
        "production_runs": [
            {
                "role": role,
                "pooling_allowed": True,
                "terminal_status": "completed",
                "terminal_step": 100,
                "source_commit": "a" * 40,
                "archive_sha256": "c" * 64,
                "local_path": f"/absolute/{role}",
            }
            for role in sorted(REQUIRED_RUN_ROLES)
        ],
        "seed_results": [
            {"seed": seed, "terminal": True} for seed in (5101, 5102, 5103)
        ],
        "excluded_runs": [{"path": "/excluded/run", "reason": "superseded"}],
        "tests": {
            "total": 1, "failures": 0, "errors": 0, "skipped": 0,
            "junit": _artifact(),
        },
        "commit_shas": ["a" * 40],
        "remaining_limitations": ["2-D surrogate"],
        "artifacts": {
            "report": _artifact(),
            "run_matrix": _artifact(),
            "analysis_summary": _artifact(),
            "plots": [_artifact()],
            "movies": [
                {**_artifact(), "role": role}
                for role in (
                    "legacy", "corrected_seed5101", "corrected_additional_seed"
                )
            ],
        },
        "canonical_preservation": {
            "verified_unchanged": True,
            "verification_manifest": _artifact(),
        },
        "pull_request": {
            "state": "OPEN", "merged": False,
            "url": "https://github.com/example/repository/pull/1",
        },
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


def test_failed_gate_or_missing_movie_role_is_rejected():
    decision = _decision()
    decision["mathematical_gates"]["equilibrium"] = False
    decision["artifacts"]["movies"] = decision["artifacts"]["movies"][:-1]

    errors = validate(decision)

    assert any("equilibrium is not recorded as passed" in error for error in errors)
    assert any("movies must contain at least 3" in error for error in errors)


def test_missing_canonical_proof_or_open_pr_is_rejected():
    decision = _decision()
    decision["canonical_preservation"]["verified_unchanged"] = False
    decision["pull_request"]["state"] = "MERGED"
    decision["pull_request"]["merged"] = True

    errors = validate(decision)

    assert any("verified_unchanged=true" in error for error in errors)
    assert any("OPEN or DRAFT" in error for error in errors)
    assert any("merged=false" in error for error in errors)
