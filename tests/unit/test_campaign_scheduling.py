import json

import yaml

import grain_growth_pf.campaign as campaign_module


def test_campaign_dispatches_one_run_per_pool_chunk(tmp_path, monkeypatch):
    spec_path = tmp_path / "campaign.yaml"
    spec_path.write_text(yaml.safe_dump({
        "regimes": {"B0": {}},
        "temperatures": [900.0],
        "seeds": [1, 2],
        "base": {
            "regime": "B0",
            "pf": {"shape": [8, 8]},
            "parameters": {"initial_grains": 2, "equilibration_steps": 0},
            "max_steps": 1,
            "termination_grains": 1,
        },
    }))
    observed_chunksizes = []

    class FakePool:
        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return False

        def map(self, _function, payloads, chunksize=None):
            observed_chunksizes.append(chunksize)
            return [{"path": payload[1], "status": "completed"} for payload in payloads]

    class FakeContext:
        @staticmethod
        def Pool(_workers):
            return FakePool()

    monkeypatch.setattr(campaign_module.mp, "get_context", lambda _method: FakeContext())
    campaign = campaign_module.launch_campaign(
        spec_path, root=tmp_path / "campaigns", processes=2
    )

    assert observed_chunksizes == [1]
    manifest = json.loads((campaign / "campaign_manifest.json").read_text())
    assert manifest["status"] == "completed"
    assert manifest["workers"] == 2
