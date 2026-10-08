#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path


EXPECTED_MODULES = {
    "B0": set(),
    "G": {"gb_compatibility", "multihit_persistent"},
    "T": {"tj_compatibility", "multihit_persistent"},
    "GT": {"gb_compatibility", "tj_compatibility", "multihit_persistent"},
    "GTC_GB": {
        "gb_compatibility", "tj_compatibility", "multihit_persistent",
        "area_loss_climb", "gb_defect_sink",
    },
    "GSC_GBTJ_Ks025": {
        "gb_compatibility", "multihit_persistent", "shear_memory",
        "area_loss_climb", "gb_defect_sink", "tj_defect_sink",
    },
    "GTSC_GBTJ_Ks025": {
        "gb_compatibility", "tj_compatibility", "multihit_persistent",
        "shear_memory", "area_loss_climb", "gb_defect_sink", "tj_defect_sink",
    },
    "GTSC_GBTJ_Ks030_LONG": {
        "gb_compatibility", "tj_compatibility", "multihit_persistent",
        "shear_memory", "area_loss_climb", "gb_defect_sink", "tj_defect_sink",
    },
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run")
    parser.add_argument("regime", choices=sorted(EXPECTED_MODULES))
    parser.add_argument("--minimum-steps", type=int, default=1)
    args = parser.parse_args()

    run = Path(args.run)
    manifest = json.loads((run / "manifest.json").read_text())
    checkpoint = json.loads((run / "checkpoint.json").read_text())
    config = manifest["config"]
    pf = config["pf"]
    modules = set(config["active_modules"])
    assertions = {
        "regime": config["regime"] == args.regime,
        "module_set": modules == EXPECTED_MODULES[args.regime],
        "anisotropy_strength": pf["anisotropy_strength"] == "A2_STRONG",
        "anisotropic_energy": pf["anisotropic_energy"] is True,
        "anisotropic_mobility": pf["anisotropic_mobility"] is True,
        "compact_support": pf["anisotropic_support_mode"] == "compact_active_set",
        "kkt_tolerance": float(pf["anisotropic_kkt_tolerance"]) == 1e-10,
        "matched_shape": pf["shape"] == [384, 384],
        "matched_temperature": float(pf["temperature"]) == 900.0,
        "physical_horizon": float(config["parameters"]["maximum_physical_time"]) == 4000.0,
        "accepted_steps": int(checkpoint["step_number"]) >= args.minimum_steps,
        "checkpoint_time_positive": float(checkpoint["time"]) > 0.0,
        "event_ledger_declared": bool(manifest.get("event_ledger")),
        "checkpoint_domain_state": "domains" in checkpoint,
    }
    if "multihit_persistent" in modules:
        domains = checkpoint.get("domains", {})
        assertions["stochastic_domain_rng_state"] = bool(domains) and all(
            "rng" in state for state in domains.values()
        )
    if "area_loss_climb" in modules:
        assertions["climb_modules_complete"] = {
            "gb_defect_sink", "area_loss_climb"
        }.issubset(modules)
    if "shear_memory" in modules:
        assertions["shear_backend"] = config["mechanics_backend"] == "local_memory"
    failed = sorted(name for name, passed in assertions.items() if not passed)
    report = {
        "classification": "PASSED" if not failed else "FAILED",
        "regime": args.regime,
        "active_modules": sorted(modules),
        "assertions": assertions,
        "failed": failed,
        "steps": int(checkpoint["step_number"]),
        "physical_time": float(checkpoint["time"]),
    }
    print(json.dumps(report, indent=2, sort_keys=True))
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
