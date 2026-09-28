#!/usr/bin/env python3
"""Freeze A2 normalization on a promoted Qiu four-reference initial network."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

from grain_growth_pf.mechanics.anisotropy import BoundaryLaw, LADDER, angular_normalization

ARCHIVE_SHA256 = "2740e0af26bd3acd24bfbbca44b8b3cd56e2e3bacf97ef5648ab69b84989bf90"


def sha256(path: Path) -> str:
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def scalar_text(value: np.ndarray) -> str:
    return str(np.asarray(value).reshape(()).item())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("network", type=Path, help="promoted initial-network NPZ")
    parser.add_argument("--template", type=Path, default=Path("configs/production/qiu_si_4ref_controls.json"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    with np.load(args.network, allow_pickle=False) as data:
        required = {
            "theta", "orientation_i", "orientation_j", "length",
            "model_identity", "archive_sha256", "native_promoted",
            "orientation_reference_sha256",
        }
        missing = required - set(data.files)
        if missing:
            raise SystemExit(f"missing normalization evidence: {sorted(missing)}")
        if scalar_text(data["model_identity"]) != "QIU_SI_REFERENCE":
            raise SystemExit("normalization input is not QIU_SI_REFERENCE")
        if scalar_text(data["archive_sha256"]) != ARCHIVE_SHA256:
            raise SystemExit("normalization archive identity mismatch")
        if not bool(np.asarray(data["native_promoted"]).reshape(())):
            raise SystemExit("native baseline has not passed its promotion gate")
        theta = data["theta"]
        orientation_i = data["orientation_i"]
        orientation_j = data["orientation_j"]
        length = data["length"]
        orientation_reference_sha256 = scalar_text(data["orientation_reference_sha256"])

    law = BoundaryLaw(strength=LADDER["A2_STRONG"], mobility0=QIU_NATIVE_MOBILITY).normalize(
        theta, orientation_i, orientation_j, length
    )
    gamma, _, _, stiffness, mobility = law.evaluate(theta, orientation_i, orientation_j)
    if np.min(gamma) <= 0 or np.min(stiffness) <= 0:
        raise SystemExit("A2 is nonpositive on the native initial network")

    manifest = json.loads(args.template.read_text())
    manifest["a2_candidate"]["energy_normalization"] = law.energy_normalization
    manifest["a2_candidate"]["mobility_normalization"] = law.mobility_normalization
    manifest["a2_candidate"]["angular_scale"] = angular_normalization(LADDER["A2_STRONG"])
    network_sha256 = sha256(args.network)
    manifest["a2_candidate"]["normalization_source"] = {
        "path": str(args.network), "sha256": network_sha256,
        "model_identity": "QIU_SI_REFERENCE", "native_promoted": True,
        "orientation_reference_sha256": orientation_reference_sha256,
    }
    mobility_stiffness = mobility * stiffness
    manifest["a2_candidate"]["sampled_distributions"] = {
        "gamma_min": float(np.min(gamma)),
        "gamma_mean_length_weighted": float(np.average(gamma, weights=length)),
        "gamma_p05": float(np.percentile(gamma, 5)),
        "gamma_p95": float(np.percentile(gamma, 95)),
        "stiffness_min": float(np.min(stiffness)),
        "stiffness_p05": float(np.percentile(stiffness, 5)),
        "stiffness_p95": float(np.percentile(stiffness, 95)),
        "mobility_mean_length_weighted": float(np.average(mobility, weights=length)),
        "mobility_p05": float(np.percentile(mobility, 5)),
        "mobility_p95": float(np.percentile(mobility, 95)),
        "mobility_times_stiffness_p05": float(np.percentile(mobility_stiffness, 5)),
        "mobility_times_stiffness_p95": float(np.percentile(mobility_stiffness, 95)),
        "gamma_mobility_pearson": float(np.corrcoef(gamma, mobility)[0, 1]),
        "gamma_mobility_rank_correlation": float(spearmanr(gamma, mobility).statistic),
        "sample_count": int(gamma.size),
    }
    manifest["a2_candidate"]["native_network_sha256"] = network_sha256
    manifest["a2_candidate"]["orientation_reference_sha256"] = orientation_reference_sha256
    manifest["release_gates"]["native_promoted"] = True
    manifest["release_gates"]["normalization_frozen"] = True
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


QIU_NATIVE_MOBILITY = 2.0 * np.pi**2 / (8.0 * 5.0)


if __name__ == "__main__":
    main()
