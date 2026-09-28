#!/usr/bin/env python3
"""Build the immutable Qiu/FFT lineage audit from verified local evidence."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
from pathlib import Path


REFS = {
    "audited_base": "9f66c8d7a5a266687284d8da35aefbc6062808f7",
    "generic_science": "1316cc89dbabfb41cb883b0d4a4c74738cc2bef6",
    "postfix_recovery": "88fd5dddd3fc425d0b920373cd3698080a3e7452",
    "qiu_si_port": "97e6a9d602792cd2642d35d64396c7513e1a6d4c",
    "qiu_qualification": "cef12338a7befed981052d16a87848c14c920365",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def classify(path: str) -> str:
    if path.startswith("results/") or path.startswith("docs/"):
        return "REPORT_OR_ANALYSIS_ONLY"
    if "fft" in path.lower() or "qiu_full_field" in path.lower():
        return "FFT_SPECIFIC"
    if "qiu_si" in path.lower():
        return "QIU_SI_SPECIFIC"
    if path.startswith("src/grain_growth_pf/pf/") or "anisotrop" in path.lower():
        return "GENERIC_ANISOTROPY_CORE"
    if path.startswith("scripts/") or path.startswith("configs/"):
        return "HPC_ORCHESTRATION_ONLY"
    return "SUPERSEDED"


def changed_files(left: str, right: str) -> list[dict[str, str]]:
    output = subprocess.check_output(
        ["git", "diff", "--name-status", left, right], text=True
    )
    rows = []
    for line in output.splitlines():
        status, path = line.split("\t", 1)
        rows.append({"left": left, "right": right, "status": status,
                     "path": path, "classification": classify(path)})
    return rows


def records(repo_root: Path) -> list[dict[str, object]]:
    qroot = Path("/Users/sdillon/PF-graingrowth/results/qiu_full_field_qualification_20260910")
    legacy1 = qroot / "hpc_terminal/20260910T232437Z-legacy-seed5101/output/qualification/QIU_LEGACY_FORENSIC-T900-s5101"
    legacy2 = qroot / "hpc_terminal/20260911T111913Z-sourcework-legacy-seed5101/output/qualification/QIU_LEGACY_FORENSIC-T900-s5101"
    fft_archive = Path("/Users/sdillon/HPC3/qiu-qualification-20260910/hpc3-results/pfgg-qiu-qualification/20260911T003946Z-nogit-c53868/results/results-single.tar.gz")
    common = {
        "initial_field_sha256": "106819770289d156b0ff6a35d356a9aa172258f428325ec50c5f62f380ea56f6",
        "orientation_reference_sha256": "9817b9d754a021ce9ff313e3e8d7bdf68582a4ea3af4641645be1e77055484c9",
        "grid": [384, 384], "domain": [384.0, 384.0], "accepted_timestep": 0.04,
    }
    return [
        {**common, "model_identity": "QIU_LEGACY_FORENSIC",
         "logical_trajectory_id": "QIU_LEGACY_FORENSIC_seed5101", "segment_index": 0,
         "run_id": "20260910T232437Z-nogit-f4803a", "slurm_job_id": "55930486",
         "source_commit": "147141b54227e7c80c457e9d80c1a2ca3475bf61",
         "source_archive_sha256": "cb9669b98eedd478fe271ac0c67f6d65bbae61cce7c59b28f73b287ba952f492",
         "resolved_configuration_sha256": "f2f744a3b7267cdb81d820bb5c45339d556bfdcff804665d66c2749626f9ae3f",
         "mechanics_backend": "qiu_full_field_legacy_point_eigenstrain", "start_step": 0,
         "end_step": 9984, "start_time": 0.0, "end_time": 399.3600000000567,
         "start_checkpoint_sha256": None, "end_checkpoint_sha256": sha256(legacy1 / "checkpoint.npz"),
         "predecessor_run_id": None, "successor_run_id": None,
         "terminal_state": "COMPLETED", "scientific_classification": "REJECTED_FORENSIC_AVALANCHE",
         "retrieval_state": "VERIFIED"},
        {**common, "model_identity": "QIU_LEGACY_FORENSIC",
         "logical_trajectory_id": "QIU_LEGACY_FORENSIC_seed5101_replay", "segment_index": 0,
         "run_id": "20260911T111913Z-nogit-06fc67", "slurm_job_id": "55948258",
         "source_commit": "bebd53b8706301715c216dafe224f2e0c4180aaa",
         "source_archive_sha256": "d024b52083f7f1c4042d4d8a99987cff8c5c17160499c8d2bbbc2d4d1f8cffb5",
         "resolved_configuration_sha256": "324568b7071d7445889ee44a9b9ad729fa05e8e278ad286fad1fa1926847c2a1",
         "mechanics_backend": "qiu_full_field_legacy_point_eigenstrain", "start_step": 9000,
         "end_step": 9984, "start_time": 360.0, "end_time": 399.3600000000567,
         "start_checkpoint_sha256": "bb7d4bdbd23c336f388ee28872155aea0ef280835bc87588edc8bae0b8c76e20",
         "end_checkpoint_sha256": sha256(legacy2 / "checkpoint.npz"),
         "predecessor_run_id": None, "successor_run_id": None,
         "terminal_state": "COMPLETED", "scientific_classification": "OVERLAPPING_DIAGNOSTIC_REPLAY",
         "retrieval_state": "VERIFIED"},
        {**common, "model_identity": "FFT_EIGENSTRAIN_V2",
         "logical_trajectory_id": "FFT_EIGENSTRAIN_V2_seed5101_dt004", "segment_index": 0,
         "run_id": "20260911T003946Z-nogit-c53868", "slurm_job_id": "55932457",
         "source_commit": "8bb7837677e6ca36fc9a952daf6aef601190b550",
         "source_archive_sha256": "e9926fc38ba41502d4deefadc57b344826924bfd2e4c0d75f4299ee1a6b8d5b0",
         "resolved_configuration_sha256": "ae96a1427ea01c33df8441dfb0cacd1ff4e7cd5c611c9c1326b3a533a3fb19ce",
         "mechanics_backend": "fft_eigenstrain_v2", "start_step": 0, "end_step": 8379,
         "start_time": 0.0, "end_time": 335.16000000002384,
         "start_checkpoint_sha256": None,
         "end_checkpoint_sha256": "d190c50b68328629e3c28779df4e5dcc8f3f80670a54a507d74ee34ee3513665",
         "predecessor_run_id": None, "successor_run_id": None,
         "terminal_state": "COMPLETED", "scientific_classification": "SEPARATE_MODEL_BASE_TIMESTEP",
         "retrieval_state": "VERIFIED", "verified_result_archive_sha256": sha256(fft_archive)},
        {**common, "accepted_timestep": 0.02, "model_identity": "FFT_EIGENSTRAIN_V2",
         "logical_trajectory_id": "FFT_EIGENSTRAIN_V2_seed5101_dt002", "segment_index": 0,
         "run_id": "20260911T110924Z-nogit-6b926f", "slurm_job_id": "55950433",
         "source_commit": "a173624e4a0b3b3e5074166a2fe0229ed2f34398",
         "source_archive_sha256": "c0ba9161a8f519436deec75c7017f03614fa04de43328221996b9efd0dd545b3",
         "resolved_configuration_sha256": None, "mechanics_backend": "fft_eigenstrain_v2",
         "start_step": 0, "end_step": None, "start_time": 0.0, "end_time": None,
         "start_checkpoint_sha256": None, "end_checkpoint_sha256": None,
         "predecessor_run_id": None, "successor_run_id": None,
         "terminal_state": "RUNNING", "scientific_classification": "SEPARATE_TIMESTEP_REFINEMENT",
         "retrieval_state": "NOT_FETCHABLE_WHILE_RUNNING"},
        {"model_identity": "QIU_SI_REFERENCE", "logical_trajectory_id": None,
         "segment_index": None, "run_id": "20260912T180419Z-nogit-a25f0c",
         "slurm_job_id": None, "source_commit": None,
         "source_archive_sha256": "2783b5c019323ee0fb04742e35226c18938b5ddc7b704a3990ecc6025097dbcd",
         "resolved_configuration_sha256": None,
         "initial_field_sha256": "2740e0af26bd3acd24bfbbca44b8b3cd56e2e3bacf97ef5648ab69b84989bf90",
         "orientation_reference_sha256": None, "mechanics_backend": "native_qiu_si_line_disconnection",
         "grid": [500, 500], "domain": [500.0, 500.0], "accepted_timestep": 0.1,
         "start_step": 0, "end_step": None, "start_time": 0.0, "end_time": None,
         "start_checkpoint_sha256": None, "end_checkpoint_sha256": None,
         "predecessor_run_id": None, "successor_run_id": None,
         "terminal_state": "PREPARED_UNSUBMITTED", "scientific_classification": "AUTHORITATIVE_NATIVE_BASELINE_PENDING",
         "retrieval_state": "NOT_SUBMITTED"},
    ]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    repo = Path(__file__).resolve().parents[1]
    rows = records(repo)
    audit = {
        "schema_version": 1,
        "classification": "QIU_ISOTROPIC_LINEAGE_INCOMPLETE",
        "family_classifications": {
            "QIU_LEGACY_FORENSIC": "QIU_ISOTROPIC_LINEAGE_DUPLICATED",
            "FFT_EIGENSTRAIN_V2": "QIU_ISOTROPIC_LINEAGE_INCOMPATIBLE",
            "QIU_SI_REFERENCE": "QIU_ISOTROPIC_LINEAGE_INCOMPLETE",
        },
        "authoritative_isotropic_baseline": "QIU_SI_REFERENCE",
        "exact_continuation_chain_exists": False,
        "reason": "The legacy diagnostic replay overlaps steps 9000-9984; FFT runs use distinct timesteps; the native SI plan is unsubmitted.",
        "source_references": REFS,
        "records": rows,
    }
    (args.output / "qiu_lineage_audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    fields = sorted({key for row in rows for key in row})
    with (args.output / "qiu_lineage_audit.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: json.dumps(v) if isinstance(v, (list, dict)) else v for k, v in row.items()})

    comparisons = []
    ordered = list(REFS.items())
    for (left_name, left), (right_name, right) in zip(ordered, ordered[1:]):
        for row in changed_files(left, right):
            row["left_name"], row["right_name"] = left_name, right_name
            comparisons.append(row)
    with (args.output / "source_file_comparison.csv").open("w", newline="") as handle:
        fields = ["left_name", "right_name", "left", "right", "status", "path", "classification"]
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader(); writer.writerows(comparisons)

    graph = """# Qiu lineage graph

```text
QIU_LEGACY_FORENSIC seed5101
  55930486: step 0 --------------------------> 9984
                    \\ 55948258 replay: 9000 -> 9984
  overlap is an exact diagnostic replay, not a continuation edge

FFT_EIGENSTRAIN_V2 seed5101
  55932457: dt=0.04, step 0 -> 8379, terminal
  55950433: dt=0.02, step 0 -> running, separate refinement

QIU_SI_REFERENCE [100] 4ref
  20260912T180419Z-nogit-a25f0c: prepared, unsubmitted
```

No pair of jobs satisfies all ten continuation requirements. The authoritative
isotropic baseline is the pristine native `QIU_SI_REFERENCE`; it has no result
yet. `QIU_LEGACY_FORENSIC` and `FFT_EIGENSTRAIN_V2` remain separate identities.
"""
    (args.output / "qiu_lineage_graph.md").write_text(graph)
    svg = """<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="520" viewBox="0 0 1200 520">
<style>text{font:18px sans-serif}.h{font-weight:bold;font-size:22px}.n{fill:#eef3f8;stroke:#24445f;stroke-width:2}.a{stroke:#24445f;stroke-width:3;marker-end:url(#m)}.r{stroke:#a34835;stroke-width:3;stroke-dasharray:8 6;marker-end:url(#r)}</style>
<defs><marker id="m" markerWidth="8" markerHeight="8" refX="7" refY="4" orient="auto"><path d="M0,0 L8,4 L0,8z" fill="#24445f"/></marker><marker id="r" markerWidth="8" markerHeight="8" refX="7" refY="4" orient="auto"><path d="M0,0 L8,4 L0,8z" fill="#a34835"/></marker></defs>
<text class="h" x="30" y="35">Immutable Qiu and FFT lineage audit</text>
<text x="30" y="85">QIU_LEGACY_FORENSIC</text><rect class="n" x="280" y="55" width="250" height="55" rx="8"/><text x="300" y="88">55930486: 0 → 9984</text>
<rect class="n" x="700" y="55" width="330" height="55" rx="8"/><text x="720" y="88">55948258: 9000 → 9984 replay</text><line class="r" x1="530" y1="82" x2="695" y2="82"/><text x="550" y="130" fill="#a34835">overlap; no continuation edge</text>
<text x="30" y="220">FFT_EIGENSTRAIN_V2</text><rect class="n" x="280" y="180" width="300" height="55" rx="8"/><text x="300" y="213">55932457: dt .04, terminal</text><rect class="n" x="700" y="180" width="320" height="55" rx="8"/><text x="720" y="213">55950433: dt .02, running</text><text x="470" y="265">separate refinement calculations</text>
<text x="30" y="350">QIU_SI_REFERENCE</text><rect class="n" x="280" y="310" width="560" height="70" rx="8"/><text x="300" y="340">20260912T180419Z-nogit-a25f0c</text><text x="300" y="367">native [100] 4ref: prepared, unsubmitted</text>
<text class="h" x="30" y="465">Classification: QIU_ISOTROPIC_LINEAGE_INCOMPLETE</text><text x="30" y="495">Authoritative baseline: pristine native QIU_SI_REFERENCE</text></svg>"""
    (args.output / "qiu_lineage_graph.svg").write_text(svg)


if __name__ == "__main__":
    main()
