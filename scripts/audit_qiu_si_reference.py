#!/usr/bin/env python3
"""Audit and stage the pristine Qiu supplementary-information code archive.

This utility never imports the repository's historical ``qiu_full_field``
backend.  It inventories the published ZIP byte-for-byte, identifies missing
native inputs/imports, and can stage one archived benchmark without modifying
the pristine extraction.  Compatibility aliases are explicit, checksummed,
and never promoted to reference evidence by this script.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import zipfile
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath


EXPECTED_ARCHIVE_MD5 = "6cd49ca72eba89210abb96e700342f12"
EXPECTED_ARCHIVE_SHA256 = "2740e0af26bd3acd24bfbbca44b8b3cd56e2e3bacf97ef5648ab69b84989bf90"

BENCHMARKS: dict[str, dict[str, object]] = {
    "four_ref": {
        "driver": "PF_Codes/Idealized_microstructure/[100]_4ref/Bicrystal-4reference-el-pf.py",
        "support": ["PF_Codes/Idealized_microstructure/[100]_4ref/functions_4ref_new.py"],
        "declared_import": "functions_4ref_new",
        "compatibility_alias": None,
        "native_parameters": {"nx": 500, "ny": 500, "number_of_grain": 17,
                              "dx": 1, "dy": 1, "dt": 0.1, "nsteps": 200000,
                              "ref": 4.0, "stress_step": 250, "gb_step": 250},
    },
    "six_ref": {
        "driver": "PF_Codes/Idealized_microstructure/[111]_6ref/Idealized-6reference.py",
        "support": ["PF_Codes/Idealized_microstructure/[111]_6ref/functions_6ref.py"],
        "declared_import": "functions_6ref_new2_1",
        "compatibility_alias": {"functions_6ref_new2_1.py": "functions_6ref.py"},
        "native_parameters": {"nx": 500, "ny": 500, "number_of_grain": 17,
                              "dx": 1, "dy": 1, "dt": 0.1, "nsteps": 200000,
                              "ref": 6, "stress_step": 500, "gb_step": 250},
    },
    "six_ref_continue": {
        "driver": "PF_Codes/Idealized_microstructure/[111]_6ref/Idealized-6reference_continue.py",
        "support": ["PF_Codes/Idealized_microstructure/[111]_6ref/functions_6ref_continue.py"],
        "declared_import": "functions_6ref_new2_2",
        "compatibility_alias": {"functions_6ref_new2_2.py": "functions_6ref_continue.py"},
        "required_checkpoint": "OP_t142500.npz",
        "native_parameters": {"nx": 500, "ny": 500, "number_of_grain": 17,
                              "dx": 1, "dy": 1, "dt": 0.1, "nsteps": 200000,
                              "ref": 6, "stress_step": 500, "gb_step": 250},
    },
    "polycrystal": {
        "driver": "PF_Codes/Polycrystals/Grain_growth_polycrystal.py",
        "support": ["PF_Codes/Polycrystals/functions_2ref.py"],
        "declared_import": "functions_2ref",
        "compatibility_alias": None,
        "required_seed": "PolycrystalSeeds_1000Grains_25_250_1.txt",
        "seed_generator": "PF_Codes/Polycrystals/Initial_configuration_Polycrystal-seeds.py",
        "native_parameters": {"nx": 250, "ny": 250, "number_of_grain": 1001,
                              "dx": 10, "dy": 10, "dt": 10, "nsteps": 300000,
                              "ref": 2.0, "stress_step": 500, "gb_step": 250},
    },
}


def digest_bytes(payload: bytes, algorithm: str = "sha256") -> str:
    digest = hashlib.new(algorithm)
    digest.update(payload)
    return digest.hexdigest()


def digest_file(path: Path, algorithm: str = "sha256") -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def declared_local_imports(source: str) -> list[str]:
    tree = ast.parse(source)
    result: list[str] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and node.module:
            if node.module.startswith("functions_"):
                result.append(node.module)
        elif isinstance(node, ast.Import):
            result.extend(alias.name for alias in node.names if alias.name.startswith("functions_"))
    return sorted(set(result))


def safe_member_name(name: str) -> bool:
    path = PurePosixPath(name)
    return not path.is_absolute() and ".." not in path.parts


def dependency_probe(python: Path) -> dict[str, object]:
    code = (
        "import importlib, json, platform; out={'python':platform.python_version()}; "
        "mods=['numpy','scipy','matplotlib','numba','sparse']; "
        "out.update({m:getattr(importlib.import_module(m),'__version__','unknown') for m in mods}); "
        "print(json.dumps(out, sort_keys=True))"
    )
    completed = subprocess.run([str(python), "-c", code], text=True, capture_output=True)
    return {
        "python_executable": str(python.resolve()),
        "returncode": completed.returncode,
        "versions": json.loads(completed.stdout) if completed.returncode == 0 else {},
        "stderr": completed.stderr.strip(),
    }


def audit_archive(archive: Path, python: Path) -> dict[str, object]:
    archive = archive.resolve()
    sha256 = digest_file(archive)
    md5 = digest_file(archive, "md5")
    with zipfile.ZipFile(archive) as stream:
        names = set(stream.namelist())
        scientific = sorted(
            name for name in names
            if not name.startswith("__MACOSX/") and not name.endswith("/")
        )
        unsafe = sorted(name for name in names if not safe_member_name(name))
        members = {
            name: {"size_bytes": stream.getinfo(name).file_size,
                   "sha256": digest_bytes(stream.read(name))}
            for name in scientific
        }
        benchmark_audits: dict[str, object] = {}
        for label, specification in BENCHMARKS.items():
            driver = str(specification["driver"])
            source = stream.read(driver).decode("utf-8")
            imports = declared_local_imports(source)
            imported_members = [
                str(PurePosixPath(driver).parent / f"{module}.py") for module in imports
            ]
            missing_imports = [member for member in imported_members if member not in names]
            required_seed = specification.get("required_seed")
            seed_member = (
                str(PurePosixPath(driver).parent / str(required_seed))
                if required_seed else None
            )
            required_checkpoint = specification.get("required_checkpoint")
            checkpoint_member = (
                str(PurePosixPath(driver).parent / str(required_checkpoint))
                if required_checkpoint else None
            )
            blockers: list[str] = []
            if missing_imports:
                blockers.append("published_driver_import_target_absent")
            if seed_member and seed_member not in names:
                blockers.append("published_seed_realization_absent")
            if checkpoint_member and checkpoint_member not in names:
                blockers.append("continuation_checkpoint_absent")
            status = "READY_NATIVE_ARCHIVE" if not blockers else "BLOCKED_OR_COMPATIBILITY_REQUIRED"
            benchmark_audits[label] = {
                "status": status,
                "driver": driver,
                "driver_sha256": members[driver]["sha256"],
                "declared_local_imports": imports,
                "missing_import_members": missing_imports,
                "compatibility_alias": specification.get("compatibility_alias"),
                "required_seed_member": seed_member,
                "required_seed_present": seed_member in names if seed_member else None,
                "required_checkpoint_member": checkpoint_member,
                "required_checkpoint_present": checkpoint_member in names if checkpoint_member else None,
                "native_parameters": specification["native_parameters"],
                "blockers": blockers,
                "scientific_baseline_eligible": False,
            }
    return {
        "schema_version": 1,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "workstream": "QIU_SI_REFERENCE",
        "archive": str(archive),
        "archive_sha256": sha256,
        "archive_md5": md5,
        "archive_identity_passed": sha256 == EXPECTED_ARCHIVE_SHA256 and md5 == EXPECTED_ARCHIVE_MD5,
        "expected_archive_sha256": EXPECTED_ARCHIVE_SHA256,
        "expected_archive_md5": EXPECTED_ARCHIVE_MD5,
        "unsafe_archive_members": unsafe,
        "members": members,
        "dependencies": dependency_probe(python),
        "benchmarks": benchmark_audits,
        "terminology": {
            "QIU_SI_REFERENCE": "archived current-geometry line-disconnection implementation",
            "FFT_EIGENSTRAIN_V2": "independent corrected periodic FFT/eigenstrain model",
            "QIU_LEGACY_FORENSIC": "rejected deterministic legacy in-house surrogate",
        },
        "scientific_baseline_eligible": False,
        "promotion_gate": (
            "Requires a completed native benchmark plus hashed intermediate beta, "
            "line/disconnection-density, stress, elastic-force, and phase-field-increment observables."
        ),
    }


def stage_benchmark(archive: Path, benchmark: str, destination: Path) -> dict[str, object]:
    specification = BENCHMARKS[benchmark]
    destination = destination.resolve()
    destination.mkdir(parents=True, exist_ok=False)
    stage = destination / "stage"
    stage.mkdir()
    members = [str(specification["driver"]), *map(str, specification["support"])]
    seed_generator = specification.get("seed_generator")
    if seed_generator:
        members.append(str(seed_generator))
    source_hashes: dict[str, str] = {}
    with zipfile.ZipFile(archive) as stream:
        for member in members:
            payload = stream.read(member)
            target = stage / PurePosixPath(member).name
            target.write_bytes(payload)
            source_hashes[target.name] = digest_bytes(payload)
    aliases: dict[str, object] = {}
    for alias, target_name in (specification.get("compatibility_alias") or {}).items():
        source = stage / target_name
        target = stage / alias
        os.link(source, target)
        aliases[alias] = {
            "target": target_name,
            "sha256": digest_file(target),
            "kind": "explicit_hardlink_compatibility_alias",
        }
    manifest = {
        "schema_version": 1,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "workstream": "QIU_SI_REFERENCE",
        "benchmark": benchmark,
        "archive": str(archive.resolve()),
        "archive_sha256": digest_file(archive),
        "driver": PurePosixPath(str(specification["driver"])).name,
        "source_files_sha256": source_hashes,
        "compatibility_aliases": aliases,
        "native_parameters": specification["native_parameters"],
        "required_seed": specification.get("required_seed"),
        "required_checkpoint": specification.get("required_checkpoint"),
        "execution_state": "STAGED_NOT_STARTED",
        "scientific_baseline_eligible": False,
        "pooling_allowed": False,
    }
    (destination / "stage_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--python", type=Path, default=Path(sys.executable))
    parser.add_argument("--audit-output", type=Path, required=True)
    parser.add_argument("--stage-benchmark", choices=sorted(BENCHMARKS))
    parser.add_argument("--stage-output", type=Path)
    args = parser.parse_args()
    audit = audit_archive(args.archive, args.python)
    args.audit_output.parent.mkdir(parents=True, exist_ok=True)
    args.audit_output.write_text(json.dumps(audit, indent=2) + "\n")
    print(args.audit_output.resolve())
    if args.stage_benchmark:
        if args.stage_output is None:
            parser.error("--stage-output is required with --stage-benchmark")
        manifest = stage_benchmark(args.archive, args.stage_benchmark, args.stage_output)
        print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
