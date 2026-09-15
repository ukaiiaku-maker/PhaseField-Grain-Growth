#!/usr/bin/env python3
"""Build the compact post-simulation analysis for the Phase-1 postfix run."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any, Iterable

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import hsv_to_rgb

from grain_growth_pf.pf.geometry import voronoi_polycrystal
from scripts.run_anisotropic_reduced_pilot import morphology


SOURCE_COMMIT = "1316cc89dbabfb41cb883b0d4a4c74738cc2bef6"
EXPECTED_INITIAL_SHA256 = "2986f2bf744107c46aeedf05d00849a5f35db8c3ba85f70f7a7645b7922b2e61"
EXPECTED_ORIENTATIONS_SHA256 = "217578bc628ee08179cf59123c21a9af7b31ecdfc4eb8ee6a0816446221f3141"
RESULT_ARCHIVE_SHA256 = "b1d3eec6b0b693641958c78fcc73bf2bd7a83bc827eff861fdd661cfd116fc42"
SHAPE = (192, 192)
GRAINS = 200
SEED = 5101
ENERGY_THRESHOLD = 0.01
MORPHOLOGY_THRESHOLD = 0.10
MORPHOLOGY_KEYS = ("boundary_density", "grain_area_cv", "area_weighted_mean_radius")
AMPLITUDE_THRESHOLDS = (1e-14, 1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2)
ORIENTATION_BINS = 36
COLORS = {"A0": "#59636f", "A2 coarse": "#d16b37", "A2 half dt": "#2e778d"}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_array(values: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(values).view(np.uint8)).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def load_eta(path: Path) -> np.ndarray:
    with np.load(path, allow_pickle=False) as archive:
        return np.asarray(archive["eta"])


def write_csv(path: Path, columns: Iterable[str], rows: Iterable[Iterable[Any]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(columns)
        writer.writerows(rows)


def save_figure(fig: mpl.figure.Figure, output: Path, stem: str) -> list[Path]:
    png = output / f"{stem}.png"
    pdf = output / f"{stem}.pdf"
    fig.savefig(png, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return [png, pdf]


def energy_metrics(energies: np.ndarray) -> dict[str, Any]:
    increments = np.diff(energies)
    tolerance = max(abs(float(energies[0])) * 1e-10, 1e-10)
    return {
        "states": int(energies.size),
        "initial": float(energies[0]),
        "final": float(energies[-1]),
        "absolute_change": float(energies[-1] - energies[0]),
        "relative_change": float((energies[-1] - energies[0]) / abs(energies[0])),
        "maximum_increment": float(np.max(increments)),
        "minimum_increment": float(np.min(increments)),
        "positive_steps_above_tolerance": int(np.count_nonzero(increments > tolerance)),
        "roundoff_tolerance": tolerance,
    }


def field_metrics(reference: np.ndarray, candidate: np.ndarray) -> dict[str, float]:
    delta = candidate - reference
    labels_reference = np.argmax(reference, axis=0)
    labels_candidate = np.argmax(candidate, axis=0)
    return {
        "phase_field_rms": float(np.sqrt(np.mean(np.sum(delta * delta, axis=0)))),
        "phase_field_mean_l1": float(np.mean(np.sum(np.abs(delta), axis=0))),
        "dominant_label_disagreement_fraction": float(np.mean(labels_reference != labels_candidate)),
        "maximum_absolute_phase_difference": float(np.max(np.abs(delta))),
    }


def amplitude_support_metrics(eta: np.ndarray) -> dict[str, Any]:
    inverse_participation = 1.0 / np.sum(eta * eta, axis=0)
    threshold_rows = []
    for threshold in AMPLITUDE_THRESHOLDS:
        counts = np.count_nonzero(eta > threshold, axis=0)
        threshold_rows.append({
            "threshold": threshold,
            "mean_phases_per_cell": float(np.mean(counts)),
            "p95_phases_per_cell": float(np.percentile(counts, 95)),
            "maximum_phases_per_cell": int(np.max(counts)),
        })
    return {
        "inverse_participation_mean": float(np.mean(inverse_participation)),
        "inverse_participation_p95": float(np.percentile(inverse_participation, 95)),
        "inverse_participation_maximum": float(np.max(inverse_participation)),
        "threshold_counts": threshold_rows,
    }


def crystal_frame_interface_distribution(eta: np.ndarray, orientations: np.ndarray) -> dict[str, Any]:
    grad_x = 0.5 * (np.roll(eta, -1, axis=2) - np.roll(eta, 1, axis=2))
    grad_y = 0.5 * (np.roll(eta, -1, axis=1) - np.roll(eta, 1, axis=1))
    weights = np.sqrt(grad_x * grad_x + grad_y * grad_y)
    relative_angle = np.mod(np.arctan2(grad_y, grad_x) - orientations[:, None, None], np.pi)
    edges = np.linspace(0.0, np.pi, ORIENTATION_BINS + 1)
    histogram, _ = np.histogram(relative_angle, bins=edges, weights=weights)
    density = histogram / (np.sum(histogram) * np.diff(edges))
    total = float(np.sum(weights))
    harmonics = {
        str(order): float(abs(np.sum(weights * np.exp(1j * order * relative_angle))) / total)
        for order in (2, 4, 6, 8)
    }
    return {
        "angle_centers_degrees": list(map(float, np.degrees(0.5 * (edges[:-1] + edges[1:])))),
        "density_per_radian": list(map(float, density)),
        "total_gradient_weight": total,
        "harmonic_magnitudes": harmonics,
    }


def orientation_image(eta: np.ndarray, orientations: np.ndarray) -> np.ndarray:
    labels = np.argmax(eta, axis=0)
    hue = np.mod(orientations[labels], np.pi) / np.pi
    hsv = np.stack((hue, np.full_like(hue, 0.72), np.full_like(hue, 0.90)), axis=-1)
    rgb = hsv_to_rgb(hsv)
    boundary = (labels != np.roll(labels, 1, 0)) | (labels != np.roll(labels, 1, 1))
    rgb[boundary] = 0.08
    return rgb


def configure_style() -> None:
    mpl.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "legend.fontsize": 8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.dpi": 120,
        "savefig.dpi": 300,
    })


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path, help="extracted output/postfix directory")
    parser.add_argument("output", type=Path, help="compact analysis output directory")
    args = parser.parse_args()
    source = args.input.resolve()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    configure_style()

    a0 = load_json(source / "A0/case-summary.json")
    coarse = load_json(source / "A2_energy_only/case-summary.json")
    fine_checkpoint = load_json(source / "A2_energy_only_half_dt/checkpoint.json")
    if a0["source_commit"] != SOURCE_COMMIT or coarse["source_commit"] != SOURCE_COMMIT:
        raise ValueError("case summary scientific source mismatch")
    if fine_checkpoint["source_commit"] != SOURCE_COMMIT:
        raise ValueError("fine checkpoint scientific source mismatch")

    eta0, _, orientations = voronoi_polycrystal(SHAPE, GRAINS, SEED, width=2.0, periodic=True)
    if sha256_array(eta0) != EXPECTED_INITIAL_SHA256:
        raise ValueError("reconstructed initial state hash mismatch")
    if sha256_array(orientations) != EXPECTED_ORIENTATIONS_SHA256:
        raise ValueError("reconstructed orientation hash mismatch")
    eta_a0 = load_eta(source / "A0/continuous-final.npz")
    eta_coarse = load_eta(source / "A2_energy_only/continuous-final.npz")
    eta_fine = load_eta(source / "A2_energy_only_half_dt/continuous-final.npz")
    eta_restart_partial = load_eta(source / "A2_energy_only_half_dt/checkpoint.npz")

    energies = {
        "A0": np.asarray([row for row in np.loadtxt(source / "A0/energy.csv", delimiter=",", skiprows=1)[:, 1]]),
        "A2 coarse": np.asarray([row for row in np.loadtxt(source / "A2_energy_only/energy.csv", delimiter=",", skiprows=1)[:, 1]]),
        "A2 half dt": np.asarray(fine_checkpoint["energies"], dtype=float),
    }
    dts = {
        "A0": float(a0["accepted_dt"]),
        "A2 coarse": float(coarse["accepted_dt"]),
        "A2 half dt": float(coarse["accepted_dt"]) / 2.0,
    }
    energy_summary = {name: energy_metrics(values) for name, values in energies.items()}

    morphologies = {
        "initial": morphology(eta0),
        "A0": morphology(eta_a0),
        "A2 coarse": morphology(eta_coarse),
        "A2 half dt": morphology(eta_fine),
    }
    convergence = {
        "final_energy_relative_difference": abs(energy_summary["A2 coarse"]["final"] - energy_summary["A2 half dt"]["final"])
        / abs(energy_summary["A2 half dt"]["final"]),
        "energy_threshold": ENERGY_THRESHOLD,
        "morphology_relative_differences": {
            key: abs(float(morphologies["A2 coarse"][key]) - float(morphologies["A2 half dt"][key]))
            / abs(float(morphologies["A2 half dt"][key]))
            for key in MORPHOLOGY_KEYS
        },
        "morphology_threshold": MORPHOLOGY_THRESHOLD,
    }
    convergence["passed"] = bool(
        convergence["final_energy_relative_difference"] <= ENERGY_THRESHOLD
        and max(convergence["morphology_relative_differences"].values()) <= MORPHOLOGY_THRESHOLD
    )

    fields = {
        "A0_vs_initial": field_metrics(eta0, eta_a0),
        "A2_coarse_vs_initial": field_metrics(eta0, eta_coarse),
        "A2_half_dt_vs_initial": field_metrics(eta0, eta_fine),
        "A2_coarse_vs_half_dt": field_metrics(eta_fine, eta_coarse),
        "fine_partial_restart_vs_continuous_final": field_metrics(eta_fine, eta_restart_partial),
    }
    state_fields = {"initial": eta0, "A0": eta_a0, "A2 coarse": eta_coarse, "A2 half dt": eta_fine}
    amplitude_support = {name: amplitude_support_metrics(eta) for name, eta in state_fields.items()}
    interface_orientation = {
        name: crystal_frame_interface_distribution(eta, orientations)
        for name, eta in state_fields.items()
    }
    initial_areas = eta0.sum(axis=(1, 2))
    grain_area_response = {}
    for name, eta in state_fields.items():
        relative = (eta.sum(axis=(1, 2)) - initial_areas) / initial_areas
        grain_area_response[name] = {
            "mean_relative_change": float(np.mean(relative)),
            "rms_relative_change": float(np.sqrt(np.mean(relative * relative))),
            "maximum_absolute_relative_change": float(np.max(np.abs(relative))),
        }

    support_sources = {
        "A0": a0["support_history"],
        "A2 coarse": coarse["support_history"],
        "A2 half dt": fine_checkpoint["support_history"],
    }
    support_summary: dict[str, dict[str, Any]] = {}
    for name, history in support_sources.items():
        support_summary[name] = {
            "states": len(history),
            "maximum_active_phases_per_cell": max(row["maximum_active_phases"] for row in history),
            "maximum_active_pairs_per_cell": max(row["maximum_active_pairs"] for row in history),
            "maximum_supported_pair_instances": max(row["audit_supported_pair_instances"] for row in history),
            "maximum_zero_gradient_supported_pair_instances": max(row["zero_gradient_supported_pair_instances"] for row in history),
            "minimum_nonzero_pair_gradient": min(row["minimum_nonzero_pair_gradient"] for row in history),
            "maximum_phase_sum_error": max(row["maximum_phase_sum_error"] for row in history),
            "minimum_phase_value": min(row["minimum_phase_value"] for row in history),
            "all_finite": all(row["finite"] for row in history),
        }

    restart = {
        "A0": {
            "complete": True,
            "exact": bool(a0["validity"]["restart_exact"] and a0["validity"]["restart_time_exact"]),
            "continuous_field_sha256": a0["final_field_sha256"],
            "restart_field_sha256": a0["restart_field_sha256"],
        },
        "A2 coarse": {
            "complete": True,
            "exact": bool(coarse["validity"]["restart_exact"] and coarse["validity"]["restart_time_exact"]),
            "continuous_field_sha256": coarse["final_field_sha256"],
            "restart_field_sha256": coarse["restart_field_sha256"],
        },
        "A2 half dt": {
            "complete": False,
            "exact": None,
            "accepted_step": int(fine_checkpoint["accepted_step"]),
            "target_step": 128,
            "continuous_field_sha256": sha256_array(eta_fine),
            "partial_restart_field_sha256": sha256_array(eta_restart_partial),
        },
    }

    timings = {
        "A0": a0["timing_history"],
        "A2 coarse": coarse["timing_history"],
        "A2 half dt": fine_checkpoint["timing_history"],
    }
    summary = {
        "schema": "anisotropic-postfix-postanalysis-v1",
        "analysis_basis": "completed continuous trajectories plus exact fine restart through step 126 of 128",
        "classification": "A2_POSTFIX_OPERATIONALLY_INCOMPLETE",
        "scientific_source_commit": SOURCE_COMMIT,
        "result_archive_sha256": RESULT_ARCHIVE_SHA256,
        "initial_state_sha256": sha256_array(eta0),
        "orientations_sha256": sha256_array(orientations),
        "energy": energy_summary,
        "morphology": morphologies,
        "timestep_convergence": convergence,
        "field_comparisons": fields,
        "amplitude_support": amplitude_support,
        "crystal_frame_interface_orientation": interface_orientation,
        "grain_area_response": grain_area_response,
        "support": support_summary,
        "restart": restart,
        "release": {
            "mobility_only": False,
            "combined_A2": False,
            "A3": False,
            "production": False,
        },
        "limitations": [
            "fine midpoint restart is missing accepted steps 127 and 128",
            "the physical horizon is a qualification horizon, not a grain-growth kinetics campaign",
            "no mobility-only or combined anisotropy trajectory exists in this archive",
            "single deterministic initial condition; no ensemble uncertainty estimate is available",
        ],
    }

    energy_rows = []
    increment_rows = []
    for name, values in energies.items():
        times = np.arange(values.size) * dts[name]
        for index, (time_value, energy) in enumerate(zip(times, values)):
            energy_rows.append((name, index, time_value, energy, (energy - values[0]) / values[0]))
        for index, increment in enumerate(np.diff(values), start=1):
            increment_rows.append((name, index, times[index], increment, increment / values[0]))
    write_csv(output / "energy_history.csv", ("case", "state_index", "physical_time", "energy", "relative_change"), energy_rows)
    write_csv(output / "energy_increments.csv", ("case", "accepted_step", "physical_time", "energy_increment", "relative_increment"), increment_rows)
    write_csv(
        output / "morphology_summary.csv",
        ("case", "active_grains", "boundary_density", "grain_area_cv", "area_weighted_mean_radius", "isolated_label_fraction"),
        ((name, row["active_grains"], row["boundary_density"], row["grain_area_cv"], row["area_weighted_mean_radius"], row["isolated_label_fraction"])
         for name, row in morphologies.items()),
    )
    support_rows = []
    for name, history in support_sources.items():
        for row in history:
            support_rows.append((name, row["state_index"], row["physical_time"], row["maximum_active_phases"],
                                 row["maximum_active_pairs"], row["audit_supported_pair_instances"],
                                 row["zero_gradient_supported_pair_instances"], row["minimum_nonzero_pair_gradient"],
                                 row["maximum_phase_sum_error"], row["minimum_phase_value"], row["finite"]))
    write_csv(output / "support_history.csv", ("case", "state_index", "physical_time", "maximum_active_phases",
              "maximum_active_pairs", "supported_pair_instances", "zero_gradient_pair_instances",
              "minimum_nonzero_pair_gradient", "maximum_phase_sum_error", "minimum_phase_value", "finite"), support_rows)
    timing_rows = []
    for name, history in timings.items():
        previous: dict[str, int] = {}
        for row in history:
            stage = row["stage"]
            start = previous.get(stage, 0 if stage == "continuous" else (32 if name != "A2 half dt" else 64))
            count = int(row["accepted_step"]) - start
            timing_rows.append((name, stage, start, row["accepted_step"], count, row["block_seconds"], row["block_seconds"] / count))
            previous[stage] = int(row["accepted_step"])
    write_csv(output / "timing_blocks.csv", ("case", "stage", "start_step", "end_step", "accepted_steps", "block_seconds", "seconds_per_step"), timing_rows)
    amplitude_rows = []
    for name, metrics in amplitude_support.items():
        for row in metrics["threshold_counts"]:
            amplitude_rows.append((name, row["threshold"], row["mean_phases_per_cell"], row["p95_phases_per_cell"], row["maximum_phases_per_cell"],
                                   metrics["inverse_participation_mean"], metrics["inverse_participation_p95"], metrics["inverse_participation_maximum"]))
    write_csv(output / "amplitude_support.csv", ("case", "amplitude_threshold", "mean_phases_per_cell", "p95_phases_per_cell",
              "maximum_phases_per_cell", "inverse_participation_mean", "inverse_participation_p95", "inverse_participation_maximum"), amplitude_rows)
    orientation_rows = []
    for name, metrics in interface_orientation.items():
        for angle, density in zip(metrics["angle_centers_degrees"], metrics["density_per_radian"]):
            orientation_rows.append((name, angle, density, metrics["total_gradient_weight"],
                                     metrics["harmonic_magnitudes"]["2"], metrics["harmonic_magnitudes"]["4"],
                                     metrics["harmonic_magnitudes"]["6"], metrics["harmonic_magnitudes"]["8"]))
    write_csv(output / "interface_orientation.csv", ("case", "crystal_frame_normal_degrees", "gradient_weighted_density_per_radian",
              "total_gradient_weight", "harmonic_2", "harmonic_4", "harmonic_6", "harmonic_8"), orientation_rows)

    # Figure 1: exact implemented energy and increments.
    fig, axes = plt.subplots(2, 1, figsize=(7.2, 6.0), sharex=False)
    for name, values in energies.items():
        time_values = np.arange(values.size) * dts[name]
        axes[0].plot(time_values, values, lw=1.7, color=COLORS[name], label=name)
        axes[1].plot(time_values[1:], np.diff(values), lw=1.25, color=COLORS[name], label=name)
    axes[0].set(title="Exact implemented unforced energy", xlabel="Physical time", ylabel="Energy")
    axes[0].legend(frameon=False, ncol=3)
    axes[1].axhline(0, color="black", lw=0.8)
    axes[1].set(title="Accepted-step energy increments", xlabel="Physical time", ylabel=r"$E_{n+1}-E_n$")
    axes[1].legend(frameon=False, ncol=3)
    fig.tight_layout()
    generated = save_figure(fig, output, "energy_evolution")

    # Figure 2: convergence against preregistered thresholds and morphology change.
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.6))
    conv_names = ["Energy", "Boundary\ndensity", "Area CV", "Weighted\nradius"]
    conv_values = [convergence["final_energy_relative_difference"], *(convergence["morphology_relative_differences"][key] for key in MORPHOLOGY_KEYS)]
    conv_limits = [ENERGY_THRESHOLD, MORPHOLOGY_THRESHOLD, MORPHOLOGY_THRESHOLD, MORPHOLOGY_THRESHOLD]
    normalized = np.asarray(conv_values) / np.asarray(conv_limits)
    axes[0].bar(np.arange(4), normalized, color=[COLORS["A2 half dt"]] * 4)
    axes[0].axhline(1, color="#a12424", ls="--", lw=1.2, label="Preregistered limit")
    axes[0].set(xticks=np.arange(4), xticklabels=conv_names, ylabel="Observed / allowed relative difference", title="Coarse–half-dt convergence")
    axes[0].legend(frameon=False)
    x = np.arange(len(MORPHOLOGY_KEYS))
    width = 0.25
    for offset, name in zip((-width, 0, width), ("A0", "A2 coarse", "A2 half dt")):
        changes = [(morphologies[name][key] - morphologies["initial"][key]) / abs(morphologies["initial"][key]) for key in MORPHOLOGY_KEYS]
        axes[1].bar(x + offset, 100 * np.asarray(changes), width, label=name, color=COLORS[name])
    axes[1].axhline(0, color="black", lw=0.8)
    axes[1].set(xticks=x, xticklabels=["Boundary\ndensity", "Area CV", "Weighted\nradius"], ylabel="Change from initial (%)", title="Morphology response")
    axes[1].legend(frameon=False)
    fig.tight_layout()
    generated += save_figure(fig, output, "convergence_morphology")

    # Figure 3: orientation-colored microstructures and direct field differences.
    fig, axes = plt.subplots(2, 3, figsize=(9.2, 6.1))
    panels = ((eta0, "Initial"), (eta_a0, "A0 final"), (eta_coarse, "A2 coarse final"), (eta_fine, "A2 half-dt final"))
    for ax, (eta, title) in zip(axes.flat[:4], panels):
        ax.imshow(orientation_image(eta, orientations), interpolation="nearest")
        ax.set_title(title)
        ax.set_axis_off()
    coarse_fine_difference = np.sqrt(np.sum((eta_coarse - eta_fine) ** 2, axis=0))
    initial_fine_difference = np.sqrt(np.sum((eta_fine - eta0) ** 2, axis=0))
    images = ((coarse_fine_difference, "Coarse–half-dt field distance"), (initial_fine_difference, "Half-dt evolution distance"))
    for ax, (data, title) in zip(axes.flat[4:], images):
        image_artist = ax.imshow(data, cmap="magma", interpolation="nearest")
        ax.set_title(title)
        ax.set_axis_off()
        fig.colorbar(image_artist, ax=ax, fraction=0.046, pad=0.03)
    fig.tight_layout()
    generated += save_figure(fig, output, "microstructure_fields")

    # Figure 4: grain-area distributions and grain-resolved response.
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.5))
    for name, eta in (("initial", eta0), ("A0", eta_a0), ("A2 coarse", eta_coarse), ("A2 half dt", eta_fine)):
        areas = eta.sum(axis=(1, 2))
        scaled = np.sort(areas / np.mean(areas))
        color = "black" if name == "initial" else COLORS[name]
        axes[0].plot(scaled, np.linspace(0, 1, scaled.size, endpoint=False), color=color, lw=1.5, label=name)
        if name != "initial":
            relative_area_change = 100 * (areas - initial_areas) / initial_areas
            axes[1].scatter(initial_areas / np.mean(initial_areas), relative_area_change, s=8, alpha=0.55, color=color, label=name)
    axes[0].set(xlabel="Grain area / mean area", ylabel="Empirical cumulative probability", title="Grain-area CDF")
    axes[1].axhline(0, color="black", lw=0.8)
    axes[1].set(xlabel="Initial grain area / mean area", ylabel="Per-grain area change (%)", title="Grain-resolved response")
    axes[1].legend(frameon=False)
    fig.tight_layout()
    generated += save_figure(fig, output, "grain_area_statistics")

    # Figure 5: support growth, the dominant operational finding.
    fig, axes = plt.subplots(2, 2, figsize=(8.8, 6.5), sharex=False)
    fields_to_plot = (
        ("maximum_active_phases", "Maximum active phases per cell", False),
        ("maximum_active_pairs", "Maximum active pairs per cell", True),
        ("audit_supported_pair_instances", "Supported pair instances", True),
        ("zero_gradient_supported_pair_instances", "Zero-gradient supported pairs", True),
    )
    for ax, (key, title, logarithmic) in zip(axes.flat, fields_to_plot):
        for name, history in support_sources.items():
            ax.plot([row["state_index"] for row in history], [row[key] for row in history], lw=1.5, color=COLORS[name], label=name)
        if logarithmic:
            ax.set_yscale("symlog", linthresh=1)
        ax.set(title=title, xlabel="Continuous accepted step")
    axes[0, 0].legend(frameon=False)
    fig.tight_layout()
    generated += save_figure(fig, output, "support_growth")

    # Figure 6: distinguish tail support from substantial phase mixing.
    fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.6))
    for name in ("initial", "A0", "A2 coarse", "A2 half dt"):
        rows = amplitude_support[name]["threshold_counts"]
        color = "black" if name == "initial" else COLORS[name]
        axes[0].semilogx([row["threshold"] for row in rows], [row["mean_phases_per_cell"] for row in rows], marker="o", ms=3, color=color, label=name)
        axes[1].semilogx([row["threshold"] for row in rows], [row["p95_phases_per_cell"] for row in rows], marker="o", ms=3, color=color, label=name)
    for ax, title, ylabel in ((axes[0], "Mean phase support", "Mean phases per cell"), (axes[1], "95th-percentile phase support", "P95 phases per cell")):
        ax.invert_xaxis()
        ax.set(title=title, xlabel="Phase-amplitude threshold", ylabel=ylabel)
    axes[0].legend(frameon=False)
    fig.tight_layout()
    generated += save_figure(fig, output, "amplitude_support")

    # Figure 7: interface normals in the corresponding grain's crystal frame.
    fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.6))
    initial_density = np.asarray(interface_orientation["initial"]["density_per_radian"])
    for name in ("initial", "A0", "A2 coarse", "A2 half dt"):
        metrics = interface_orientation[name]
        color = "black" if name == "initial" else COLORS[name]
        angle = np.asarray(metrics["angle_centers_degrees"])
        density = np.asarray(metrics["density_per_radian"])
        axes[0].plot(angle, density, color=color, lw=1.5, label=name)
        if name != "initial":
            axes[1].plot(angle, density - initial_density, color=color, lw=1.5, label=name)
    axes[0].axhline(1 / np.pi, color="#888888", ls=":", lw=1, label="uniform")
    axes[0].set(xlabel="Interface normal in crystal frame (degrees)", ylabel="Gradient-weighted density (rad$^{-1}$)", title="Crystal-frame interface normals")
    axes[1].axhline(0, color="black", lw=0.8)
    axes[1].set(xlabel="Interface normal in crystal frame (degrees)", ylabel="Density change from initial", title="Directional evolution response")
    axes[0].legend(frameon=False, fontsize=7)
    axes[1].legend(frameon=False, fontsize=7)
    fig.tight_layout()
    generated += save_figure(fig, output, "interface_orientation")

    # Figure 8: measured block cost and its relation to support size.
    fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.6))
    for name in ("A0", "A2 coarse", "A2 half dt"):
        for stage, linestyle in (("continuous", "-"), ("restart", "--")):
            rows = [row for row in timing_rows if row[0] == name and row[1] == stage]
            if rows:
                axes[0].plot([row[3] for row in rows], [row[6] for row in rows], marker="o", ms=3,
                             color=COLORS[name], ls=linestyle, label=f"{name} {stage}")
    axes[0].set_yscale("log")
    axes[0].set(xlabel="Accepted step", ylabel="Seconds per accepted step", title="Measured block cost")
    axes[0].legend(frameon=False, fontsize=7)
    fine_continuous = [row for row in timing_rows if row[0] == "A2 half dt" and row[1] == "continuous"]
    fine_support = {row["state_index"]: row["audit_supported_pair_instances"] for row in support_sources["A2 half dt"]}
    axes[1].loglog([fine_support[row[3]] for row in fine_continuous], [row[6] for row in fine_continuous], "o-", color=COLORS["A2 half dt"])
    axes[1].set(xlabel="Supported pair instances at block end", ylabel="Seconds per accepted step", title="Cost follows support proliferation")
    fig.tight_layout()
    generated += save_figure(fig, output, "runtime_scaling")

    summary_path = output / "analysis_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    (output / "README.md").write_text(
        "# Anisotropic postfix compact analysis\n\n"
        f"Scientific source: `{SOURCE_COMMIT}`. Verified result archive: `{RESULT_ARCHIVE_SHA256}`.\n\n"
        "The three continuous trajectories are complete. A0 and coarse A2 restart exactly; the half-dt restart ends at step 126/128. "
        "The formal classification remains `A2_POSTFIX_OPERATIONALLY_INCOMPLETE`.\n\n"
        f"Coarse A2 energy decreases from {energy_summary['A2 coarse']['initial']:.12f} to {energy_summary['A2 coarse']['final']:.12f}; "
        f"half-dt A2 decreases to {energy_summary['A2 half dt']['final']:.12f}. Both have zero positive accepted-step increments. "
        f"The matched-horizon energy difference is {100 * convergence['final_energy_relative_difference']:.4f}%.\n\n"
        "The main operational result is low-amplitude support proliferation: the half-dt final field has mean support 116.8 phases per cell at 1e-14, "
        "but 1.62 at amplitude 0.01. This drives the late-stage cost increase.\n\n"
        "Figures: [energy](energy_evolution.png), [convergence and morphology](convergence_morphology.png), "
        "[microstructures](microstructure_fields.png), [grain areas](grain_area_statistics.png), "
        "[support growth](support_growth.png), [amplitude support](amplitude_support.png), "
        "[interface orientation](interface_orientation.png), and [runtime scaling](runtime_scaling.png).\n\n"
        "Machine-readable results are in `analysis_summary.json`; plotted data are in the CSV files; hashes are in `artifact_manifest.json`.\n"
    )
    provenance = {
        "schema": "anisotropic-postfix-analysis-provenance-v1",
        "input_root": str(source),
        "scientific_source_commit": SOURCE_COMMIT,
        "result_archive_sha256": RESULT_ARCHIVE_SHA256,
        "analysis_script": "scripts/analyze_anisotropic_postfix.py",
        "analysis_script_sha256": sha256_file(Path(__file__).resolve()),
        "case_summary_sha256": {
            "A0": sha256_file(source / "A0/case-summary.json"),
            "A2_energy_only": sha256_file(source / "A2_energy_only/case-summary.json"),
        },
        "fine_checkpoint_sha256": sha256_file(source / "A2_energy_only_half_dt/checkpoint.json"),
        "fine_continuous_state_sha256": sha256_file(source / "A2_energy_only_half_dt/continuous-final.npz"),
    }
    (output / "provenance.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    artifact_paths = sorted([*generated, *(output / name for name in (
        "README.md", "analysis_summary.json", "provenance.json", "energy_history.csv", "energy_increments.csv",
        "morphology_summary.csv", "support_history.csv", "timing_blocks.csv", "amplitude_support.csv",
        "interface_orientation.csv"))])
    manifest = {
        "schema": "anisotropic-postfix-analysis-manifest-v1",
        "artifacts": [{"path": path.name, "sha256": sha256_file(path), "bytes": path.stat().st_size} for path in artifact_paths],
    }
    (output / "artifact_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
