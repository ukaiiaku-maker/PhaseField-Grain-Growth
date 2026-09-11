from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numba import njit
from numpy.typing import NDArray


@njit(cache=True)
def _beta_first(orientation_i: float, orientation_j: float) -> float:
    theta = orientation_i - orientation_j
    if theta >= np.pi:
        theta -= np.pi
    elif theta <= -np.pi:
        theta += np.pi
    magnitude = abs(theta)
    if magnitude <= 70.0 * np.pi / 180.0:
        value = 2.0 * np.tan(theta / 2.0)
    elif magnitude <= 150.0 * np.pi / 180.0:
        value = 2.0 * np.tan(
            theta / 2.0 - np.sign(theta) * 109.38 * np.pi / 360.0
        )
    else:
        value = -2.0 * np.tan(np.sign(theta) * np.pi / 2.0 - theta / 2.0)
    return 0.0333 * value


@njit(cache=True)
def _pair_tensor(
    eta: NDArray[np.float64], i: int, j: int, y: int, x: int,
    orientations: NDArray[np.float64], dx: float,
) -> tuple[float, float, float]:
    ny, nx = eta.shape[1:]
    ym, yp = (y - 1) % ny, (y + 1) % ny
    xm, xp = (x - 1) % nx, (x + 1) % nx
    # n points from canonical phase i toward canonical phase j.
    gy = (
        (eta[j, yp, x] - eta[i, yp, x])
        - (eta[j, ym, x] - eta[i, ym, x])
    ) / (2.0 * dx)
    gx = (
        (eta[j, y, xp] - eta[i, y, xp])
        - (eta[j, y, xm] - eta[i, y, xm])
    ) / (2.0 * dx)
    norm = np.sqrt(gy * gy + gx * gx)
    if norm <= 1e-14:
        return 0.0, 0.0, 0.0
    ny_value, nx_value = gy / norm, gx / norm
    ty_value, tx_value = -nx_value, ny_value
    beta = _beta_first(orientations[i], orientations[j])
    return (
        beta * ty_value * ny_value,
        beta * tx_value * nx_value,
        0.5 * beta * (ty_value * nx_value + ny_value * tx_value),
    )


@njit(cache=True)
def _pair_tensor_midpoint(
    before: NDArray[np.float64], after: NDArray[np.float64],
    i: int, j: int, y: int, x: int,
    orientations: NDArray[np.float64], dx: float,
) -> tuple[float, float, float]:
    ny, nx = before.shape[1:]
    ym, yp = (y - 1) % ny, (y + 1) % ny
    xm, xp = (x - 1) % nx, (x + 1) % nx

    def difference(phase_j: int, phase_i: int, yy: int, xx: int) -> float:
        return 0.5 * (
            before[phase_j, yy, xx] + after[phase_j, yy, xx]
            - before[phase_i, yy, xx] - after[phase_i, yy, xx]
        )

    gy = (difference(j, i, yp, x) - difference(j, i, ym, x)) / (2.0 * dx)
    gx = (difference(j, i, y, xp) - difference(j, i, y, xm)) / (2.0 * dx)
    norm = np.sqrt(gy * gy + gx * gx)
    if norm <= 1e-14:
        return 0.0, 0.0, 0.0
    ny_value, nx_value = gy / norm, gx / norm
    ty_value, tx_value = -nx_value, ny_value
    beta = _beta_first(orientations[i], orientations[j])
    return (
        beta * ty_value * ny_value,
        beta * tx_value * nx_value,
        0.5 * beta * (ty_value * nx_value + ny_value * tx_value),
    )


@njit(cache=True)
def _map_local_transfer(
    before: NDArray[np.float64], after: NDArray[np.float64],
    orientations: NDArray[np.float64], dx: float,
) -> tuple[NDArray[np.float64], float, float]:
    phases, ny, nx = before.shape
    midpoint = 0.5 * (before + after)
    increment = np.zeros((2, 2, ny, nx), dtype=np.float64)
    donors = np.empty(phases, dtype=np.int64)
    receivers = np.empty(phases, dtype=np.int64)
    donor_amount = np.empty(phases, dtype=np.float64)
    receiver_amount = np.empty(phases, dtype=np.float64)
    absolute_sweep = 0.0
    signed_transfer = 0.0
    for y in range(ny):
        for x in range(nx):
            donor_count = 0
            receiver_count = 0
            positive_total = 0.0
            for phase in range(phases):
                delta = after[phase, y, x] - before[phase, y, x]
                if delta < -1e-14:
                    donors[donor_count] = phase
                    donor_amount[donor_count] = -delta
                    donor_count += 1
                elif delta > 1e-14:
                    receivers[receiver_count] = phase
                    receiver_amount[receiver_count] = delta
                    positive_total += delta
                    receiver_count += 1
            if positive_total <= 0.0:
                continue
            absolute_sweep += positive_total * dx * dx
            for donor_index in range(donor_count):
                donor = donors[donor_index]
                for receiver_index in range(receiver_count):
                    receiver = receivers[receiver_index]
                    amount = (
                        donor_amount[donor_index]
                        * receiver_amount[receiver_index] / positive_total
                    )
                    if donor < receiver:
                        i, j, sign = donor, receiver, 1.0
                    else:
                        i, j, sign = receiver, donor, -1.0
                    b00, b11, b01 = _pair_tensor(
                        midpoint, i, j, y, x, orientations, dx
                    )
                    q = sign * amount
                    increment[0, 0, y, x] += q * b00
                    increment[1, 1, y, x] += q * b11
                    increment[0, 1, y, x] += q * b01
                    increment[1, 0, y, x] += q * b01
                    signed_transfer += q * dx * dx
    return increment, absolute_sweep, signed_transfer


@njit(cache=True)
def _map_local_transfer_sparse(
    before: NDArray[np.float64], after: NDArray[np.float64],
    orientations: NDArray[np.float64], dx: float, max_local_phases: int,
) -> tuple[NDArray[np.float64], float, float, int]:
    """Cache-friendly phase scan followed by small local pair loops."""
    phases, ny, nx = before.shape
    local_ids = np.full((max_local_phases, ny, nx), -1, dtype=np.int32)
    local_delta = np.zeros((max_local_phases, ny, nx), dtype=np.float64)
    local_count = np.zeros((ny, nx), dtype=np.int16)
    overflow = 0
    for phase in range(phases):
        for y in range(ny):
            for x in range(nx):
                delta = after[phase, y, x] - before[phase, y, x]
                if abs(delta) <= 1e-14:
                    continue
                slot = local_count[y, x]
                if slot >= max_local_phases:
                    overflow += 1
                    continue
                local_ids[slot, y, x] = phase
                local_delta[slot, y, x] = delta
                local_count[y, x] += 1
    increment = np.zeros((2, 2, ny, nx), dtype=np.float64)
    absolute_sweep = 0.0
    signed_transfer = 0.0
    for y in range(ny):
        for x in range(nx):
            count = local_count[y, x]
            positive_total = 0.0
            for slot in range(count):
                positive_total += max(local_delta[slot, y, x], 0.0)
            if positive_total <= 0.0:
                continue
            absolute_sweep += positive_total * dx * dx
            for donor_slot in range(count):
                donor_delta = local_delta[donor_slot, y, x]
                if donor_delta >= 0.0:
                    continue
                donor = local_ids[donor_slot, y, x]
                for receiver_slot in range(count):
                    receiver_delta = local_delta[receiver_slot, y, x]
                    if receiver_delta <= 0.0:
                        continue
                    receiver = local_ids[receiver_slot, y, x]
                    amount = -donor_delta * receiver_delta / positive_total
                    if donor < receiver:
                        i, j, sign = donor, receiver, 1.0
                    else:
                        i, j, sign = receiver, donor, -1.0
                    b00, b11, b01 = _pair_tensor_midpoint(
                        before, after, i, j, y, x, orientations, dx
                    )
                    q = sign * amount
                    increment[0, 0, y, x] += q * b00
                    increment[1, 1, y, x] += q * b11
                    increment[0, 1, y, x] += q * b01
                    increment[1, 0, y, x] += q * b01
                    signed_transfer += q * dx * dx
    return increment, absolute_sweep, signed_transfer, overflow


@njit(cache=True)
def _work_conjugate_driving(
    eta: NDArray[np.float64], orientations: NDArray[np.float64],
    stress: NDArray[np.float64], dx: float,
) -> NDArray[np.float64]:
    phases, ny, nx = eta.shape
    driving = np.zeros_like(eta)
    present = np.empty(phases, dtype=np.int64)
    for y in range(ny):
        ym, yp = (y - 1) % ny, (y + 1) % ny
        for x in range(nx):
            xm, xp = (x - 1) % nx, (x + 1) % nx
            count = 0
            for phase in range(phases):
                if (
                    eta[phase, y, x] > 1e-14
                    or eta[phase, ym, x] > 1e-14
                    or eta[phase, yp, x] > 1e-14
                    or eta[phase, y, xm] > 1e-14
                    or eta[phase, y, xp] > 1e-14
                ):
                    present[count] = phase
                    count += 1
            for left in range(count):
                i = present[left]
                for right in range(left + 1, count):
                    j = present[right]
                    b00, b11, b01 = _pair_tensor(
                        eta, i, j, y, x, orientations, dx
                    )
                    work = (
                        stress[0, 0, y, x] * b00
                        + stress[1, 1, y, x] * b11
                        + (stress[0, 1, y, x] + stress[1, 0, y, x]) * b01
                    )
                    # f=-delta E/deta, with q positive for i -> j.
                    driving[i, y, x] -= 0.5 * work
                    driving[j, y, x] += 0.5 * work
    return driving


@njit(cache=True)
def _work_conjugate_driving_sparse(
    eta: NDArray[np.float64], orientations: NDArray[np.float64],
    stress: NDArray[np.float64], dx: float, max_local_phases: int,
) -> tuple[NDArray[np.float64], int]:
    phases, ny, nx = eta.shape
    local_ids = np.full((max_local_phases, ny, nx), -1, dtype=np.int32)
    local_count = np.zeros((ny, nx), dtype=np.int16)
    overflow = 0
    for phase in range(phases):
        for y in range(ny):
            ym, yp = (y - 1) % ny, (y + 1) % ny
            for x in range(nx):
                xm, xp = (x - 1) % nx, (x + 1) % nx
                if (
                    eta[phase, y, x] <= 1e-14
                    and eta[phase, ym, x] <= 1e-14
                    and eta[phase, yp, x] <= 1e-14
                    and eta[phase, y, xm] <= 1e-14
                    and eta[phase, y, xp] <= 1e-14
                ):
                    continue
                slot = local_count[y, x]
                if slot >= max_local_phases:
                    overflow += 1
                    continue
                local_ids[slot, y, x] = phase
                local_count[y, x] += 1
    driving = np.zeros_like(eta)
    for y in range(ny):
        for x in range(nx):
            count = local_count[y, x]
            for left in range(count):
                i = local_ids[left, y, x]
                for right in range(left + 1, count):
                    j = local_ids[right, y, x]
                    b00, b11, b01 = _pair_tensor(
                        eta, i, j, y, x, orientations, dx
                    )
                    work = (
                        stress[0, 0, y, x] * b00
                        + stress[1, 1, y, x] * b11
                        + (stress[0, 1, y, x] + stress[1, 0, y, x]) * b01
                    )
                    driving[i, y, x] -= 0.5 * work
                    driving[j, y, x] += 0.5 * work
    return driving, overflow


def _is_exact_periodic_translation(
    before: NDArray[np.float64], after: NDArray[np.float64], max_shift: int
) -> bool:
    # Reject almost every evolving phase-field state using a small exact probe.
    # A true roll must pass every probe, after which the full array comparison
    # below remains the authoritative test.  This avoids repeatedly forming two
    # dense label maps and up to 24 rolled copies of a production-sized field.
    phases, ny, nx = before.shape
    phase_probes = sorted({0, phases // 3, (2 * phases) // 3, phases - 1})
    y_probes = sorted({0, ny // 3, (2 * ny) // 3, ny - 1})
    x_probes = sorted({0, nx // 3, (2 * nx) // 3, nx - 1})
    for shift_y in range(-max_shift, max_shift + 1):
        for shift_x in range(-max_shift, max_shift + 1):
            if shift_y == 0 and shift_x == 0:
                continue
            matches_probe = all(
                before[phase, (y - shift_y) % ny, (x - shift_x) % nx]
                == after[phase, y, x]
                for phase in phase_probes for y in y_probes for x in x_probes
            )
            if matches_probe and np.array_equal(
                np.roll(before, (shift_y, shift_x), axis=(-2, -1)), after
            ):
                return True
    return False


@dataclass(frozen=True)
class SweepDiagnostics:
    absolute_swept_area: float
    signed_pair_transfer: float
    integrated_source: tuple[tuple[float, float], tuple[float, float]]
    rigid_translation_removed: bool


class LocalInterfaceSweepCoupling:
    """Unique phase-transfer source and its discrete work-conjugate force."""

    def __init__(self, grid_spacing: float = 1.0, *, rigid_shift_search: int = 2,
                 max_local_phases: int = 8):
        if grid_spacing <= 0.0 or rigid_shift_search < 0 or max_local_phases < 2:
            raise ValueError("invalid local-sweep discretization")
        self.grid_spacing = float(grid_spacing)
        self.rigid_shift_search = int(rigid_shift_search)
        self.max_local_phases = int(max_local_phases)
        self.last_diagnostics = SweepDiagnostics(0.0, 0.0, ((0.0, 0.0), (0.0, 0.0)), False)

    def source_increment(
        self, before: NDArray[np.float64], after: NDArray[np.float64],
        orientations: NDArray[np.float64],
    ) -> NDArray[np.float64]:
        before_value = np.asarray(before, dtype=float)
        after_value = np.asarray(after, dtype=float)
        orientation_value = np.asarray(orientations, dtype=float)
        if before_value.shape != after_value.shape or before_value.ndim != 3:
            raise ValueError("before and after must have the same phase-field shape")
        if orientation_value.shape != (before_value.shape[0],):
            raise ValueError("one orientation is required for every phase")
        rigid = (
            self.rigid_shift_search > 0
            and _is_exact_periodic_translation(
                before_value, after_value, self.rigid_shift_search
            )
        )
        if rigid:
            increment = np.zeros((2, 2, *before_value.shape[1:]), dtype=float)
            swept, signed = 0.0, 0.0
        else:
            increment, swept, signed, overflow = _map_local_transfer_sparse(
                before_value, after_value, orientation_value, self.grid_spacing,
                self.max_local_phases,
            )
            if overflow:
                raise RuntimeError(
                    f"local source support exceeded {self.max_local_phases} phases at {overflow} insertions"
                )
        integrated = np.sum(increment, axis=(-2, -1)) * self.grid_spacing**2
        self.last_diagnostics = SweepDiagnostics(
            float(swept), float(signed),
            ((float(integrated[0, 0]), float(integrated[0, 1])),
             (float(integrated[1, 0]), float(integrated[1, 1]))),
            rigid,
        )
        return increment

    def driving_field(
        self, eta: NDArray[np.float64], orientations: NDArray[np.float64],
        stress: NDArray[np.float64],
    ) -> NDArray[np.float64]:
        driving, overflow = _work_conjugate_driving_sparse(
            np.asarray(eta, dtype=float), np.asarray(orientations, dtype=float),
            np.asarray(stress, dtype=float), self.grid_spacing,
            self.max_local_phases,
        )
        if overflow:
            raise RuntimeError(
                f"local driving support exceeded {self.max_local_phases} phases at {overflow} insertions"
            )
        return driving
