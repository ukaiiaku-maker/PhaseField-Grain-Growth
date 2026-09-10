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


def _is_exact_periodic_translation(
    before: NDArray[np.float64], after: NDArray[np.float64], max_shift: int
) -> bool:
    before_labels = np.argmax(before, axis=0)
    after_labels = np.argmax(after, axis=0)
    for shift_y in range(-max_shift, max_shift + 1):
        for shift_x in range(-max_shift, max_shift + 1):
            if shift_y == 0 and shift_x == 0:
                continue
            if np.array_equal(
                np.roll(before_labels, (shift_y, shift_x), axis=(0, 1)), after_labels
            ) and np.array_equal(
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

    def __init__(self, grid_spacing: float = 1.0, *, rigid_shift_search: int = 2):
        if grid_spacing <= 0.0 or rigid_shift_search < 0:
            raise ValueError("invalid local-sweep discretization")
        self.grid_spacing = float(grid_spacing)
        self.rigid_shift_search = int(rigid_shift_search)
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
            increment, swept, signed = _map_local_transfer(
                before_value, after_value, orientation_value, self.grid_spacing
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
        return _work_conjugate_driving(
            np.asarray(eta, dtype=float), np.asarray(orientations, dtype=float),
            np.asarray(stress, dtype=float), self.grid_spacing,
        )
