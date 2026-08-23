from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Iterable

import numpy as np


def gb_excess_volume_density(
    *,
    excess_volume_per_area: float | None = None,
    delta_gb: float | None = None,
    rho_lattice: float | None = None,
    rho_gb: float | None = None,
) -> float:
    """Return GB excess volume per physical GB measure.

    In the two-dimensional solver the GB measure is length per unit
    out-of-plane thickness.  A directly supplied value is preferred.  The
    density form is retained for configurations that explicitly provide all
    three quantities.
    """
    if excess_volume_per_area is not None:
        value = float(excess_volume_per_area)
    else:
        if delta_gb is None or rho_lattice is None or rho_gb is None:
            raise ValueError(
                "supply excess_volume_per_area or delta_gb/rho_lattice/rho_gb"
            )
        if float(rho_gb) == 0.0:
            raise ValueError("rho_gb must be nonzero")
        value = float(delta_gb) * (1.0 - float(rho_lattice) / float(rho_gb))
    if not np.isfinite(value) or value < 0.0:
        raise ValueError("GB excess-volume density must be finite and nonnegative")
    return value


def released_excess_volume(
    previous_gb_measure: float,
    current_gb_measure: float,
    excess_volume_density: float,
    alpha: float = 1.0,
) -> float:
    """Excess volume released only by a decrease in total GB measure."""
    values = np.asarray(
        (previous_gb_measure, current_gb_measure, excess_volume_density, alpha),
        dtype=float,
    )
    if not np.all(np.isfinite(values)) or np.any(values < 0.0):
        raise ValueError("GB measures, excess-volume density, and alpha must be finite and nonnegative")
    return float(alpha * excess_volume_density * max(previous_gb_measure - current_gb_measure, 0.0))


def released_point_defect_quota(
    previous_gb_measure: float,
    current_gb_measure: float,
    excess_volume_density: float,
    formation_volume: float,
    alpha: float = 1.0,
) -> float:
    """Positive vacancy-equivalent quota released by GB measure loss."""
    if not np.isfinite(formation_volume) or formation_volume <= 0.0:
        raise ValueError("point-defect formation volume must be finite and positive")
    return released_excess_volume(
        previous_gb_measure, current_gb_measure, excess_volume_density, alpha
    ) / float(formation_volume)


def best_compatible_tj_velocity(
    normals: np.ndarray,
    desired_normal_velocities: np.ndarray,
    weights: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Least-squares common TJ velocity and signed branch residuals."""
    n = np.asarray(normals, dtype=float)
    desired = np.asarray(desired_normal_velocities, dtype=float)
    if n.ndim != 2 or n.shape[1] != 2 or desired.shape != (n.shape[0],):
        raise ValueError("normals must be (branches,2) and velocities must match")
    if n.shape[0] < 2 or not np.all(np.isfinite(n)) or not np.all(np.isfinite(desired)):
        raise ValueError("finite data for at least two TJ branches are required")
    norms = np.linalg.norm(n, axis=1)
    if np.any(norms <= np.finfo(float).tiny):
        raise ValueError("TJ branch normals must be nonzero")
    n = n / norms[:, None]
    w = np.ones(n.shape[0], dtype=float) if weights is None else np.asarray(weights, dtype=float)
    if w.shape != (n.shape[0],) or np.any(~np.isfinite(w)) or np.any(w <= 0.0):
        raise ValueError("TJ compatibility weights must be finite and positive")
    root_w = np.sqrt(w)
    velocity, *_ = np.linalg.lstsq(n * root_w[:, None], desired * root_w, rcond=None)
    residual = desired - n @ velocity
    norm = float(np.sqrt(np.sum(w * residual * residual)))
    return velocity.astype(float), residual.astype(float), norm


@dataclass(frozen=True)
class TJSinkDecision:
    allowed: bool
    reason: str
    flux_sign: int
    compatibility_norm: float
    residual_energy_change: float


def validate_tj_sink_candidate(
    *,
    available_signed_quota: float,
    requested_sign: int,
    compatibility_norm: float,
    compatibility_tolerance: float,
    residual_burgers: np.ndarray,
    burgers_increment: np.ndarray,
    residual_stiffness: float,
    strict_burgers_tolerance: float | None = None,
) -> TJSinkDecision:
    """Validate flux direction, geometry, and Burgers accounting for a TJ sink."""
    sign = int(np.sign(requested_sign))
    if sign == 0 or available_signed_quota * sign <= 0.0:
        return TJSinkDecision(False, "wrong_sign_or_no_flux", sign, float(compatibility_norm), 0.0)
    if compatibility_norm > compatibility_tolerance:
        return TJSinkDecision(False, "incompatible_geometry", sign, float(compatibility_norm), 0.0)
    before = np.asarray(residual_burgers, dtype=float)
    increment = np.asarray(burgers_increment, dtype=float)
    if before.shape != (2,) or increment.shape != (2,):
        raise ValueError("TJ residual and Burgers increment must be 2-vectors")
    if not np.isfinite(residual_stiffness) or residual_stiffness < 0.0:
        raise ValueError("residual stiffness must be finite and nonnegative")
    after = before + increment
    if strict_burgers_tolerance is not None and np.linalg.norm(after) > strict_burgers_tolerance:
        return TJSinkDecision(False, "burgers_incompatible", sign, float(compatibility_norm), np.inf)
    energy = 0.5 * residual_stiffness * float(after @ after - before @ before)
    return TJSinkDecision(True, "allowed_finite_residual", sign, float(compatibility_norm), energy)


@dataclass
class ConservedDefectInventory:
    """Material-wide vacancy/interstitial ledger with topology-safe storage."""

    required_vacancy: float = 0.0
    required_interstitial: float = 0.0
    accommodated_gb_vacancy: float = 0.0
    accommodated_gb_interstitial: float = 0.0
    accommodated_tj_vacancy: float = 0.0
    accommodated_tj_interstitial: float = 0.0
    accommodated_external_vacancy: float = 0.0
    accommodated_external_interstitial: float = 0.0
    active_vacancy: dict[str, float] = field(default_factory=dict)
    active_interstitial: dict[str, float] = field(default_factory=dict)
    retired_vacancy: float = 0.0
    retired_interstitial: float = 0.0

    @property
    def required_total(self) -> float:
        return self.required_vacancy + self.required_interstitial

    @property
    def accommodated_gb(self) -> float:
        return self.accommodated_gb_vacancy + self.accommodated_gb_interstitial

    @property
    def accommodated_tj(self) -> float:
        return self.accommodated_tj_vacancy + self.accommodated_tj_interstitial

    @property
    def accommodated_external(self) -> float:
        return self.accommodated_external_vacancy + self.accommodated_external_interstitial

    @property
    def active_deficit(self) -> float:
        return sum(self.active_vacancy.values()) + sum(self.active_interstitial.values())

    @property
    def retired_inventory(self) -> float:
        return self.retired_vacancy + self.retired_interstitial

    @property
    def stored_total(self) -> float:
        return self.active_deficit + self.retired_inventory

    @property
    def stored_signed(self) -> float:
        return (
            sum(self.active_vacancy.values()) + self.retired_vacancy
            - sum(self.active_interstitial.values()) - self.retired_interstitial
        )

    @property
    def conservation_residual(self) -> float:
        return self.required_total - (
            self.accommodated_gb + self.accommodated_tj
            + self.accommodated_external + self.stored_total
        )

    @property
    def vacancy_residual(self) -> float:
        return self.required_vacancy - (
            self.accommodated_gb_vacancy + self.accommodated_tj_vacancy
            + self.accommodated_external_vacancy + sum(self.active_vacancy.values())
            + self.retired_vacancy
        )

    @property
    def interstitial_residual(self) -> float:
        return self.required_interstitial - (
            self.accommodated_gb_interstitial + self.accommodated_tj_interstitial
            + self.accommodated_external_interstitial + sum(self.active_interstitial.values())
            + self.retired_interstitial
        )

    def require(self, signed_quota: float, entity_id: str | None = None) -> float:
        value = float(signed_quota)
        if not np.isfinite(value):
            raise ValueError("required point-defect quota must be finite")
        amount = abs(value)
        if amount == 0.0:
            return 0.0
        key = entity_id or "material-reservoir"
        if value > 0.0:
            self.required_vacancy += amount
            self.active_vacancy[key] = self.active_vacancy.get(key, 0.0) + amount
        else:
            self.required_interstitial += amount
            self.active_interstitial[key] = self.active_interstitial.get(key, 0.0) + amount
        self.assert_conserved()
        return value

    def retire_missing(self, active_entity_ids: Iterable[str]) -> None:
        active = set(active_entity_ids)
        for mapping, retired_name in (
            (self.active_vacancy, "retired_vacancy"),
            (self.active_interstitial, "retired_interstitial"),
        ):
            for key in sorted(set(mapping) - active - {"material-reservoir"}):
                setattr(self, retired_name, getattr(self, retired_name) + mapping.pop(key))
        self.assert_conserved()

    def split(self, parent: str, children: Iterable[str], weights: Iterable[float] | None = None) -> None:
        child_ids = list(children)
        if not child_ids:
            raise ValueError("a split requires at least one child")
        fractions = np.ones(len(child_ids), dtype=float) if weights is None else np.asarray(list(weights), dtype=float)
        if fractions.shape != (len(child_ids),) or np.any(fractions < 0.0) or not np.any(fractions > 0.0):
            raise ValueError("split weights must be nonnegative with positive sum")
        fractions /= fractions.sum()
        for mapping in (self.active_vacancy, self.active_interstitial):
            amount = mapping.pop(parent, 0.0)
            for child, fraction in zip(child_ids, fractions, strict=True):
                mapping[child] = mapping.get(child, 0.0) + amount * float(fraction)
        self.assert_conserved()

    def merge(self, parents: Iterable[str], child: str) -> None:
        for mapping in (self.active_vacancy, self.active_interstitial):
            amount = sum(mapping.pop(parent, 0.0) for parent in parents)
            mapping[child] = mapping.get(child, 0.0) + amount
        self.assert_conserved()

    def _consume_species(self, amount: float, mapping: dict[str, float], retired_name: str,
                         preferred_entity: str | None) -> float:
        remaining = amount
        keys = sorted(mapping)
        if preferred_entity in mapping:
            keys.remove(preferred_entity)
            keys.insert(0, preferred_entity)
        for key in keys:
            accepted = min(mapping[key], remaining)
            mapping[key] -= accepted
            remaining -= accepted
            if mapping[key] <= 16 * np.finfo(float).eps:
                mapping.pop(key)
            if remaining <= 16 * np.finfo(float).eps:
                return amount
        retired = getattr(self, retired_name)
        accepted = min(retired, remaining)
        setattr(self, retired_name, retired - accepted)
        return amount - (remaining - accepted)

    def accommodate(self, signed_quota: float, sink: str,
                    preferred_entity: str | None = None) -> float:
        if sink not in {"gb", "tj", "external"}:
            raise ValueError("sink must be gb, tj, or external")
        requested = float(signed_quota)
        if not np.isfinite(requested):
            raise ValueError("sink quota must be finite")
        if requested == 0.0:
            return 0.0
        vacancy = requested > 0.0
        mapping = self.active_vacancy if vacancy else self.active_interstitial
        retired_name = "retired_vacancy" if vacancy else "retired_interstitial"
        accepted = self._consume_species(abs(requested), mapping, retired_name, preferred_entity)
        field_name = f"accommodated_{sink}_{'vacancy' if vacancy else 'interstitial'}"
        setattr(self, field_name, getattr(self, field_name) + accepted)
        self.assert_conserved()
        return accepted if vacancy else -accepted

    def assert_conserved(self, tolerance: float = 1e-10) -> None:
        scale = max(1.0, self.required_total)
        if (
            abs(self.conservation_residual) > tolerance * scale
            or abs(self.vacancy_residual) > tolerance * scale
            or abs(self.interstitial_residual) > tolerance * scale
        ):
            raise RuntimeError(
                "point-defect conservation failure: "
                f"total={self.conservation_residual:g}, vacancy={self.vacancy_residual:g}, "
                f"interstitial={self.interstitial_residual:g}"
            )

    def state_dict(self) -> dict[str, Any]:
        return {
            name: value
            for name, value in vars(self).items()
        }

    @classmethod
    def from_state_dict(cls, state: dict[str, Any]) -> "ConservedDefectInventory":
        result = cls(**state)
        result.assert_conserved()
        return result

    def diagnostics(self) -> dict[str, float]:
        accommodated = self.accommodated_gb + self.accommodated_tj + self.accommodated_external
        return {
            "N_required": self.required_total,
            "N_required_signed": self.required_vacancy - self.required_interstitial,
            "N_accommodated_GB": self.accommodated_gb,
            "N_accommodated_TJ": self.accommodated_tj,
            "N_external": self.accommodated_external,
            "N_active_deficit": self.active_deficit,
            "N_retired": self.retired_inventory,
            "N_stored_signed": self.stored_signed,
            "conservation_residual": self.conservation_residual,
            "vacancy_conservation_residual": self.vacancy_residual,
            "interstitial_conservation_residual": self.interstitial_residual,
            "GB_accommodation_fraction": self.accommodated_gb / accommodated if accommodated else 0.0,
            "TJ_accommodation_fraction": self.accommodated_tj / accommodated if accommodated else 0.0,
            "external_accommodation_fraction": self.accommodated_external / accommodated if accommodated else 0.0,
            "active_deficit_fraction": self.active_deficit / self.required_total if self.required_total else 0.0,
            "retired_inventory_fraction": self.retired_inventory / self.required_total if self.required_total else 0.0,
        }
