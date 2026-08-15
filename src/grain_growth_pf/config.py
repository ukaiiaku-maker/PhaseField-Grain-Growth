from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import yaml


@dataclass(frozen=True)
class PFConfig:
    simulation_dimension: int = 2
    physics_dimension: int = 2
    shape: tuple[int, int] = (96, 96)
    grid_spacing: float = 1.0
    interface_width: float = 4.0
    time_step: float = 0.02
    gb_energy: float = 1.0
    intrinsic_mobility: float = 1.0
    boundary_conditions: str = "periodic"
    temperature: float = 900.0
    adaptive_stepping: bool = False
    grain_extinction_threshold: float = 0.5

    def __post_init__(self) -> None:
        if self.simulation_dimension != 2:
            raise ValueError("Only two-dimensional PF geometry is implemented")
        if self.physics_dimension not in (2, 3):
            raise ValueError("physics_dimension must be 2 or 3")
        if self.grid_spacing <= 0 or self.interface_width <= 0 or self.time_step <= 0:
            raise ValueError("dx, interface width, and dt must be positive")
        if self.boundary_conditions not in {"periodic", "neumann"}:
            raise ValueError("boundary_conditions must be periodic or neumann")
        if not 0 < self.grain_extinction_threshold < 1:
            raise ValueError("grain_extinction_threshold must lie in (0,1)")


@dataclass(frozen=True)
class ModelConfig:
    regime: str = "B0"
    seed: int = 1
    pf: PFConfig = field(default_factory=PFConfig)
    mechanics_backend: str = "none"
    compatibility_model: str = "off"
    active_modules: tuple[str, ...] = ()
    output_cadence: int = 20
    max_steps: int = 1000
    termination_grains: int = 8
    parameters: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.mechanics_backend not in {"none", "local_memory", "qiu_full_field"}:
            raise ValueError("unknown mechanics backend")
        if self.compatibility_model not in {"off", "explicit_modes", "geometric_surrogate"}:
            raise ValueError("unknown compatibility model")

    def to_dict(self) -> dict[str, Any]:
        data = asdict(self)
        data["pf"]["shape"] = list(self.pf.shape)
        data["active_modules"] = list(self.active_modules)
        return data

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "ModelConfig":
        d = dict(data)
        pf_data = dict(d.pop("pf", {}))
        if "shape" in pf_data:
            pf_data["shape"] = tuple(pf_data["shape"])
        if "active_modules" in d:
            d["active_modules"] = tuple(d["active_modules"])
        if d.get("compatibility_model") is False:
            d["compatibility_model"] = "off"
        return cls(pf=PFConfig(**pf_data), **d)

    @classmethod
    def load(cls, path: str | Path) -> "ModelConfig":
        with Path(path).open(encoding="utf-8") as handle:
            return cls.from_dict(yaml.safe_load(handle))
