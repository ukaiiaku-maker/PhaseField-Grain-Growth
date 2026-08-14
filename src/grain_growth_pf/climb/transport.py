from grain_growth_pf.disconnections.mode import K_B_EV

import numpy as np


def diffusivity(temperature: float, prefactor: float, activation_energy_ev: float) -> float:
    if temperature <= 0 or prefactor < 0 or activation_energy_ev < 0:
        raise ValueError("invalid diffusivity parameters")
    return float(prefactor * np.exp(-activation_energy_ev / (K_B_EV * temperature)))


def transport_time(length: float, diffusion_coefficient: float, geometry_factor: float = 1.0) -> float:
    if length < 0 or diffusion_coefficient <= 0 or geometry_factor <= 0:
        raise ValueError("invalid transport parameters")
    return geometry_factor * length**2 / diffusion_coefficient

