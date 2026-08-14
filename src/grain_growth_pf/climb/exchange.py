from __future__ import annotations

import numpy as np

from grain_growth_pf.disconnections.mode import K_B_EV


def butler_volmer_flux(delta_mu_ev: float, temperature: float, exchange_rate: float,
                       alpha: float = 0.5) -> float:
    if temperature <= 0 or exchange_rate < 0 or not 0 <= alpha <= 1:
        raise ValueError("invalid Butler-Volmer parameters")
    x = delta_mu_ev / (K_B_EV * temperature)
    # Defined high-driving limit: each exponential is evaluated up to the
    # floating-point logarithmic limit and reported as finite saturated flux.
    limit = np.log(np.finfo(float).max / max(exchange_rate, 1.0)) - 2.0
    forward = np.exp(min(alpha * x, limit))
    reverse = np.exp(min(-(1 - alpha) * x, limit))
    return float(exchange_rate * (forward - reverse))


def linear_onsager_coefficient(temperature: float, exchange_rate: float) -> float:
    return exchange_rate / (K_B_EV * temperature)

