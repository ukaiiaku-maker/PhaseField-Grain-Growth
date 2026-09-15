from __future__ import annotations

from dataclasses import dataclass
import math


@dataclass(frozen=True)
class NormalForceBalance:
    """Signed sharp-interface diagnostic for the local normal driving.

    ``p_cap`` is the sharp-interface estimate of the capillary contribution.
    The PF solver evaluates capillarity from its diffuse free-energy functional;
    this scalar is therefore diagnostic rather than an extra applied force.
    ``p_event`` is kept separate because the migration-closure release impulse
    is not a climb chemical pressure.
    """

    p_cap: float
    p_chem: float
    p_shear: float
    p_event: float
    p_net: float
    p_applied_total: float
    chi_s: float


def normal_force_balance(
    *,
    capillary_pressure: float,
    chemical_pressure: float,
    beta: float,
    resolved_shear: float,
    event_pressure: float = 0.0,
) -> NormalForceBalance:
    """Compose the signed local driving without changing the PF evolution.

    The mechanical arrest ratio is ``-p_shear / (p_cap + p_chem)``.  It is
    positive when shear opposes the non-shear drive, equals one at force
    balance, and exceeds one after local drive reversal.  A NaN is returned
    when the non-shear denominator vanishes.
    """

    values = (
        capillary_pressure, chemical_pressure, beta, resolved_shear,
        event_pressure,
    )
    if not all(math.isfinite(float(value)) for value in values):
        raise ValueError("normal-force inputs must be finite")
    p_cap = float(capillary_pressure)
    p_chem = float(chemical_pressure)
    p_shear = float(beta) * float(resolved_shear)
    p_event = float(event_pressure)
    non_shear = p_cap + p_chem
    chi_s = -p_shear / non_shear if abs(non_shear) > 1e-14 else math.nan
    return NormalForceBalance(
        p_cap=p_cap,
        p_chem=p_chem,
        p_shear=p_shear,
        p_event=p_event,
        p_net=non_shear + p_shear,
        p_applied_total=non_shear + p_shear + p_event,
        chi_s=chi_s,
    )
