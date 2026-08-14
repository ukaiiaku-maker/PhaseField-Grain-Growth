from grain_growth_pf.disconnections.mode import DisconnectionMode, ModeDriving


def nucleation_rate(mode: DisconnectionMode, temperature: float, driving: ModeDriving) -> float:
    return mode.rate(temperature, driving)

