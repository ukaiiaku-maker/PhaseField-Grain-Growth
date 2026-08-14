def boundary_measure_change(previous: float, current: float) -> float:
    return abs(float(current) - float(previous))


def point_defect_requirement(area_change: float, excess_volume_per_area: float,
                             point_defect_formation_volume: float) -> float:
    if point_defect_formation_volume <= 0 or excess_volume_per_area < 0:
        raise ValueError("formation volume must be positive and excess volume nonnegative")
    return excess_volume_per_area * abs(area_change) / point_defect_formation_volume

