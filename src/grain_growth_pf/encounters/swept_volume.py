def swept_measure(interface_measure: float, normal_displacement: float,
                  physics_dimension: int = 2, out_of_plane_thickness: float = 1.0) -> float:
    if physics_dimension == 2:
        return abs(interface_measure * normal_displacement)
    if physics_dimension == 3:
        return abs(interface_measure * normal_displacement * out_of_plane_thickness)
    raise ValueError("physics_dimension must be 2 or 3")

