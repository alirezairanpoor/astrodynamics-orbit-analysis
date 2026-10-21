from collections.abc import Callable


def drag_srp_balance_altitude(
    density_function: Callable[[float], float],
    h_lower: float,
    h_upper: float,
    C_D: float,
    C_R: float,
    solar_pressure: float,
    mu: float,
    earth_radius: float,
    tolerance_km: float = 1e-3,
) -> float:
    """
    Find the altitude where simplified drag and SRP magnitudes
    are equal, using bisection.

    Assumptions:
        Circular orbit.
        Full sunlight.
        Atmospheric rotation and winds neglected.
        Same effective area in both acceleration models.

    Inputs:
        density_function: altitude [km] -> density [kg/m^3]
        h_lower, h_upper: search interval [km]
        C_D, C_R: dimensionless coefficients
        solar_pressure: radiation pressure [N/m^2]
        mu: Earth gravitational parameter [km^3/s^2]
        earth_radius: Earth radius [km]
        tolerance_km: final interval width [km]

    Output:
        Balance altitude [km].
    """

    if h_upper <= h_lower:
        raise ValueError("Require h_upper > h_lower.")

    if tolerance_km <= 0.0:
        raise ValueError("tolerance_km must be positive.")

    # Use SI internally because density and pressure are in SI.
    mu_si = mu * 1e9

    def residual(h: float) -> float:
        r_m = (earth_radius + h) * 1000.0

        required_density = (
            2.0 * C_R * solar_pressure * r_m
            / (C_D * mu_si)
        )

        return density_function(h) - required_density

    f_lower = residual(h_lower)
    f_upper = residual(h_upper)

    if f_lower == 0.0:
        return h_lower

    if f_upper == 0.0:
        return h_upper

    if f_lower * f_upper > 0.0:
        raise ValueError(
            "Search interval must bracket a balance altitude."
        )

    while h_upper - h_lower > tolerance_km:
        h_mid = 0.5 * (h_lower + h_upper)
        f_mid = residual(h_mid)

        if f_mid == 0.0:
            return h_mid

        if f_lower * f_mid < 0.0:
            h_upper = h_mid
        else:
            h_lower = h_mid
            f_lower = f_mid

    return 0.5 * (h_lower + h_upper)
