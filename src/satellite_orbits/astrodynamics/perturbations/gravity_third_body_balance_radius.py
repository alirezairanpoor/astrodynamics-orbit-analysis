from math import hypot


def gravity_third_body_balance_radius(
    n: int,
    C_normalized: float,
    S_normalized: float,
    mu_earth: float,
    earth_radius: float,
    mu_third_body: float,
    third_body_distance: float,
) -> float:
    """
    Estimate the geocentric radius where one gravity harmonic
    and the third-body acceleration have equal magnitudes.

    Uses the simplified relations from Exercise 3.3.

    Inputs:
        n: harmonic degree
        C_normalized, S_normalized: normalized coefficients
        mu_earth, mu_third_body: km^3/s^2
        earth_radius, third_body_distance: km

    Returns:
        Geocentric balance radius [km], not altitude.
    """

    coefficient_magnitude = hypot(
        C_normalized,
        S_normalized,
    )

    r = (
        (n + 1) / 2.0
        * (mu_earth / mu_third_body)
        * earth_radius**n
        * third_body_distance**3
        * coefficient_magnitude
    ) ** (1.0 / (n + 3))

    return r
