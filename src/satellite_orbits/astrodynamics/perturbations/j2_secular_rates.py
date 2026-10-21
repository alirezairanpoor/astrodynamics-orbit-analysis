from math import sqrt, cos


def j2_secular_rates(
    a: float,
    i: float,
    mu: float,
    earth_radius: float,
    J2: float,
) -> tuple[float, float, float, float]:
    """
    Compute the first-order J2 secular rates
    for a circular Earth orbit.

    Parameters
    ----------
    a : float
        Semi-major axis [km].

    i : float
        Inclination [rad].

    mu : float
        Gravitational parameter [km^3/s^2].

    earth_radius : float
        Earth radius [km].

    J2 : float
        Earth's second zonal harmonic coefficient [-].

    Returns
    -------
    tuple
        n0        : Keplerian mean motion [rad/s].
        RAAN_rate : Secular RAAN rate [rad/s].
        AOP_rate  : Secular AOP rate [rad/s].
        delta_n   : J2 correction to mean motion [rad/s].
    """

    n0 = sqrt(
        mu / a ** 3
    )

    radius_ratio_squared = (
        earth_radius / a
    ) ** 2

    RAAN_rate = (
        -3 / 2
        * n0
        * J2
        * radius_ratio_squared
        * cos(i)
    )

    AOP_rate = (
        -3 / 4
        * n0
        * J2
        * radius_ratio_squared
        * (
            1
            - 5 * cos(i) ** 2
        )
    )

    delta_n = (
        -3 / 4
        * n0
        * J2
        * radius_ratio_squared
        * (
            1
            - 3 * cos(i) ** 2
        )
    )

    return (
        n0,
        RAAN_rate,
        AOP_rate,
        delta_n,
    )
