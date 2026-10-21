from math import sqrt, acos


def sun_synchronous_inclination(
    a: float,
    RAAN_rate_target: float,
    mu: float,
    earth_radius: float,
    J2: float,
) -> float:
    """
    Compute the inclination required for a
    sun-synchronous circular orbit.

    Parameters
    ----------
    a : float
        Semi-major axis [km].

    RAAN_rate_target : float
        Required RAAN rate [rad/s].

    mu : float
        Gravitational parameter [km^3/s^2].

    earth_radius : float
        Earth radius [km].

    J2 : float
        Earth's second zonal harmonic coefficient [-].

    Returns
    -------
    float
        Required inclination i [rad].
    """

    n0 = sqrt(
        mu / a ** 3
    )

    coefficient = (
        -3 / 2
        * n0
        * J2
        * (earth_radius / a) ** 2
    )

    cos_i = (
        RAAN_rate_target
        / coefficient
    )

    i = acos(cos_i)

    return i
