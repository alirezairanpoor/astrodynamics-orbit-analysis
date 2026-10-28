from math import sqrt

from satellite_orbits.astrodynamics.keplerian.two_body import circular_velocity


def hohmann_transfer(mu: float, r1: float, r2: float) -> tuple[float, float, float]:
    """
    Compute the delta-v requirements for a Hohmann transfer
    between two circular orbits.

    Parameters
    ----------
    mu : float
        Gravitational parameter [km^3/s^2].
    r1 : float
        Initial orbital radius [km].
    r2 : float
        Final orbital radius [km].

    Returns
    -------
    tuple[float, float, float]
        delta_v1 : First burn [km/s].
        delta_v2 : Second burn [km/s].
        delta_v_total : Total delta-v [km/s].
    """

    a_transfer = (r1 + r2) / 2

    v1 = circular_velocity(mu, r1)
    v2 = circular_velocity(mu, r2)

    v_perigee = sqrt(mu * (2 / r1 - 1 / a_transfer))
    v_apogee = sqrt(mu * (2 / r2 - 1 / a_transfer))

    delta_v1 = v_perigee - v1
    delta_v2 = v2 - v_apogee
    delta_v_total = delta_v1 + delta_v2

    return delta_v1, delta_v2, delta_v_total