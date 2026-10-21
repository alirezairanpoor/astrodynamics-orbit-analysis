from math import atan2, sqrt, sin, cos, pi

import numpy as np


def states_to_elements(
    r_eci: np.ndarray,
    v_eci: np.ndarray,
    mu: float,
) -> tuple[float, float, float, float, float, float]:
    """
    Convert ECI position and velocity vectors
    to classical Keplerian orbital elements.

    Returns
    -------
    tuple
        a    : Semi-major axis [km]
        e    : Eccentricity [-]
        i    : Inclination [rad]
        RAAN : Right ascension of ascending node [rad]
        AOP  : Argument of perigee [rad]
        M    : Mean anomaly [rad]
    """

    r_mag = np.linalg.norm(r_eci)
    v_mag = np.linalg.norm(v_eci)

    h_vec = np.cross(r_eci, v_eci)
    h = np.linalg.norm(h_vec)

    W = h_vec / h

    i = atan2(
        sqrt(W[0] ** 2 + W[1] ** 2),
        W[2],
    )

    RAAN = atan2(
        W[0],
        -W[1],
    ) % (2 * pi)

    p = h ** 2 / mu

    a = 1 / (
        2 / r_mag
        - v_mag ** 2 / mu
    )

    e = sqrt(
        1 - p / a
    )

    n = sqrt(
        mu / a ** 3
    )

    E = atan2(
        np.dot(r_eci, v_eci) / (a ** 2 * n),
        1 - r_mag / a,
    ) % (2 * pi)

    M = (
        E
        - e * sin(E)
    ) % (2 * pi)

    argument_of_latitude = atan2(
        r_eci[2] / sin(i),
        r_eci[0] * cos(RAAN)
        + r_eci[1] * sin(RAAN),
    ) % (2 * pi)

    TA = atan2(
        sqrt(1 - e ** 2) * sin(E),
        cos(E) - e,
    ) % (2 * pi)

    AOP = (
        argument_of_latitude
        - TA
    ) % (2 * pi)

    return (
        a,
        e,
        i,
        RAAN,
        AOP,
        M,
    )
