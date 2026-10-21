from math import atan2, sqrt, sin, cos, pi

import numpy as np

from .gauss_sector_triangle_ratio import (
    gauss_sector_triangle_ratio,
)


def gauss_two_position_elements(
    r_a_eci: np.ndarray,
    r_b_eci: np.ndarray,
    delta_t: float,
    mu: float,
) -> tuple[float, float, float, float, float, float]:
    """
    Determine Keplerian orbital elements from two
    ECI position vectors using Gauss's method.

    Parameters
    ----------
    r_a_eci : np.ndarray
        First ECI position vector [km].

    r_b_eci : np.ndarray
        Second ECI position vector [km].

    delta_t : float
        Time from the first to the second position [s].

    mu : float
        Gravitational parameter [km^3/s^2].

    Returns
    -------
    tuple
        a    : Semi-major axis [km].
        e    : Eccentricity [-].
        i    : Inclination [rad].
        RAAN : Right ascension of ascending node [rad].
        AOP  : Argument of perigee [rad].
        M    : Mean anomaly at the first position [rad].
    """

    # ========================================================
    # POSITION MAGNITUDES
    # ========================================================

    r_a = np.linalg.norm(r_a_eci)
    r_b = np.linalg.norm(r_b_eci)


    # ========================================================
    # ORBITAL PLANE
    # ========================================================

    e_a = (
        r_a_eci
        / r_a
    )

    r_0_eci = (
        r_b_eci
        - np.dot(
            r_b_eci,
            e_a,
        )
        * e_a
    )

    r_0 = np.linalg.norm(
        r_0_eci
    )

    e_0 = (
        r_0_eci
        / r_0
    )

    W = np.cross(
        e_a,
        e_0,
    )


    # ========================================================
    # INCLINATION AND RAAN
    # ========================================================

    i = atan2(
        sqrt(
            W[0] ** 2
            + W[1] ** 2
        ),
        W[2],
    )

    RAAN = atan2(
        W[0],
        -W[1],
    ) % (2 * pi)


    # ========================================================
    # ARGUMENT OF LATITUDE AT FIRST POSITION
    # ========================================================

    argument_of_latitude = atan2(
        r_a_eci[2],
        (
            -r_a_eci[0] * W[1]
            + r_a_eci[1] * W[0]
        ),
    ) % (2 * pi)


    # ========================================================
    # GAUSS SECTOR-TO-TRIANGLE RATIO
    # ========================================================

    eta = gauss_sector_triangle_ratio(
        r_a_eci,
        r_b_eci,
        delta_t,
        mu,
    )


    # ========================================================
    # SEMI-LATUS RECTUM
    # ========================================================

    tau = (
        sqrt(mu)
        * delta_t
    )

    triangle_area = (
        0.5
        * r_a
        * r_0
    )

    p = (
        2
        * triangle_area
        * eta
        / tau
    ) ** 2


    # ========================================================
    # ECCENTRICITY AND TRUE ANOMALY
    # AT FIRST POSITION
    # ========================================================

    e_cos_TA = (
        p / r_a
        - 1
    )

    cos_delta_TA = (
        np.dot(
            r_b_eci,
            e_a,
        )
        / r_b
    )

    sin_delta_TA = (
        r_0
        / r_b
    )

    e_cos_TA_b = (
        p / r_b
        - 1
    )

    e_sin_TA = (
        (
            e_cos_TA
            * cos_delta_TA
            - e_cos_TA_b
        )
        / sin_delta_TA
    )

    e = sqrt(
        e_cos_TA ** 2
        + e_sin_TA ** 2
    )

    TA = atan2(
        e_sin_TA,
        e_cos_TA,
    ) % (2 * pi)


    # ========================================================
    # ARGUMENT OF PERIGEE
    # ========================================================

    AOP = (
        argument_of_latitude
        - TA
    ) % (2 * pi)


    # ========================================================
    # SEMI-MAJOR AXIS
    # ========================================================

    a = (
        p
        / (1 - e ** 2)
    )


    # ========================================================
    # ECCENTRIC ANOMALY
    # ========================================================

    E = atan2(
        sqrt(1 - e ** 2)
        * sin(TA),
        cos(TA) + e,
    ) % (2 * pi)


    # ========================================================
    # MEAN ANOMALY
    # ========================================================

    M = (
        E
        - e * sin(E)
    ) % (2 * pi)


    return (
        a,
        e,
        i,
        RAAN,
        AOP,
        M,
    )
