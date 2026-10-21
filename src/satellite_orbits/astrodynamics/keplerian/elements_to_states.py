from math import cos, sin, sqrt

import numpy as np

from .kepler import solve_kepler_newton


def elements_to_states(
    a: float,
    e: float,
    i: float,
    RAAN: float,
    AOP: float,
    M: float,
    mu: float,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Convert classical Keplerian orbital elements
    to ECI position and velocity vectors.
    """

    E, _, _ = solve_kepler_newton(
        M,
        e,
    )

    P = np.array([
        cos(AOP) * cos(RAAN)
        - sin(AOP) * cos(i) * sin(RAAN),

        cos(AOP) * sin(RAAN)
        + sin(AOP) * cos(i) * cos(RAAN),

        sin(AOP) * sin(i),
    ])

    Q = np.array([
        -sin(AOP) * cos(RAAN)
        - cos(AOP) * cos(i) * sin(RAAN),

        -sin(AOP) * sin(RAAN)
        + cos(AOP) * cos(i) * cos(RAAN),

        cos(AOP) * sin(i),
    ])

    r_mag = a * (
        1 - e * cos(E)
    )

    r_eci = (
        a * (cos(E) - e) * P
        + a * sqrt(1 - e ** 2) * sin(E) * Q
    )

    v_eci = (
        sqrt(mu * a) / r_mag
        * (
            -sin(E) * P
            + sqrt(1 - e ** 2) * cos(E) * Q
        )
    )

    return r_eci, v_eci
