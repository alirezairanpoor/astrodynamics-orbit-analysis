from math import sqrt, asin, sin

import numpy as np


def gauss_sector_triangle_ratio(
    r_a_eci: np.ndarray,
    r_b_eci: np.ndarray,
    delta_t: float,
    mu: float,
    tolerance: float = 1e-12,
) -> float:
    """
    Compute Gauss's sector-to-triangle ratio eta
    from two ECI position vectors and their time interval.

    Parameters
    ----------
    r_a_eci : np.ndarray
        First ECI position vector [km].

    r_b_eci : np.ndarray
        Second ECI position vector [km].

    delta_t : float
        Time from the first to the second observation [s].

    mu : float
        Gravitational parameter [km^3/s^2].

    tolerance : float
        Convergence tolerance.

    Returns
    -------
    float
        Sector-to-triangle ratio eta.
    """

    r_a = np.linalg.norm(r_a_eci)
    r_b = np.linalg.norm(r_b_eci)

    tau = (
        sqrt(mu)
        * delta_t
    )

    geometry_term = sqrt(
        2
        * (
            r_a * r_b
            + np.dot(r_a_eci, r_b_eci)
        )
    )

    m = (
        tau ** 2
        / geometry_term ** 3
    )

    l = (
        (r_a + r_b)
        / (2 * geometry_term)
        - 0.5
    )

    # Hansen approximation
    eta_0 = (
        12 / 22
        + 10 / 22
        * sqrt(
            1
            + 44 / 9
            * m
            / (l + 5 / 6)
        )
    )

    def W(w: float) -> float:

        g = (
            2
            * asin(
                sqrt(w)
            )
        )

        return (
            2 * g
            - sin(2 * g)
        ) / (
            sin(g) ** 3
        )

    def f(eta: float) -> float:

        w = (
            m / eta ** 2
            - l
        )

        return (
            1
            - eta
            + m / eta ** 2
            * W(w)
        )

    # Secant starting values from the book
    eta_previous = (
        eta_0 + 0.1
    )

    eta = eta_0

    while True:

        f_previous = f(
            eta_previous
        )

        f_current = f(
            eta
        )

        eta_new = (
            eta
            - f_current
            * (
                eta
                - eta_previous
            )
            / (
                f_current
                - f_previous
            )
        )

        if abs(
            eta_new - eta
        ) < tolerance:
            return eta_new

        eta_previous = eta
        eta = eta_new
