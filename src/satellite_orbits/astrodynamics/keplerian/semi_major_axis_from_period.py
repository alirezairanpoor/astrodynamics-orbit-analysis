from math import pi


def semi_major_axis_from_period(
    period: float,
    mu: float,
) -> float:
    """
    Compute the semi-major axis from the
    Keplerian orbital period.

    Parameters
    ----------
    period : float
        Orbital period [s].

    mu : float
        Gravitational parameter [km^3/s^2].

    Returns
    -------
    float
        Semi-major axis a [km].
    """

    a = (
        mu
        * (period / (2 * pi)) ** 2
    ) ** (1 / 3)

    return a
