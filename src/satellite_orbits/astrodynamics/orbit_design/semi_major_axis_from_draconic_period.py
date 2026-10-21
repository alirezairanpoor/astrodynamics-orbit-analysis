from math import pi


def semi_major_axis_from_draconic_period(
    TN: float,
    delta_n: float,
    AOP_rate: float,
    mu: float,
) -> float:
    """
    Compute a refined semi-major axis from the
    draconic period including secular J2 corrections.

    Parameters
    ----------
    TN : float
        Draconic period [s].

    delta_n : float
        J2 correction to mean motion [rad/s].

    AOP_rate : float
        Secular argument-of-perigee rate [rad/s].

    mu : float
        Gravitational parameter [km^3/s^2].

    Returns
    -------
    float
        Refined semi-major axis a [km].
    """

    n0 = (
        2 * pi / TN
        - delta_n
        - AOP_rate
    )

    a = (
        mu / n0 ** 2
    ) ** (1 / 3)

    return a
