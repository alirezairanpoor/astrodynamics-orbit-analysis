from math import sqrt


def circular_velocity(mu: float, radius: float) -> float:
    """
    Compute the circular orbital velocity.

    Parameters
    ----------
    mu : float
        Gravitational parameter [km^3/s^2].
    radius : float
        Orbital radius from the center of the central body [km].

    Returns
    -------
    float
        Circular orbital velocity [km/s].
    """
    return sqrt(mu / radius)
