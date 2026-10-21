from math import sqrt, pi


def propagate_mean_anomaly(
    a: float,
    M0: float,
    elapsed_time: float,
    mu: float,
) -> float:
    """
    Propagate the mean anomaly in an unperturbed Keplerian orbit.

    Parameters
    ----------
    a : float
        Semi-major axis [km].

    M0 : float
        Mean anomaly at the reference epoch [rad].

    elapsed_time : float
        Time elapsed since the reference epoch [s].

    mu : float
        Gravitational parameter [km^3/s^2].

    Returns
    -------
    float
        Propagated mean anomaly M [rad].
    """

    mean_motion = sqrt(
        mu / a ** 3
    )

    M = (
        M0
        + mean_motion * elapsed_time
    ) % (2 * pi)

    return M
