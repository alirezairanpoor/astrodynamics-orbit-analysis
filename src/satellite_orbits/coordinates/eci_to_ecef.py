from math import cos, sin

import numpy as np


def eci_to_ecef(
    r_eci: np.ndarray,
    greenwich_hour_angle: float,
) -> np.ndarray:
    """
    Transform a position vector from ECI to ECEF.

    Parameters
    ----------
    r_eci : np.ndarray
        Position vector in ECI [km].

    greenwich_hour_angle : float
        Greenwich Hour Angle [rad].

    Returns
    -------
    np.ndarray
        Position vector in ECEF [km].
    """

    theta = greenwich_hour_angle

    rotation_matrix = np.array([
        [cos(theta),  sin(theta), 0.0],
        [-sin(theta), cos(theta), 0.0],
        [0.0,         0.0,        1.0],
    ])

    r_ecef = rotation_matrix @ r_eci

    return r_ecef
