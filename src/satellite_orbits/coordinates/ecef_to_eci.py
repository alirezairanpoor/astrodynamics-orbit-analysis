from math import cos, sin

import numpy as np


def ecef_to_eci(
    r_ecef: np.ndarray,
    greenwich_hour_angle: float,
) -> np.ndarray:
    """
    Transform a position vector from ECEF to ECI.

    Parameters
    ----------
    r_ecef : np.ndarray
        Position vector in ECEF [km].

    greenwich_hour_angle : float
        Greenwich Hour Angle [rad].

    Returns
    -------
    np.ndarray
        Position vector in ECI [km].
    """

    theta = greenwich_hour_angle

    rotation_matrix = np.array([
        [cos(theta),  sin(theta), 0.0],
        [-sin(theta), cos(theta), 0.0],
        [0.0,         0.0,        1.0],
    ])

    r_eci = (
        rotation_matrix.T
        @ r_ecef
    )

    return r_eci
