from math import cos, sin

import numpy as np


def enz_to_ecef(
    s_enz: np.ndarray,
    station_longitude: float,
    station_latitude: float,
) -> np.ndarray:
    """
    Transform a station-to-satellite vector
    from the local ENZ frame to ECEF.

    Parameters
    ----------
    s_enz : np.ndarray
        Station-to-satellite vector in ENZ [km].

    station_longitude : float
        Ground-station longitude [rad].

    station_latitude : float
        Ground-station latitude [rad].

    Returns
    -------
    np.ndarray
        Station-to-satellite vector in ECEF [km].
    """

    longitude = station_longitude
    latitude = station_latitude

    east = np.array([
        -sin(longitude),
        cos(longitude),
        0.0,
    ])

    north = np.array([
        -sin(latitude) * cos(longitude),
        -sin(latitude) * sin(longitude),
        cos(latitude),
    ])

    zenith = np.array([
        cos(latitude) * cos(longitude),
        cos(latitude) * sin(longitude),
        sin(latitude),
    ])

    rotation_matrix = np.array([
        east,
        north,
        zenith,
    ])

    s_ecef = rotation_matrix.T @ s_enz

    return s_ecef
