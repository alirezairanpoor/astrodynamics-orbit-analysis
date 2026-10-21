from math import cos, sin

import numpy as np


def ecef_to_enz(
    r_ecef: np.ndarray,
    station_longitude: float,
    station_latitude: float,
    earth_radius: float,
) -> np.ndarray:
    """
    Transform a satellite position from ECEF
    to the local East-North-Zenith frame
    of a ground station.

    Parameters
    ----------
    r_ecef : np.ndarray
        Satellite position vector in ECEF [km].

    station_longitude : float
        Ground-station longitude [rad].

    station_latitude : float
        Ground-station latitude [rad].

    earth_radius : float
        Earth radius [km].

    Returns
    -------
    np.ndarray
        Station-to-satellite vector in ENZ [km].
    """

    longitude = station_longitude
    latitude = station_latitude

    # Ground-station position in ECEF
    r_station_ecef = earth_radius * np.array([
        cos(latitude) * cos(longitude),
        cos(latitude) * sin(longitude),
        sin(latitude),
    ])

    # Vector from ground station to satellite
    s_ecef = (
        r_ecef
        - r_station_ecef
    )

    # Local ENZ unit vectors expressed in ECEF
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

    s_enz = rotation_matrix @ s_ecef

    return s_enz
