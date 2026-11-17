from math import cos, sin

import numpy as np


def ecef_to_enz_rotation(
    s_ecef: np.ndarray,
    station_longitude: float,
    station_latitude: float,
) -> np.ndarray:
    """
    Rotate a vector from ECEF to the local East-North-Zenith frame.

    This is a pure rotation with no translation: the input is
    already a station-to-satellite vector. It is the exact inverse
    of enz_to_ecef.

    The rotation projects the input onto the east, north and zenith
    unit vectors of the station, which are defined by the geodetic
    latitude and longitude.

    Keeping the translation outside the function lets the caller
    build the station position with the geometry the problem needs.
    ecef_to_enz builds it on a sphere, which is adequate for many
    purposes but displaces a mid-latitude station by about 20 km
    relative to the reference ellipsoid.

    Source: equations 6.26 to 6.28.

    Parameters
    ----------
    s_ecef : np.ndarray
        Station-to-satellite vector in ECEF [km].

    station_longitude : float
        Ground-station geodetic longitude [rad].

    station_latitude : float
        Ground-station geodetic latitude [rad].

    Returns
    -------
    np.ndarray
        Station-to-satellite vector in ENZ [km].
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

    s_enz = rotation_matrix @ s_ecef

    return s_enz
