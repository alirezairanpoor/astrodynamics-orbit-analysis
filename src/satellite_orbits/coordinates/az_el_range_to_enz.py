from math import sin, cos

import numpy as np


def az_el_range_to_enz(
    azimuth: float,
    elevation: float,
    range_km: float,
) -> np.ndarray:
    """
    Convert azimuth, elevation, and range
    to a local ENZ vector.
    """

    s_east = (
        range_km
        * cos(elevation)
        * sin(azimuth)
    )

    s_north = (
        range_km
        * cos(elevation)
        * cos(azimuth)
    )

    s_zenith = (
        range_km
        * sin(elevation)
    )

    s_enz = np.array([
        s_east,
        s_north,
        s_zenith,
    ])

    return s_enz
