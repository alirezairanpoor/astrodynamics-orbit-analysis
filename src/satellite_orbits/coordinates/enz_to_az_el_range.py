from math import atan2, sqrt, pi

import numpy as np


def enz_to_az_el_range(
    s_enz: np.ndarray,
) -> tuple[float, float, float]:
    """
    Convert a local ENZ vector
    to azimuth, elevation, and range.
    """

    s_east = s_enz[0]
    s_north = s_enz[1]
    s_zenith = s_enz[2]

    azimuth = atan2(
        s_east,
        s_north,
    ) % (2 * pi)

    horizontal_distance = sqrt(
        s_east ** 2
        + s_north ** 2
    )

    elevation = atan2(
        s_zenith,
        horizontal_distance,
    )

    range_km = np.linalg.norm(
        s_enz
    )

    return (
        azimuth,
        elevation,
        range_km,
    )
