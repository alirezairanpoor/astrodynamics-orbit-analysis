from math import atan2, sqrt

import numpy as np


def ecef_to_geodetic(
    r_ecef: np.ndarray,
    earth_radius: float,
    flattening: float,
    tolerance: float = 1e-12,
) -> tuple[float, float, float]:
    """
    Convert an ECEF position vector to
    geodetic longitude, latitude, and height.

    Parameters
    ----------
    r_ecef : np.ndarray
        Position vector in ECEF [km].

    earth_radius : float
        Equatorial Earth radius [km].

    flattening : float
        Earth flattening [-].

    tolerance : float
        Iteration convergence tolerance [km].

    Returns
    -------
    tuple
        longitude : Geodetic longitude [rad].
        latitude  : Geodetic latitude [rad].
        height    : Height above reference ellipsoid [km].
    """

    x = r_ecef[0]
    y = r_ecef[1]
    z = r_ecef[2]

    horizontal_distance = sqrt(
        x ** 2
        + y ** 2
    )

    eccentricity_squared = (
        flattening
        * (2 - flattening)
    )

    delta_z = (
        eccentricity_squared
        * z
    )

    while True:

        sin_latitude = (
            (z + delta_z)
            / sqrt(
                horizontal_distance ** 2
                + (z + delta_z) ** 2
            )
        )

        N = (
            earth_radius
            / sqrt(
                1
                - eccentricity_squared
                * sin_latitude ** 2
            )
        )

        new_delta_z = (
            N
            * eccentricity_squared
            * sin_latitude
        )

        if abs(
            new_delta_z - delta_z
        ) < tolerance:
            delta_z = new_delta_z
            break

        delta_z = new_delta_z

    longitude = atan2(
        y,
        x,
    )

    latitude = atan2(
        z + delta_z,
        horizontal_distance,
    )

    height = (
        sqrt(
            horizontal_distance ** 2
            + (z + delta_z) ** 2
        )
        - N
    )

    return (
        longitude,
        latitude,
        height,
    )
