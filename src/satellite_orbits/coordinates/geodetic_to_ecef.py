from math import cos, sin, sqrt

import numpy as np


def geodetic_to_ecef(
    longitude: float,
    latitude: float,
    height: float,
    earth_radius: float,
    flattening: float,
) -> np.ndarray:
    """
    Convert geodetic coordinates to an ECEF position vector.

    Inverse of ecef_to_geodetic. Geodetic latitude is measured
    from the equatorial plane to the ellipsoid normal, not to the
    geocentric radius vector.

    Parameters
    ----------
    longitude : float
        Geodetic longitude, positive east [rad].

    latitude : float
        Geodetic latitude, positive north [rad].

    height : float
        Height above the reference ellipsoid [km].

    earth_radius : float
        Equatorial Earth radius [km].

    flattening : float
        Earth flattening [-].

    Returns
    -------
    np.ndarray
        Position vector in ECEF [km].
    """

    eccentricity_squared = (
        flattening
        * (2 - flattening)
    )

    # Radius of curvature in the prime vertical.
    N = (
        earth_radius
        / sqrt(
            1
            - eccentricity_squared
            * sin(latitude) ** 2
        )
    )

    # The height is added along the ellipsoid normal, while the
    # (1 - e^2) factor accounts for the normal not passing
    # through the geocenter.
    x = (
        (N + height)
        * cos(latitude)
        * cos(longitude)
    )

    y = (
        (N + height)
        * cos(latitude)
        * sin(longitude)
    )

    z = (
        (
            (1 - eccentricity_squared) * N
            + height
        )
        * sin(latitude)
    )

    r_ecef = np.array([
        x,
        y,
        z,
    ])

    return r_ecef
