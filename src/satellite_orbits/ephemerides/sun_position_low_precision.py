from math import sin, cos, radians

import numpy as np


def sun_position_low_precision(mjd_tt: float) -> np.ndarray:
    """
    Compute the approximate geocentric position of the Sun.

    Based on Montenbruck & Gill, equations 3.43–3.46.

    Input:
        mjd_tt: Modified Julian Date expressed in TT.

    Output:
        Earth-to-Sun position in EME2000 [km].

    Intended for low-precision calculations near J2000.
    """

    # --------------------------------------------------------
    # JULIAN CENTURIES SINCE J2000
    # --------------------------------------------------------

    T = (mjd_tt - 51544.5) / 36525.0

    # --------------------------------------------------------
    # MEAN ANOMALY
    # --------------------------------------------------------

    M = radians(
        (357.5256 + 35999.049 * T) % 360.0
    )

    # --------------------------------------------------------
    # ECLIPTIC LONGITUDE
    # --------------------------------------------------------

    arcsec_to_rad = radians(1.0 / 3600.0)

    longitude = (
        radians(282.9400)
        + M
        + arcsec_to_rad * (
            6892.0 * sin(M)
            + 72.0 * sin(2.0 * M)
        )
    )

    # --------------------------------------------------------
    # EARTH–SUN DISTANCE
    # --------------------------------------------------------

    distance = (
        149.619
        - 2.499 * cos(M)
        - 0.021 * cos(2.0 * M)
    ) * 1e6  # km

    # --------------------------------------------------------
    # EQUATORIAL J2000 CARTESIAN POSITION
    # --------------------------------------------------------

    # Ecliptic latitude is approximated as zero.
    epsilon = radians(23.43929111)

    r_sun_eci = distance * np.array([
        cos(longitude),
        sin(longitude) * cos(epsilon),
        sin(longitude) * sin(epsilon),
    ])

    return r_sun_eci
