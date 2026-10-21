from math import sin, cos, radians

import numpy as np


def moon_position_low_precision(mjd_tt: float) -> np.ndarray:
    """
    Compute the approximate geocentric position of the Moon.

    Based on Montenbruck & Gill, equations 3.47–3.51.

    Input
    -----
    mjd_tt : Modified Julian Date expressed in TT

    Output
    ------
    r_moon_eci : Moon position in EME2000 [km]
    """

    # --------------------------------------------------------
    # TIME: JULIAN CENTURIES SINCE J2000
    # --------------------------------------------------------

    T = (mjd_tt - 51544.5) / 36525.0

    # --------------------------------------------------------
    # FUNDAMENTAL ANGLES — EQUATION 3.47
    # --------------------------------------------------------

    L0 = radians(
        (
            218.31617
            + 481267.88088 * T
            - 1.3972 * T
        ) % 360.0
    )

    M_moon = radians(
        (134.96292 + 477198.86753 * T) % 360.0
    )

    M_sun = radians(
        (357.52543 + 35999.04944 * T) % 360.0
    )

    F = radians(
        (93.27283 + 483202.01873 * T) % 360.0
    )

    D = radians(
        (297.85027 + 445267.11135 * T) % 360.0
    )

    arcsec_to_rad = radians(1.0 / 3600.0)

    # --------------------------------------------------------
    # ECLIPTIC LONGITUDE — EQUATION 3.48
    # --------------------------------------------------------

    delta_longitude = arcsec_to_rad * (
        22640 * sin(M_moon)
        + 769 * sin(2 * M_moon)
        - 4586 * sin(M_moon - 2 * D)
        + 2370 * sin(2 * D)
        - 668 * sin(M_sun)
        - 412 * sin(2 * F)
        - 212 * sin(2 * M_moon - 2 * D)
        - 206 * sin(M_moon + M_sun - 2 * D)
        + 192 * sin(M_moon + 2 * D)
        - 165 * sin(M_sun - 2 * D)
        + 148 * sin(M_moon - M_sun)
        - 125 * sin(D)
        - 110 * sin(M_moon + M_sun)
        - 55 * sin(2 * F - 2 * D)
    )

    longitude = L0 + delta_longitude

    # --------------------------------------------------------
    # ECLIPTIC LATITUDE — EQUATION 3.49
    # --------------------------------------------------------

    latitude = arcsec_to_rad * (
        18520 * sin(
            F
            + delta_longitude
            + arcsec_to_rad * (
                412 * sin(2 * F)
                + 541 * sin(M_sun)
            )
        )
        - 526 * sin(F - 2 * D)
        + 44 * sin(M_moon + F - 2 * D)
        - 31 * sin(-M_moon + F - 2 * D)
        - 25 * sin(-2 * M_moon + F)
        - 23 * sin(M_sun + F - 2 * D)
        + 21 * sin(-M_moon + F)
        + 11 * sin(-M_sun + F - 2 * D)
    )

    # --------------------------------------------------------
    # GEOCENTRIC DISTANCE — EQUATION 3.50
    # --------------------------------------------------------

    distance = (
        385000
        - 20905 * cos(M_moon)
        - 3699 * cos(2 * D - M_moon)
        - 2956 * cos(2 * D)
        - 570 * cos(2 * M_moon)
        + 246 * cos(2 * M_moon - 2 * D)
        - 205 * cos(M_sun - 2 * D)
        - 171 * cos(M_moon + 2 * D)
        - 152 * cos(M_moon + M_sun - 2 * D)
    )  # km

    # --------------------------------------------------------
    # ECLIPTIC CARTESIAN POSITION
    # --------------------------------------------------------

    x_ecliptic = distance * cos(longitude) * cos(latitude)
    y_ecliptic = distance * sin(longitude) * cos(latitude)
    z_ecliptic = distance * sin(latitude)

    # --------------------------------------------------------
    # ECLIPTIC J2000 -> EQUATORIAL J2000 — EQUATION 3.51
    # --------------------------------------------------------

    epsilon = radians(23.43929111)

    r_moon_eci = np.array([
        x_ecliptic,
        (
            cos(epsilon) * y_ecliptic
            - sin(epsilon) * z_ecliptic
        ),
        (
            sin(epsilon) * y_ecliptic
            + cos(epsilon) * z_ecliptic
        ),
    ])

    return r_moon_eci
