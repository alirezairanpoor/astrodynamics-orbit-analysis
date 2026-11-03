from math import cos, pi, sin

import erfa
import numpy as np


ARCSEC_TO_RAD = pi / (180.0 * 3600.0)
SECONDS_PER_DAY = 86400.0
MJD_J2000 = 51544.5


def precession_matrix_iau1976(mjd_tt: float) -> np.ndarray:
    """
    Input:
        mjd_tt: Modified Julian Date in TT [days].

    Output:
        P: Dimensionless 3x3 matrix mapping J2000 mean-equatorial
           column vectors to mean-equatorial coordinates of date.
    """
    T = (mjd_tt - MJD_J2000) / 36525.0

    zeta = (
        2306.2181 * T
        + 0.30188 * T**2
        + 0.017998 * T**3
    ) * ARCSEC_TO_RAD

    theta = (
        2004.3109 * T
        - 0.42665 * T**2
        - 0.041833 * T**3
    ) * ARCSEC_TO_RAD

    z = zeta + (
        0.79280 * T**2
        + 0.000205 * T**3
    ) * ARCSEC_TO_RAD

    cz, sz = cos(z), sin(z)
    ct, st = cos(theta), sin(theta)
    cze, sze = cos(zeta), sin(zeta)

    return np.array([
        [
            -sz * sze + cz * ct * cze,
            -sz * cze - cz * ct * sze,
            -cz * st,
        ],
        [
            cz * sze + sz * ct * cze,
            cz * cze - sz * ct * sze,
            -sz * st,
        ],
        [st * cze, -st * sze, ct],
    ])


def nutation_matrix_iau1980(
    mjd_tt: float,
) -> tuple[np.ndarray, float, float, float]:
    """
    Input:
        mjd_tt: Modified Julian Date in TT [days].

    Outputs:
        N: Dimensionless 3x3 mean-of-date to true-of-date matrix.
        delta_psi: Nutation in longitude [rad].
        delta_epsilon: Nutation in obliquity [rad].
        epsilon: Mean obliquity [rad].

    The complete IAU 1980 series is evaluated by ERFA.
    Observed corrections to the nutation angles are not included.
    """
    T = (mjd_tt - MJD_J2000) / 36525.0

    delta_psi, delta_epsilon = erfa.nut80(2400000.5, mjd_tt)
    delta_psi = float(delta_psi)
    delta_epsilon = float(delta_epsilon)

    epsilon = (
        23.43929111 * 3600.0
        - 46.8150 * T
        - 0.00059 * T**2
        + 0.001813 * T**3
    ) * ARCSEC_TO_RAD

    epsilon_true = epsilon + delta_epsilon

    ce, se = cos(epsilon), sin(epsilon)
    cet, set_ = cos(epsilon_true), sin(epsilon_true)
    cp, sp = cos(delta_psi), sin(delta_psi)

    N = np.array([
        [cp, -ce * sp, -se * sp],
        [
            cet * sp,
            ce * cet * cp + se * set_,
            se * cet * cp - ce * set_,
        ],
        [
            set_ * sp,
            ce * set_ * cp - se * cet,
            se * set_ * cp + ce * cet,
        ],
    ])

    return N, delta_psi, delta_epsilon, epsilon


def earth_rotation_matrix(
    mjd_ut1: float,
    delta_psi: float,
    epsilon: float,
) -> np.ndarray:
    """
    Inputs:
        mjd_ut1: Modified Julian Date in UT1 [days].
        delta_psi: Nutation in longitude at the same epoch [rad].
        epsilon: Mean obliquity at the same epoch [rad].

    Output:
        Theta: Dimensionless 3x3 rotation from true-of-date
               coordinates to the rotating Greenwich frame.

    Uses GMST82 and the first-order equation of the equinoxes.
    """
    gmst = float(erfa.gmst82(2400000.5, mjd_ut1))
    gast = (gmst + delta_psi * cos(epsilon)) % (2.0 * pi)

    c, s = cos(gast), sin(gast)

    return np.array([
        [c, s, 0.0],
        [-s, c, 0.0],
        [0.0, 0.0, 1.0],
    ])


def polar_motion_matrix(
    xp_arcsec: float,
    yp_arcsec: float,
) -> np.ndarray:
    """
    Inputs:
        xp_arcsec: IERS pole coordinate xp [arcsec].
        yp_arcsec: IERS pole coordinate yp [arcsec].

    Output:
        Pi: Dimensionless 3x3 first-order polar-motion matrix
            mapping the rotating Greenwich frame to ITRS.
    """
    xp = xp_arcsec * ARCSEC_TO_RAD
    yp = yp_arcsec * ARCSEC_TO_RAD

    return np.array([
        [1.0, 0.0, xp],
        [0.0, 1.0, -yp],
        [-xp, yp, 1.0],
    ])


def icrs_to_itrs_matrix(
    mjd_utc: float,
    tt_minus_utc: float,
    ut1_minus_utc: float,
    xp_arcsec: float,
    yp_arcsec: float,
) -> np.ndarray:
    """
    Inputs:
        mjd_utc: Modified Julian Date in UTC [days].
        tt_minus_utc: TT minus UTC at the requested epoch [s].
        ut1_minus_utc: UT1 minus UTC at the requested epoch [s].
        xp_arcsec: IERS pole coordinate xp at that epoch [arcsec].
        yp_arcsec: IERS pole coordinate yp at that epoch [arcsec].

    Output:
        U: Dimensionless 3x3 coordinate transformation matrix.
           For a position column vector:
               r_itrs = U @ r_icrs
           The position units are preserved.

    Model:
        Textbook J2000 mean-equator/equinox approximation to ICRS.
        IAU 1976 precession and IAU 1980 nutation.
        GMST82 with the first-order equation of the equinoxes.
        First-order polar motion.

    Limitations:
        No frame-bias matrix or observed nutation corrections.
        EOP values must be supplied for the requested epoch.
        This matrix alone is not a velocity transformation.
        Leap-second instants require dedicated time handling.
    """
    inputs = np.array([
        mjd_utc,
        tt_minus_utc,
        ut1_minus_utc,
        xp_arcsec,
        yp_arcsec,
    ], dtype=float)

    if not np.all(np.isfinite(inputs)):
        raise ValueError("All inputs must be finite.")

    mjd_tt = mjd_utc + tt_minus_utc / SECONDS_PER_DAY
    mjd_ut1 = mjd_utc + ut1_minus_utc / SECONDS_PER_DAY

    P = precession_matrix_iau1976(mjd_tt)
    N, delta_psi, _, epsilon = nutation_matrix_iau1980(mjd_tt)
    Theta = earth_rotation_matrix(mjd_ut1, delta_psi, epsilon)
    Pi = polar_motion_matrix(xp_arcsec, yp_arcsec)

    return Pi @ Theta @ N @ P
