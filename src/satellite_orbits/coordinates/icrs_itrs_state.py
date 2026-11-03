from math import pi

import numpy as np

from satellite_orbits.coordinates.icrs_to_itrs import (
    earth_rotation_matrix,
    nutation_matrix_iau1980,
    polar_motion_matrix,
    precession_matrix_iau1976,
)


EARTH_ROTATION_RATE = (
    1.002737909350795 * 2.0 * pi / 86400.0
)


def icrs_to_itrs_matrix_and_rate(
    mjd_utc: float,
    tt_minus_utc: float,
    ut1_minus_utc: float,
    xp_arcsec: float,
    yp_arcsec: float,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute the celestial-to-terrestrial matrix and its time derivative.

    Inputs:
        mjd_utc:
            Modified Julian Date in UTC [days].
        tt_minus_utc:
            TT minus UTC at the requested epoch [s].
        ut1_minus_utc:
            UT1 minus UTC at the requested epoch [s].
        xp_arcsec, yp_arcsec:
            IERS pole coordinates at that epoch [arcsec].

    Outputs:
        U:
            Dimensionless 3x3 J2000-to-ITRS transformation matrix.
        U_dot:
            Time derivative of U [1/s].

    Conventions:
        Column vectors.
        J2000 mean equator and equinox approximate ICRS.
        The transformation model is inherited from icrs_to_itrs.py.

    Derivative approximation:
        Precession, nutation, and polar motion are held constant.
        Only Earth rotation contributes to U_dot.
        This is the approximation specified in Exercise 5.2.

    Time limitation:
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
        raise ValueError("All time and EOP inputs must be finite.")

    mjd_tt = mjd_utc + tt_minus_utc / 86400.0
    mjd_ut1 = mjd_utc + ut1_minus_utc / 86400.0

    P = precession_matrix_iau1976(mjd_tt)
    N, delta_psi, _, epsilon = nutation_matrix_iau1980(mjd_tt)
    Theta = earth_rotation_matrix(mjd_ut1, delta_psi, epsilon)
    Pi = polar_motion_matrix(xp_arcsec, yp_arcsec)

    rotation_generator = np.array([
        [0.0, 1.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
    ])

    Theta_dot = (
        EARTH_ROTATION_RATE * rotation_generator @ Theta
    )

    U = Pi @ Theta @ N @ P
    U_dot = Pi @ Theta_dot @ N @ P

    return U, U_dot


def _validate_inputs(
    position: np.ndarray,
    velocity: np.ndarray,
    U: np.ndarray,
    U_dot: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Validate vector and matrix shapes and finite values."""
    position = np.asarray(position, dtype=float)
    velocity = np.asarray(velocity, dtype=float)
    U = np.asarray(U, dtype=float)
    U_dot = np.asarray(U_dot, dtype=float)

    if position.shape != (3,) or velocity.shape != (3,):
        raise ValueError("Position and velocity must have shape (3,).")

    if U.shape != (3, 3) or U_dot.shape != (3, 3):
        raise ValueError("U and U_dot must have shape (3, 3).")

    if not all(
        np.all(np.isfinite(value))
        for value in (position, velocity, U, U_dot)
    ):
        raise ValueError("All vector and matrix entries must be finite.")

    return position, velocity, U, U_dot


def icrs_to_itrs_state(
    r_icrs: np.ndarray,
    v_icrs: np.ndarray,
    U: np.ndarray,
    U_dot: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Transform a J2000 position and velocity to ITRS.

    Inputs:
        r_icrs:
            Geocentric J2000 position [km].
        v_icrs:
            Geocentric velocity relative to J2000 axes [km/s].
        U:
            J2000-to-ITRS matrix at the state epoch [dimensionless].
        U_dot:
            Time derivative of that matrix [1/s].

    Outputs:
        r_itrs:
            Geocentric ITRS position [km].
        v_itrs:
            Time derivative of the ITRS coordinates [km/s].

    All inputs must refer to the same physical epoch.
    """
    r_icrs, v_icrs, U, U_dot = _validate_inputs(
        r_icrs, v_icrs, U, U_dot
    )

    r_itrs = U @ r_icrs
    v_itrs = U @ v_icrs + U_dot @ r_icrs

    return r_itrs, v_itrs


def itrs_to_icrs_state(
    r_itrs: np.ndarray,
    v_itrs: np.ndarray,
    U: np.ndarray,
    U_dot: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Transform an ITRS position and velocity to J2000.

    Inputs:
        r_itrs:
            Geocentric ITRS position [km].
        v_itrs:
            Time derivative of the ITRS coordinates [km/s].
        U:
            J2000-to-ITRS matrix at the state epoch [dimensionless].
        U_dot:
            Time derivative of that matrix [1/s].

    Outputs:
        r_icrs:
            Geocentric J2000 position [km].
        v_icrs:
            Geocentric velocity relative to J2000 axes [km/s].

    All inputs must refer to the same physical epoch.

    Linear solves invert the implemented matrix directly.
    This also accommodates the linearized polar-motion matrix,
    which is not exactly orthogonal.
    """
    r_itrs, v_itrs, U, U_dot = _validate_inputs(
        r_itrs, v_itrs, U, U_dot
    )

    r_icrs = np.linalg.solve(U, r_itrs)
    v_icrs = np.linalg.solve(
        U,
        v_itrs - U_dot @ r_icrs,
    )

    return r_icrs, v_icrs
