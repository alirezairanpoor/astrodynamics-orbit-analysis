from math import cos, sin
from typing import Callable

import numpy as np


def light_time_earth_fixed(
    satellite_position_ecef: Callable[[float], np.ndarray],
    receiver_position_ecef: np.ndarray,
    reception_time: float,
    c_light: float,
    omega_earth: float,
    tolerance: float = 1.0e-9,
    max_iterations: int = 10,
) -> tuple[float, float]:
    """
    Solve the light-time equation directly in the Earth-fixed frame.

    Tabulated GPS ephemerides are given in an Earth-fixed frame, so
    transforming to an inertial frame would require precession,
    nutation and UT1 data that a pseudorange model at the metre
    level does not need.

    The difficulty is that the satellite position at the emission
    time refers to the Earth-fixed frame of that earlier epoch,
    while the receiver position refers to the frame at reception.
    Subtracting them directly mixes two frames. Using the
    approximation

        U(t) = Rz(omega * tau) U(t - tau)

    which ignores precession and nutation over the signal travel
    time, the emission-epoch satellite position is re-expressed in
    the reception-epoch frame, giving equation 6.77 of the
    textbook.

    Rz here is the coordinate rotation, which expresses a fixed
    vector in a frame that has itself rotated. It is the transpose
    of the matrix that rotates a vector within a fixed frame, and
    it matches the convention of eci_to_ecef in this package.

    Parameters
    ----------
    satellite_position_ecef : Callable
        Function returning the Earth-fixed satellite position for a
        time expressed in the same scale as reception_time.

    receiver_position_ecef : np.ndarray
        Receiver position in the Earth-fixed frame at the reception
        time.

    reception_time : float
        Time of signal reception at the receiver [s].

    c_light : float
        Speed of light, in the same length unit per second as the
        position vectors.

    omega_earth : float
        Earth rotation rate [rad/s].

    tolerance : float
        Convergence threshold on the signal path length, in the
        same length unit as the position vectors.

    max_iterations : int
        Maximum number of fixed-point iterations.

    Returns
    -------
    tuple
        travel_time : Signal travel time [s].
        distance : Signal path length, in the length unit of the
            position vectors.
    """

    travel_time = 0.0
    distance = 0.0
    previous_distance = None

    for _ in range(max_iterations):

        r_satellite = satellite_position_ecef(
            reception_time - travel_time
        )

        # Re-express the emission-epoch position in the Earth-fixed
        # frame of the reception epoch.
        rotation_angle = omega_earth * travel_time

        rotation_matrix = np.array([
            [
                cos(rotation_angle),
                sin(rotation_angle),
                0.0,
            ],
            [
                -sin(rotation_angle),
                cos(rotation_angle),
                0.0,
            ],
            [
                0.0,
                0.0,
                1.0,
            ],
        ])

        r_satellite_rotated = (
            rotation_matrix
            @ r_satellite
        )

        distance = float(
            np.linalg.norm(
                r_satellite_rotated
                - receiver_position_ecef
            )
        )

        travel_time = distance / c_light

        if previous_distance is not None:
            if abs(distance - previous_distance) < tolerance:
                break

        previous_distance = distance

    return (
        travel_time,
        distance,
    )
