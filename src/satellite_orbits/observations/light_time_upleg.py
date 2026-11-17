from typing import Callable

import numpy as np


def light_time_upleg(
    station_position_eci: Callable[[float], np.ndarray],
    satellite_position_eci: np.ndarray,
    transponding_time: float,
    c_light: float,
    tolerance: float = 1.0e-9,
    max_iterations: int = 10,
) -> tuple[float, list[float]]:
    """
    Solve the upleg light-time equation by fixed-point iteration.

    The signal is transponded by the satellite at transponding_time
    and was emitted by the station at transponding_time - travel_time.
    The satellite position is evaluated at the transponding time, so
    it is passed as a fixed vector, while the station moves during
    the signal travel because of the Earth rotation.

    Parameters
    ----------
    station_position_eci : Callable
        Function returning the inertial station position [km] for a
        time expressed in the same scale as transponding_time.

    satellite_position_eci : np.ndarray
        Satellite position in the inertial frame at the transponding
        time [km].

    transponding_time : float
        Time at which the satellite transponds the signal [s].

    c_light : float
        Speed of light [km/s].

    tolerance : float
        Convergence threshold on the signal path length [km].

    max_iterations : int
        Maximum number of fixed-point iterations.

    Returns
    -------
    tuple
        travel_time : Upleg signal travel time [s].
        distances : Signal path length [km] after each iteration,
            starting with the geometric distance obtained for a
            zero initial travel time.
    """

    travel_time = 0.0
    distances = []

    for _ in range(max_iterations):

        r_station = station_position_eci(
            transponding_time - travel_time
        )

        distance = np.linalg.norm(
            satellite_position_eci - r_station
        )

        distances.append(float(distance))

        travel_time = float(distance) / c_light

        if len(distances) > 1:
            if abs(distances[-1] - distances[-2]) < tolerance:
                break

    return (
        travel_time,
        distances,
    )
