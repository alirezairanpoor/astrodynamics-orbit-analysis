from typing import Callable

import numpy as np


def light_time_downleg(
    satellite_position: Callable[[float], np.ndarray],
    station_position_eci: np.ndarray,
    reception_time: float,
    c_light: float,
    tolerance: float = 1.0e-9,
    max_iterations: int = 10,
) -> tuple[float, list[float]]:
    """
    Solve the downleg light-time equation by fixed-point iteration.

    The signal is received at the station at reception_time and was
    emitted by the satellite at reception_time - travel_time. The
    station position is evaluated at the reception time, so it is
    passed as a fixed vector.

    Parameters
    ----------
    satellite_position : Callable
        Function returning the inertial satellite position [km]
        for a time expressed in the same scale as reception_time.

    station_position_eci : np.ndarray
        Station position in the inertial frame at the reception
        time [km].

    reception_time : float
        Time of signal reception at the station [s].

    c_light : float
        Speed of light [km/s].

    tolerance : float
        Convergence threshold on the signal path length [km].

    max_iterations : int
        Maximum number of fixed-point iterations.

    Returns
    -------
    tuple
        travel_time : Downleg signal travel time [s].
        distances : Signal path length [km] after each iteration,
            starting with the geometric distance obtained for a
            zero initial travel time.
    """

    travel_time = 0.0
    distances = []

    for _ in range(max_iterations):

        r_satellite = satellite_position(
            reception_time - travel_time
        )

        distance = np.linalg.norm(
            r_satellite - station_position_eci
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
