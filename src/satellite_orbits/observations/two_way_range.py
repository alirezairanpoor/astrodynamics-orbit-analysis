from typing import Callable

import numpy as np

from .light_time_downleg import light_time_downleg
from .light_time_upleg import light_time_upleg


def two_way_range(
    satellite_position: Callable[[float], np.ndarray],
    station_position: Callable[[float], np.ndarray],
    reception_time: float,
    c_light: float,
    tolerance: float = 1.0e-9,
    max_iterations: int = 10,
) -> tuple[float, float, float]:
    """
    Model a two-way range measurement taken at the ground station.

    The signal is emitted by the station, transponded by the
    satellite, and received again at the same station at
    reception_time. The downleg is solved first, which fixes the
    transponding epoch, and the upleg is solved from there. The
    two-way range is the mean of the two one-way path lengths.

    Parameters
    ----------
    satellite_position : Callable
        Function returning the inertial satellite position [km].

    station_position : Callable
        Function returning the inertial station position [km].

    reception_time : float
        Time of signal reception at the station [s].

    c_light : float
        Speed of light [km/s].

    tolerance : float
        Convergence threshold on the signal path length [km].

    max_iterations : int
        Maximum number of fixed-point iterations per leg.

    Returns
    -------
    tuple
        range_two_way : Mean of the uplink and downlink path
            lengths [km].
        travel_time_down : Downleg signal travel time [s].
        travel_time_up : Upleg signal travel time [s].
    """

    travel_time_down, _ = light_time_downleg(
        satellite_position,
        station_position(reception_time),
        reception_time,
        c_light,
        tolerance,
        max_iterations,
    )

    transponding_time = reception_time - travel_time_down

    travel_time_up, _ = light_time_upleg(
        station_position,
        satellite_position(transponding_time),
        transponding_time,
        c_light,
        tolerance,
        max_iterations,
    )

    range_two_way = (
        0.5
        * c_light
        * (travel_time_down + travel_time_up)
    )

    return (
        range_two_way,
        travel_time_down,
        travel_time_up,
    )