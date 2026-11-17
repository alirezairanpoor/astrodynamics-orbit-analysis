from typing import Callable

import numpy as np

from .two_way_range import two_way_range


def average_two_way_range_rate(
    satellite_position: Callable[[float], np.ndarray],
    station_position: Callable[[float], np.ndarray],
    reception_time: float,
    count_interval: float,
    c_light: float,
    tolerance: float = 1.0e-9,
    max_iterations: int = 10,
) -> float:
    """
    Model an average two-way range rate measurement.

    A Doppler counter accumulates cycles over a count interval, so
    the measurement is the difference of two rigorous two-way range
    values divided by that interval. The time tag refers to the end
    of the count interval.

    Four light-time iterations are involved, two per range value.

    Parameters
    ----------
    satellite_position : Callable
        Function returning the inertial satellite position [km].

    station_position : Callable
        Function returning the inertial station position [km].

    reception_time : float
        Time tag of the measurement, at the end of the count
        interval [s].

    count_interval : float
        Doppler count interval [s].

    c_light : float
        Speed of light [km/s].

    tolerance : float
        Convergence threshold on the signal path length [km].

    max_iterations : int
        Maximum number of fixed-point iterations per leg.

    Returns
    -------
    float
        Average two-way range rate over the count interval [km/s].
    """

    range_end, _, _ = two_way_range(
        satellite_position,
        station_position,
        reception_time,
        c_light,
        tolerance,
        max_iterations,
    )

    range_start, _, _ = two_way_range(
        satellite_position,
        station_position,
        reception_time - count_interval,
        c_light,
        tolerance,
        max_iterations,
    )

    average_range_rate = (
        (range_end - range_start)
        / count_interval
    )

    return float(average_range_rate)