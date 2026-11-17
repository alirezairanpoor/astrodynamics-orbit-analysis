import numpy as np


def instantaneous_range_rate(
    r_satellite_eci: np.ndarray,
    v_satellite_eci: np.ndarray,
    r_station_eci: np.ndarray,
    v_station_eci: np.ndarray,
) -> float:
    """
    Compute the instantaneous geometric range rate.

    The geometric distance is d0 = |r - R|, so its time derivative
    is the projection of the relative velocity onto the line of
    sight. Light-time effects are deliberately excluded, since this
    function represents the simplified Doppler model.

    A positive value means the satellite is receding.

    Parameters
    ----------
    r_satellite_eci : np.ndarray
        Satellite position in the inertial frame [km].

    v_satellite_eci : np.ndarray
        Satellite velocity in the inertial frame [km/s].

    r_station_eci : np.ndarray
        Station position in the inertial frame [km].

    v_station_eci : np.ndarray
        Station velocity in the inertial frame [km/s].

    Returns
    -------
    float
        Instantaneous range rate [km/s].
    """

    relative_position = (
        r_satellite_eci
        - r_station_eci
    )

    relative_velocity = (
        v_satellite_eci
        - v_station_eci
    )

    distance = np.linalg.norm(
        relative_position
    )

    range_rate = (
        np.dot(
            relative_position,
            relative_velocity,
        )
        / distance
    )

    return float(range_rate)
