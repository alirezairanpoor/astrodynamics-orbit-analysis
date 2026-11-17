import numpy as np

from satellite_orbits.coordinates.ecef_to_eci import ecef_to_eci


def ground_station_state(
    r_ecef: np.ndarray,
    greenwich_hour_angle: float,
    omega_earth: float,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute the inertial position and velocity of a ground station.

    The station is fixed in the Earth-fixed frame, so its inertial
    velocity is produced entirely by the Earth rotation.

    Parameters
    ----------
    r_ecef : np.ndarray
        Station position in ECEF [km].

    greenwich_hour_angle : float
        Greenwich Hour Angle at the epoch of interest [rad].

    omega_earth : float
        Earth rotation rate [rad/s].

    Returns
    -------
    tuple
        r_eci : Station position in ECI [km].
        v_eci : Station velocity in ECI [km/s].
    """

    r_eci = ecef_to_eci(
        r_ecef,
        greenwich_hour_angle,
    )

    omega_vector = np.array([
        0.0,
        0.0,
        omega_earth,
    ])

    v_eci = np.cross(
        omega_vector,
        r_eci,
    )

    return (
        r_eci,
        v_eci,
    )
