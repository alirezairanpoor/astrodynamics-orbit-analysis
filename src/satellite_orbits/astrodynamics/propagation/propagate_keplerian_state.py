import numpy as np

from satellite_orbits.astrodynamics.keplerian.propagate_mean_anomaly import (
    propagate_mean_anomaly,
)
from satellite_orbits.astrodynamics.keplerian.elements_to_states import (
    elements_to_states,
)


def propagate_keplerian_state(
    a: float,
    e: float,
    i: float,
    RAAN: float,
    AOP: float,
    M0: float,
    elapsed_time: float,
    mu: float,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Propagate an unperturbed Keplerian orbit and return the
    inertial state at the requested time.

    This is the analytical counterpart of the numerical
    integrators in this package: both advance a state in time,
    one in closed form and one step by step.

    Parameters
    ----------
    a : float
        Semi-major axis [km].

    e : float
        Eccentricity [-].

    i : float
        Inclination [rad].

    RAAN : float
        Right ascension of the ascending node [rad].

    AOP : float
        Argument of perigee [rad].

    M0 : float
        Mean anomaly at the reference epoch [rad].

    elapsed_time : float
        Time elapsed since the reference epoch [s].

    mu : float
        Gravitational parameter [km^3/s^2].

    Returns
    -------
    tuple
        r_eci : Position in the inertial frame [km].
        v_eci : Velocity in the inertial frame [km/s].
    """

    M = propagate_mean_anomaly(
        a,
        M0,
        elapsed_time,
        mu,
    )

    r_eci, v_eci = elements_to_states(
        a,
        e,
        i,
        RAAN,
        AOP,
        M,
        mu,
    )

    return (
        r_eci,
        v_eci,
    )
