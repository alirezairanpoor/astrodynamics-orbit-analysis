import numpy as np


def third_body_acceleration(
    r_eci: np.ndarray,
    r_body_eci: np.ndarray,
    mu_body: float,
) -> np.ndarray:
    """
    Compute third-body acceleration relative to Earth.

    Inputs:
        r_eci:
            Earth-to-satellite position [km].
        r_body_eci:
            Earth-to-third-body position [km].
        mu_body:
            Third-body gravitational parameter [km^3/s^2].

    Both positions must refer to the same epoch and ECI frame.

    Output:
        Relative perturbing acceleration in ECI [km/s^2].
    """

    r_eci = np.asarray(r_eci, dtype=float)
    r_body_eci = np.asarray(r_body_eci, dtype=float)

    if r_eci.shape != (3,) or r_body_eci.shape != (3,):
        raise ValueError("Position vectors must have shape (3,).")

    # Satellite-to-third-body vector.
    satellite_to_body_eci = r_body_eci - r_eci

    distance_satellite_body = np.linalg.norm(
        satellite_to_body_eci
    )

    distance_earth_body = np.linalg.norm(r_body_eci)

    if (
        distance_satellite_body == 0.0
        or distance_earth_body == 0.0
    ):
        raise ValueError("Distances to the third body must be nonzero.")

    a_eci = mu_body * (
        satellite_to_body_eci / distance_satellite_body**3
        - r_body_eci / distance_earth_body**3
    )

    return a_eci
