import numpy as np


def solar_radiation_pressure_acceleration(
    r_eci: np.ndarray,
    r_sun_eci: np.ndarray,
    area: float,
    mass: float,
    C_R: float,
    illumination: float,
    solar_pressure: float = 4.560e-6,
    astronomical_unit: float = 149597870.691,
) -> np.ndarray:
    """
    Compute direct solar radiation pressure acceleration.

    Uses a constant effective area and radiation coefficient.

    Inputs:
        r_eci:
            Earth-to-satellite position [km].
        r_sun_eci:
            Earth-to-Sun position [km].
        area:
            Effective illuminated cross-sectional area [m^2].
        mass:
            Satellite mass [kg].
        C_R:
            Radiation pressure coefficient.
        illumination:
            Visible fraction of the solar disk, from 0 to 1.
        solar_pressure:
            Radiation pressure at 1 AU [N/m^2].
        astronomical_unit:
            Reference distance [km].

    Both position vectors must use the same ECI frame and epoch.

    Output:
        SRP acceleration in ECI [km/s^2].
    """

    r_eci = np.asarray(r_eci, dtype=float)
    r_sun_eci = np.asarray(r_sun_eci, dtype=float)

    if r_eci.shape != (3,) or r_sun_eci.shape != (3,):
        raise ValueError("Position vectors must have shape (3,).")

    if mass <= 0.0:
        raise ValueError("Mass must be positive.")

    if area < 0.0 or C_R < 0.0:
        raise ValueError("Area and C_R must be nonnegative.")

    if not 0.0 <= illumination <= 1.0:
        raise ValueError("Illumination must be between 0 and 1.")

    # Direction away from the Sun.
    sun_to_satellite_eci = r_eci - r_sun_eci

    distance = np.linalg.norm(sun_to_satellite_eci)

    if distance == 0.0:
        raise ValueError("Sun-to-satellite distance must be nonzero.")

    direction_eci = sun_to_satellite_eci / distance

    # Pressure scaled from 1 AU to the actual distance.
    pressure = solar_pressure * (
        astronomical_unit / distance
    ) ** 2

    # N/m^2 * m^2/kg = m/s^2
    acceleration_m_s2 = (
        illumination
        * C_R
        * pressure
        * area
        / mass
    )

    # Convert m/s^2 to km/s^2.
    a_eci = (
        acceleration_m_s2
        / 1000.0
        * direction_eci
    )

    return a_eci
