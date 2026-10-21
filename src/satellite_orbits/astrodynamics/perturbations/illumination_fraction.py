from math import asin, acos, sqrt, pi

import numpy as np


def illumination_fraction(
    r_eci: np.ndarray,
    r_sun_eci: np.ndarray,
    earth_radius: float,
    sun_radius: float = 696000.0,
) -> float:
    """
    Compute the visible fraction of the solar disk.

    Uses the apparent-disk overlap model of Section 3.4.2.

    Inputs:
        r_eci:
            Earth-to-satellite position [km].
        r_sun_eci:
            Earth-to-Sun position [km].
        earth_radius:
            Spherical Earth radius [km].
        sun_radius:
            Spherical Sun radius [km].

    Both vectors must use the same frame and epoch.

    Output:
        Illumination fraction between 0 and 1.

    Assumptions:
        Spherical Earth and Sun.
        Uniform solar disk brightness.
        Atmospheric refraction and attenuation neglected.
    """

    r_eci = np.asarray(r_eci, dtype=float)
    r_sun_eci = np.asarray(r_sun_eci, dtype=float)

    if r_eci.shape != (3,) or r_sun_eci.shape != (3,):
        raise ValueError("Position vectors must have shape (3,).")

    if earth_radius <= 0.0 or sun_radius <= 0.0:
        raise ValueError("Radii must be positive.")

    # Directions as seen from the satellite.
    satellite_to_earth = -r_eci
    satellite_to_sun = r_sun_eci - r_eci

    distance_earth = np.linalg.norm(satellite_to_earth)
    distance_sun = np.linalg.norm(satellite_to_sun)

    if distance_earth <= earth_radius:
        raise ValueError("Satellite must be outside the Earth.")

    if distance_sun <= sun_radius:
        raise ValueError("Satellite must be outside the Sun.")

    # --------------------------------------------------------
    # APPARENT ANGULAR RADII AND CENTER SEPARATION
    # --------------------------------------------------------

    sun_angle = asin(sun_radius / distance_sun)
    earth_angle = asin(earth_radius / distance_earth)

    cos_separation = np.dot(
        satellite_to_earth / distance_earth,
        satellite_to_sun / distance_sun,
    )

    separation = acos(
        float(np.clip(cos_separation, -1.0, 1.0))
    )

    # --------------------------------------------------------
    # NO OVERLAP: FULL SUNLIGHT
    # --------------------------------------------------------

    if separation >= sun_angle + earth_angle:
        return 1.0

    # --------------------------------------------------------
    # EARTH COMPLETELY COVERS THE SUN
    # --------------------------------------------------------

    if earth_angle >= separation + sun_angle:
        return 0.0

    # --------------------------------------------------------
    # EARTH DISK ENTIRELY INSIDE THE SUN DISK
    # --------------------------------------------------------

    if sun_angle >= separation + earth_angle:
        return float(
            np.clip(
                1.0 - (earth_angle / sun_angle) ** 2,
                0.0,
                1.0,
            )
        )

    # --------------------------------------------------------
    # PARTIAL OVERLAP
    # --------------------------------------------------------

    x = (
        separation**2
        + sun_angle**2
        - earth_angle**2
    ) / (2.0 * separation)

    y = sqrt(max(0.0, sun_angle**2 - x**2))

    sun_sector_angle = acos(
        float(np.clip(x / sun_angle, -1.0, 1.0))
    )

    earth_sector_angle = acos(
        float(
            np.clip(
                (separation - x) / earth_angle,
                -1.0,
                1.0,
            )
        )
    )

    overlap_area = (
        sun_angle**2 * sun_sector_angle
        + earth_angle**2 * earth_sector_angle
        - separation * y
    )

    visible_fraction = (
        1.0
        - overlap_area / (pi * sun_angle**2)
    )

    return float(np.clip(visible_fraction, 0.0, 1.0))
