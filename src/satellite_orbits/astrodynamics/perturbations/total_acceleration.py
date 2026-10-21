import numpy as np

from satellite_orbits.constants import (
    R_EARTH,
    FLATTENING_EARTH,
)
from satellite_orbits.coordinates.eci_to_ecef import eci_to_ecef
from satellite_orbits.coordinates.ecef_to_eci import ecef_to_eci
from satellite_orbits.coordinates.ecef_to_geodetic import (
    ecef_to_geodetic,
)

from .gravity_acceleration_cunningham import (
    gravity_acceleration_cunningham,
)
from .third_body_acceleration import third_body_acceleration
from .solar_radiation_pressure_acceleration import (
    solar_radiation_pressure_acceleration,
)
from .illumination_fraction import illumination_fraction
from .atmospheric_drag_acceleration import (
    atmospheric_drag_acceleration,
)
from .harris_priester_density import harris_priester_density


def total_acceleration(
    r_eci,
    v_eci,
    *,
    greenwich_hour_angle,
    r_sun_eci,
    r_moon_eci,
    C,
    S,
    n_max,
    m_max,
    mu_earth,
    gravity_reference_radius,
    mu_sun,
    mu_moon,
    area,
    mass,
    C_R,
    C_D,
    hp_exponent,
    include_sun=True,
    include_moon=True,
    include_srp=True,
    include_drag=True,
):
    """Return total acceleration in ECI [km/s^2].

    Positions: km.
    Velocity: km/s.
    Gravitational parameters: km^3/s^2.
    Greenwich hour angle: rad.
    Spacecraft area: m^2.
    Spacecraft mass: kg.

    C and S must be unnormalized gravity coefficients.
    C[0, 0] must equal 1 to include central gravity.

    Earth orientation uses the existing rotation-only
    ECI/ECEF transformation.
    """
    r_eci = np.asarray(r_eci, dtype=float)
    v_eci = np.asarray(v_eci, dtype=float)

    if r_eci.shape != (3,) or v_eci.shape != (3,):
        raise ValueError("State vectors must have shape (3,).")

    if not (
        np.all(np.isfinite(r_eci))
        and np.all(np.isfinite(v_eci))
    ):
        raise ValueError("State vectors must be finite.")

    if mass <= 0 or area < 0:
        raise ValueError("Mass must be positive and area nonnegative.")

    # 1. Earth gravity: evaluate in ECEF, then rotate to ECI.
    r_ecef = eci_to_ecef(r_eci, greenwich_hour_angle)

    a_earth_ecef = gravity_acceleration_cunningham(
        r_ecef,
        C,
        S,
        n_max,
        m_max,
        mu_earth,
        gravity_reference_radius,
    )

    a_total_eci = np.array(
        ecef_to_eci(a_earth_ecef, greenwich_hour_angle),
        dtype=float,
        copy=True,
    )

    # 2. Differential gravity due to the Sun.
    if include_sun:
        a_total_eci += third_body_acceleration(
            r_eci,
            r_sun_eci,
            mu_sun,
        )

    # 3. Differential gravity due to the Moon.
    if include_moon:
        a_total_eci += third_body_acceleration(
            r_eci,
            r_moon_eci,
            mu_moon,
        )

    # 4. Solar radiation pressure, including Earth's shadow.
    if include_srp:
        illumination = illumination_fraction(
            r_eci,
            r_sun_eci,
            R_EARTH,
        )

        a_total_eci += solar_radiation_pressure_acceleration(
            r_eci,
            r_sun_eci,
            area,
            mass,
            C_R,
            illumination,
        )

    # 5. Atmospheric drag using Harris-Priester density.
    if include_drag:
        _, _, h = ecef_to_geodetic(
            r_ecef,
            R_EARTH,
            FLATTENING_EARTH,
        )

        density = harris_priester_density(
            h,
            r_eci,
            r_sun_eci,
            exponent=hp_exponent,
        )

        a_total_eci += atmospheric_drag_acceleration(
            r_eci,
            v_eci,
            density,
            area,
            mass,
            C_D,
        )

    return a_total_eci
