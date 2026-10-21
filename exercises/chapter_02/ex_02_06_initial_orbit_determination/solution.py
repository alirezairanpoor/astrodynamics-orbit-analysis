from math import radians

import numpy as np

from math import radians, degrees

from satellite_orbits.constants import (
    MU_EARTH,
    R_EARTH,
    FLATTENING_EARTH,
)
from satellite_orbits.coordinates import (
    az_el_range_to_enz,
    enz_to_ecef,
    ecef_to_geodetic,
    ecef_to_eci,
)
from satellite_orbits.time_systems import (
    calendar_to_mjd,
    mjd_to_greenwich_hour_angle,
)
from satellite_orbits.astrodynamics import (
    gauss_two_position_elements,
)

# ============================================================
# GROUND STATION
# ============================================================

r_station_ecef = np.array([
    1344.143,
    6068.601,
    1429.311,
])  # km


# Convert station ECEF coordinates
# to geodetic longitude and latitude
(
    station_longitude,
    station_latitude,
    station_height,
) = ecef_to_geodetic(
    r_station_ecef,
    R_EARTH,
    FLATTENING_EARTH,
)


# ============================================================
# OBSERVATION 1
# ============================================================

azimuth_1 = radians(132.67)
elevation_1 = radians(32.44)
range_1 = 16945.450  # km

mjd_1 = calendar_to_mjd(
    1999,
    4,
    2,
    0,
    30,
    0.0,
)

theta_1 = mjd_to_greenwich_hour_angle(
    mjd_1
)


# Local observation -> ENZ
s_enz_1 = az_el_range_to_enz(
    azimuth_1,
    elevation_1,
    range_1,
)

# ENZ -> ECEF relative vector
s_ecef_1 = enz_to_ecef(
    s_enz_1,
    station_longitude,
    station_latitude,
)

# Station + relative vector
# -> satellite ECEF position
r_ecef_1 = (
    r_station_ecef
    + s_ecef_1
)

# ECEF -> ECI
r_eci_1 = ecef_to_eci(
    r_ecef_1,
    theta_1,
)


# ============================================================
# OBSERVATION 2
# ============================================================

azimuth_2 = radians(123.08)
elevation_2 = radians(50.06)
range_2 = 37350.340  # km

mjd_2 = calendar_to_mjd(
    1999,
    4,
    2,
    3,
    0,
    0.0,
)

theta_2 = mjd_to_greenwich_hour_angle(
    mjd_2
)


# Local observation -> ENZ
s_enz_2 = az_el_range_to_enz(
    azimuth_2,
    elevation_2,
    range_2,
)

# ENZ -> ECEF relative vector
s_ecef_2 = enz_to_ecef(
    s_enz_2,
    station_longitude,
    station_latitude,
)

# Station + relative vector
# -> satellite ECEF position
r_ecef_2 = (
    r_station_ecef
    + s_ecef_2
)

# ECEF -> ECI
r_eci_2 = ecef_to_eci(
    r_ecef_2,
    theta_2,
)

# ============================================================
# GAUSS INITIAL ORBIT DETERMINATION
# ============================================================

delta_t = (
    2.5
    * 3600
)

(
    a,
    e,
    i,
    RAAN,
    AOP,
    M,
) = gauss_two_position_elements(
    r_eci_1,
    r_eci_2,
    delta_t,
    MU_EARTH,
)
# ============================================================
# OUTPUT
# ============================================================

print("Exercise 2.6 — Initial Orbit Determination")
print("--------------------------------------------")

print("\nECI position at observation 1:")
print(r_eci_1)

print("\nECI position at observation 2:")
print(r_eci_2)

print("\nOrbital elements:")
print(f"a    = {a:.3f} km")
print(f"e    = {e:.7f}")
print(f"i    = {degrees(i):.3f} deg")
print(f"RAAN = {degrees(RAAN):.3f} deg")
print(f"AOP  = {degrees(AOP):.3f} deg")
print(f"M    = {degrees(M):.3f} deg")
