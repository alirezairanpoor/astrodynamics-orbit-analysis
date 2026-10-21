from math import radians, degrees

from satellite_orbits.constants import MU_EARTH, R_EARTH
from satellite_orbits.astrodynamics import (
    propagate_mean_anomaly,
    elements_to_states,
)
from satellite_orbits.coordinates import (
    eci_to_ecef,
    ecef_to_enz,
    enz_to_az_el_range,
)
from satellite_orbits.time_systems import (
    mjd_to_greenwich_hour_angle,
)


# ============================================================
# USER INPUTS
# ============================================================

# Orbital elements
a = R_EARTH + 960.0  # km
e = 0.0
i = radians(97.0)
RAAN = radians(130.7)

# Circular orbit
AOP = 0.0
M0 = 0.0

# Ground station
station_longitude = radians(11.0)
station_latitude = radians(48.0)

# Reference epoch
reference_mjd = 50449.0

# Simulation
duration = 30 * 60  # s
time_step = 60      # s


# ============================================================
# SOLUTION
# ============================================================

results = []

for elapsed_time in range(
    0,
    duration + time_step,
    time_step,
):

    # --------------------------------------------------------
    # 1. Propagate mean anomaly
    # --------------------------------------------------------

    M = propagate_mean_anomaly(
        a,
        M0,
        elapsed_time,
        MU_EARTH,
    )

    # --------------------------------------------------------
    # 2. Orbital elements -> ECI state
    # --------------------------------------------------------

    r_eci, v_eci = elements_to_states(
        a,
        e,
        i,
        RAAN,
        AOP,
        M,
        MU_EARTH,
    )

    # --------------------------------------------------------
    # 3. Current time -> Greenwich Hour Angle
    # --------------------------------------------------------

    current_mjd = (
        reference_mjd
        + elapsed_time / 86400.0
    )

    greenwich_hour_angle = (
        mjd_to_greenwich_hour_angle(
            current_mjd
        )
    )

    # --------------------------------------------------------
    # 4. ECI -> ECEF
    # --------------------------------------------------------

    r_ecef = eci_to_ecef(
        r_eci,
        greenwich_hour_angle,
    )

    # --------------------------------------------------------
    # 5. ECEF -> Local ENZ
    # --------------------------------------------------------

    s_enz = ecef_to_enz(
        r_ecef,
        station_longitude,
        station_latitude,
        R_EARTH,
    )

    # --------------------------------------------------------
    # 6. ENZ -> Azimuth / Elevation / Range
    # --------------------------------------------------------

    azimuth, elevation, range_km = (
        enz_to_az_el_range(
            s_enz
        )
    )

    results.append(
        (
            elapsed_time,
            azimuth,
            elevation,
            range_km,
        )
    )


# ============================================================
# VISIBILITY
# ============================================================

rise_time = None
set_time = None

for index in range(len(results) - 1):

    time_1 = results[index][0]
    elevation_1 = results[index][2]

    time_2 = results[index + 1][0]
    elevation_2 = results[index + 1][2]

    # Satellite crosses upward through the horizon
    if elevation_1 <= 0 < elevation_2:

        rise_time = (
            time_1
            + (0 - elevation_1)
            / (elevation_2 - elevation_1)
            * (time_2 - time_1)
        )

    # Satellite crosses downward through the horizon
    if elevation_1 >= 0 > elevation_2:

        set_time = (
            time_1
            + (0 - elevation_1)
            / (elevation_2 - elevation_1)
            * (time_2 - time_1)
        )


# ============================================================
# OUTPUT
# ============================================================

print("Exercise 2.4 — Topocentric Satellite Motion")
print("------------------------------------------------")
print("Time [min] | Azimuth [deg] | Elevation [deg] | Range [km]")
print("------------------------------------------------------------")

for (
    elapsed_time,
    azimuth,
    elevation,
    range_km,
) in results:

    print(
        f"{elapsed_time / 60:10.1f} | "
        f"{degrees(azimuth):13.1f} | "
        f"{degrees(elevation):15.1f} | "
        f"{range_km:10.1f}"
    )


# Maximum elevation
max_elevation = max(
    elevation
    for _, _, elevation, _ in results
)

print(
    "\nMaximum elevation: "
    f"{degrees(max_elevation):.1f} deg"
)


# Visibility duration
if rise_time is not None and set_time is not None:

    visibility_duration = (
        set_time - rise_time
    )

    print(
        f"Rise time: "
        f"{rise_time / 60:.2f} min"
    )

    print(
        f"Set time: "
        f"{set_time / 60:.2f} min"
    )

    print(
        f"Visibility duration: "
        f"{visibility_duration / 60:.2f} min"
    )
