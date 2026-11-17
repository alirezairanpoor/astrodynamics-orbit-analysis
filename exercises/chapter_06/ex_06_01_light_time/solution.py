import numpy as np

from satellite_orbits.constants import (
    C_LIGHT,
    FLATTENING_EARTH,
    MU_EARTH,
    OMEGA_EARTH,
    R_EARTH,
    SECONDS_PER_DAY,
)
from satellite_orbits.astrodynamics.propagation.propagate_keplerian_state import (
    propagate_keplerian_state,
)
from satellite_orbits.coordinates import geodetic_to_ecef
from satellite_orbits.observations import (
    ground_station_state,
    light_time_downleg,
    light_time_upleg,
)
from satellite_orbits.time_systems import mjd_to_greenwich_hour_angle


def main() -> None:
    # Inputs: Exercise 6.1 scenario.
    # Source: Montenbruck and Gill, Exercise 6.1, page 229.
    #
    # Polar satellite in a circular orbit at 960 km altitude.
    # The ascending node passage is taken as the element epoch,
    # so the mean anomaly and the argument of perigee are zero.
    semi_major_axis = R_EARTH + 960.0
    eccentricity = 0.0
    inclination = np.radians(97.0)
    raan = np.radians(130.7)
    argument_of_perigee = 0.0
    mean_anomaly_epoch = 0.0

    # Ascending node passage: 1997/01/01 0h UTC.
    mjd_epoch_utc = 50449.0

    # Ground station in central Europe. No height is specified in
    # the exercise, so the station is placed on the ellipsoid.
    station_longitude = np.radians(11.0)
    station_latitude = np.radians(48.0)
    station_height = 0.0

    # Evaluation times past the ascending node passage.
    evaluation_minutes = (6.0, 15.0, 24.0)

    station_r_ecef = geodetic_to_ecef(
        station_longitude,
        station_latitude,
        station_height,
        R_EARTH,
        FLATTENING_EARTH,
    )

    def satellite_position(t):
        r_eci, _ = propagate_keplerian_state(
            semi_major_axis,
            eccentricity,
            inclination,
            raan,
            argument_of_perigee,
            mean_anomaly_epoch,
            t,
            MU_EARTH,
        )
        return r_eci

    def station_position(t):
        greenwich_hour_angle = mjd_to_greenwich_hour_angle(
            mjd_epoch_utc + t / SECONDS_PER_DAY
        )
        r_eci, _ = ground_station_state(
            station_r_ecef,
            greenwich_hour_angle,
            OMEGA_EARTH,
        )
        return r_eci

    # Reference outputs: rounded values printed in the textbook.
    book_distance_m = (3644878.6, 1167049.9, 4110978.7)
    book_down_first_m = (75.8, -6.4, -87.0)
    book_down_second_mm = (1.6, 0.0, 1.8)
    book_up_first_m = (1.8, 0.5, -0.7)
    book_range_m = (3644956.2, 1167044.0, 4110891.0)

    computed_distances_m = []
    computed_down_first_m = []
    computed_down_second_mm = []
    computed_up_first_m = []
    computed_ranges_m = []

    print("Exercise 6.1 - Light-Time Iteration")
    print("Satellite: circular polar orbit, 960 km altitude")
    print("Station: 11 deg East, 48 deg North, on the ellipsoid")
    print("Model: Keplerian motion, no precession or nutation")
    print()
    print(
        "UTC        Distance     Down I-1  Down I-2   Up I-1"
        "      Range"
    )
    print(
        "hh:mm:ss        [m]          [m]      [mm]      [m]"
        "        [m]"
    )
    print("-" * 70)

    for minutes in evaluation_minutes:

        reception_time = minutes * 60.0

        # Downleg: the signal arrives at the station at the
        # reception time, so the station position is fixed and the
        # satellite is traced back along its orbit.
        travel_time_down, down_distances = light_time_downleg(
            satellite_position,
            station_position(reception_time),
            reception_time,
            C_LIGHT,
        )

        transponding_time = reception_time - travel_time_down

        # Upleg: the satellite transponds at the epoch fixed by the
        # downleg, so the satellite position is fixed and the
        # station is traced back along the Earth rotation.
        travel_time_up, up_distances = light_time_upleg(
            station_position,
            satellite_position(transponding_time),
            transponding_time,
            C_LIGHT,
        )

        range_two_way = (
            0.5
            * C_LIGHT
            * (travel_time_down + travel_time_up)
        )

        geometric_distance_m = down_distances[0] * 1000.0
        down_first_m = (
            down_distances[1] - down_distances[0]
        ) * 1000.0
        down_second_mm = (
            down_distances[2] - down_distances[1]
        ) * 1.0e6
        up_first_m = (
            up_distances[1] - up_distances[0]
        ) * 1000.0
        range_m = range_two_way * 1000.0

        computed_distances_m.append(geometric_distance_m)
        computed_down_first_m.append(down_first_m)
        computed_down_second_mm.append(down_second_mm)
        computed_up_first_m.append(up_first_m)
        computed_ranges_m.append(range_m)

        print(
            f"00:{int(minutes):02d}:00.0 "
            f"{geometric_distance_m:12.1f} "
            f"{down_first_m:10.1f} "
            f"{down_second_mm:9.1f} "
            f"{up_first_m:8.1f} "
            f"{range_m:12.1f}"
        )

    # Verification: the light-time corrections are differences of
    # two evaluations a few milliseconds apart, so a constant
    # rotation offset cancels out of them. They are therefore
    # compared against the textbook at its printed resolution,
    # unlike the absolute distances below.
    np.testing.assert_allclose(
        computed_down_first_m,
        book_down_first_m,
        rtol=0.0,
        atol=0.05,
    )

    np.testing.assert_allclose(
        computed_down_second_mm,
        book_down_second_mm,
        rtol=0.0,
        atol=0.05,
    )

    np.testing.assert_allclose(
        computed_up_first_m,
        book_up_first_m,
        rtol=0.0,
        atol=0.05,
    )

    # The absolute distances and ranges carry a systematic offset
    # with respect to the textbook, discussed in the README. It
    # stays below 2 m and is dominated by the sidereal-time
    # convention, not by the light-time algorithm.
    distance_residuals_m = (
        np.array(computed_distances_m)
        - np.array(book_distance_m)
    )
    range_residuals_m = (
        np.array(computed_ranges_m)
        - np.array(book_range_m)
    )

    print()
    print("Residuals with respect to the textbook values [m]:")
    print(
        "Geometric distance: "
        + np.array2string(distance_residuals_m, precision=2)
    )
    print(
        "Two-way range:      "
        + np.array2string(range_residuals_m, precision=2)
    )

    if np.max(np.abs(distance_residuals_m)) > 2.5:
        raise AssertionError("Distance residual check failed.")

    if np.max(np.abs(range_residuals_m)) > 2.5:
        raise AssertionError("Range residual check failed.")

    print()
    print("Light-time correction comparison: passed")
    print("Absolute distance residual bound: passed")


if __name__ == "__main__":
    main()
