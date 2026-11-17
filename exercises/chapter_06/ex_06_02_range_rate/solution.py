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
    average_two_way_range_rate,
    ground_station_state,
    instantaneous_range_rate,
)
from satellite_orbits.time_systems import mjd_to_greenwich_hour_angle


def main() -> None:
    # Inputs: same scenario as Exercise 6.1.
    # Source: Montenbruck and Gill, Exercise 6.2, page 229.
    semi_major_axis = R_EARTH + 960.0
    eccentricity = 0.0
    inclination = np.radians(97.0)
    raan = np.radians(130.7)
    argument_of_perigee = 0.0
    mean_anomaly_epoch = 0.0

    mjd_epoch_utc = 50449.0

    station_longitude = np.radians(11.0)
    station_latitude = np.radians(48.0)
    station_height = 0.0

    # Doppler count interval specified by the exercise.
    count_interval = 1.0

    evaluation_minutes = (6.0, 9.0, 12.0, 15.0, 18.0, 21.0, 24.0)

    station_r_ecef = geodetic_to_ecef(
        station_longitude,
        station_latitude,
        station_height,
        R_EARTH,
        FLATTENING_EARTH,
    )

    def satellite_state(t):
        return propagate_keplerian_state(
            semi_major_axis,
            eccentricity,
            inclination,
            raan,
            argument_of_perigee,
            mean_anomaly_epoch,
            t,
            MU_EARTH,
        )

    def station_state(t):
        greenwich_hour_angle = mjd_to_greenwich_hour_angle(
            mjd_epoch_utc + t / SECONDS_PER_DAY
        )
        return ground_station_state(
            station_r_ecef,
            greenwich_hour_angle,
            OMEGA_EARTH,
        )

    def satellite_position(t):
        return satellite_state(t)[0]

    def station_position(t):
        return station_state(t)[0]

    # Reference outputs: rounded values printed in the textbook.
    book_average_ms = (
        -6385.694, -6113.312, -4575.872, 1488.564,
        5500.474, 6286.690, 6392.679,
    )
    book_instantaneous_ms = (
        -6385.554, -6113.161, -4575.720, 1488.732,
        5500.635, 6286.839, 6392.813,
    )
    book_difference_ms = (
        -0.140, -0.151, -0.153, -0.168,
        -0.161, -0.149, -0.134,
    )

    computed_average_ms = []
    computed_instantaneous_ms = []
    computed_difference_ms = []

    print("Exercise 6.2 - Range Rate Modeling")
    print(f"Doppler count interval: {count_interval:.1f} s")
    print("Average: rigorous two-way range difference")
    print("Instantaneous: geometric, at the mid of the interval")
    print()
    print(
        "UTC        Average  Instantaneous      Delta"
    )
    print(
        "hh:mm:ss     [m/s]          [m/s]      [m/s]"
    )
    print("-" * 46)

    for minutes in evaluation_minutes:

        reception_time = minutes * 60.0

        # The measurement time tag refers to the end of the count
        # interval, and the full light-time model is used.
        average_ms = average_two_way_range_rate(
            satellite_position,
            station_position,
            reception_time,
            count_interval,
            C_LIGHT,
        ) * 1000.0

        # The instantaneous value is evaluated at the mid of the
        # count interval, where the Taylor expansion of the
        # average is centred.
        mid_time = reception_time - 0.5 * count_interval

        r_satellite, v_satellite = satellite_state(mid_time)
        r_station, v_station = station_state(mid_time)

        instantaneous_ms = instantaneous_range_rate(
            r_satellite,
            v_satellite,
            r_station,
            v_station,
        ) * 1000.0

        difference_ms = average_ms - instantaneous_ms

        computed_average_ms.append(average_ms)
        computed_instantaneous_ms.append(instantaneous_ms)
        computed_difference_ms.append(difference_ms)

        print(
            f"00:{int(minutes):02d}:00.0 "
            f"{average_ms:12.3f} "
            f"{instantaneous_ms:14.3f} "
            f"{difference_ms:10.3f}"
        )

    # Verification: the difference between the two models is the
    # quantity the exercise asks for. It is a difference of two
    # values sharing the same reference-frame bias, so it is
    # compared at the printed resolution of the textbook.
    np.testing.assert_allclose(
        computed_difference_ms,
        book_difference_ms,
        rtol=0.0,
        atol=0.002,
    )

    # The absolute range rates inherit the systematic offset
    # discussed in the Exercise 6.1 README. It stays far below the
    # effect under study here.
    average_residuals_ms = (
        np.array(computed_average_ms)
        - np.array(book_average_ms)
    )
    instantaneous_residuals_ms = (
        np.array(computed_instantaneous_ms)
        - np.array(book_instantaneous_ms)
    )

    print()
    print("Residuals with respect to the textbook values [m/s]:")
    print(
        "Average:       "
        + np.array2string(average_residuals_ms, precision=4)
    )
    print(
        "Instantaneous: "
        + np.array2string(instantaneous_residuals_ms, precision=4)
    )

    if np.max(np.abs(average_residuals_ms)) > 0.02:
        raise AssertionError("Average range rate check failed.")

    if np.max(np.abs(instantaneous_residuals_ms)) > 0.02:
        raise AssertionError("Instantaneous range rate check failed.")

    print()
    print("Model difference comparison: passed")
    print("Absolute range rate residual bound: passed")
    print()
    print(
        "The instantaneous model is in error by 0.13-0.17 m/s, "
        "which is two orders of magnitude above the 1 mm/s "
        "accuracy of real Doppler systems."
    )


if __name__ == "__main__":
    main()
