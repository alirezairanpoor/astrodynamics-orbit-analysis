import numpy as np

from satellite_orbits.constants import (
    FLATTENING_EARTH,
    MU_EARTH,
    R_EARTH,
    SECONDS_PER_DAY,
)
from satellite_orbits.astrodynamics.propagation.propagate_keplerian_state import (
    propagate_keplerian_state,
)
from satellite_orbits.coordinates import (
    ecef_to_enz_rotation,
    eci_to_ecef,
    enz_to_az_el_range,
    geodetic_to_ecef,
)
from satellite_orbits.observations import (
    elevation_refraction_correction,
    tropospheric_refractivity,
)
from satellite_orbits.time_systems import mjd_to_greenwich_hour_angle


def main() -> None:
    # Inputs: near-geostationary satellite.
    # Source: Montenbruck and Gill, Exercise 6.4, page 232.
    #
    # The argument of perigee and the mean anomaly are not given by
    # the exercise and are taken as zero at the epoch.
    semi_major_axis = 42164.0
    eccentricity = 0.00296
    inclination = np.radians(0.05)
    raan = np.radians(150.7)
    argument_of_perigee = 0.0
    mean_anomaly_epoch = 0.0

    # Epoch: 1997/01/01 0h UTC.
    mjd_epoch_utc = 50449.0

    # Ground station of Exercise 6.1.
    station_longitude = np.radians(11.0)
    station_latitude = np.radians(48.0)
    station_height = 0.0

    # Meteorological conditions specified by the exercise.
    air_pressure = 1024.0
    relative_humidity = 0.7
    temperature_1 = 283.0
    temperature_2 = 303.0

    evaluation_hours = np.arange(0.0, 24.1, 3.0)

    # Reference outputs printed in the textbook.
    #
    # The textbook reports E0 = 23.363 deg and a peak-to-peak
    # variation of about 0.058 deg. That variation is reproduced
    # here only when the eccentricity is set to zero, which
    # indicates that its daily motion was modelled without the
    # eccentric libration. See the README.
    book_elevation_epoch_deg = 23.363
    book_circular_spread_deg = 0.058

    station_r_ecef = geodetic_to_ecef(
        station_longitude,
        station_latitude,
        station_height,
        R_EARTH,
        FLATTENING_EARTH,
    )

    def geometric_elevation(elapsed_seconds, orbit_eccentricity):
        r_eci, _ = propagate_keplerian_state(
            semi_major_axis,
            orbit_eccentricity,
            inclination,
            raan,
            argument_of_perigee,
            mean_anomaly_epoch,
            elapsed_seconds,
            MU_EARTH,
        )

        greenwich_hour_angle = mjd_to_greenwich_hour_angle(
            mjd_epoch_utc
            + elapsed_seconds / SECONDS_PER_DAY
        )

        r_ecef = eci_to_ecef(
            r_eci,
            greenwich_hour_angle,
        )

        s_enz = ecef_to_enz_rotation(
            r_ecef - station_r_ecef,
            station_longitude,
            station_latitude,
        )

        _, elevation, _ = enz_to_az_el_range(s_enz)

        return elevation

    elevation_epoch = geometric_elevation(0.0, eccentricity)
    elevation_epoch_deg = np.degrees(elevation_epoch)

    print("Exercise 6.4 - Tropospheric Refraction")
    print("Station: 11 deg East, 48 deg North, on the ellipsoid")
    print(
        f"Conditions: p = {air_pressure:.0f} hPa, "
        f"fh = {relative_humidity:.1f}"
    )
    print(f"Elevation at the epoch: {elevation_epoch_deg:.3f} deg")
    print()

    for temperature in (temperature_1, temperature_2):
        dry, wet = tropospheric_refractivity(
            air_pressure,
            temperature,
            relative_humidity,
        )
        print(
            f"T = {temperature:5.1f} K  ->  "
            f"N1 = {dry:7.3f}, N2 = {wet:7.3f}, "
            f"N = {dry + wet:7.3f}"
        )

    print()
    print("UTC        E - E0     dE(T1)     dE(T2)")
    print("hh:mm:ss    [deg]      [deg]      [deg]")
    print("-" * 40)

    variation_deg = []
    correction_1_deg = []
    correction_2_deg = []

    for hours in evaluation_hours:

        elapsed_seconds = hours * 3600.0

        elevation = geometric_elevation(
            elapsed_seconds,
            eccentricity,
        )

        variation = (
            np.degrees(elevation)
            - elevation_epoch_deg
        )

        correction_1 = np.degrees(
            elevation_refraction_correction(
                elevation,
                air_pressure,
                temperature_1,
                relative_humidity,
            )
        )

        correction_2 = np.degrees(
            elevation_refraction_correction(
                elevation,
                air_pressure,
                temperature_2,
                relative_humidity,
            )
        )

        variation_deg.append(variation)
        correction_1_deg.append(correction_1)
        correction_2_deg.append(correction_2)

        print(
            f"{int(hours):02d}:00:00.0 "
            f"{variation:9.3f} "
            f"{correction_1:10.3f} "
            f"{correction_2:10.3f}"
        )

    variation_deg = np.array(variation_deg)
    correction_1_deg = np.array(correction_1_deg)
    correction_2_deg = np.array(correction_2_deg)

    maximum_index = int(np.argmax(variation_deg))
    minimum_index = int(np.argmin(variation_deg))

    print()
    print(
        "Peak-to-peak elevation variation [deg]: "
        f"{np.ptp(variation_deg):.3f}"
    )
    print(
        "Maximum elevation at "
        f"{int(evaluation_hours[maximum_index]):02d} h, "
        "minimum at "
        f"{int(evaluation_hours[minimum_index]):02d} h"
    )
    print(
        "Refraction variation over the day [deg]: "
        f"{np.ptp(correction_1_deg):.4f}"
    )
    print(
        "Effect of the 20 K temperature difference [deg]: "
        f"{abs(correction_2_deg[0] - correction_1_deg[0]):.4f}"
    )

    # The eccentricity makes the satellite librate in longitude
    # with an amplitude of 2e radians, which dominates the daily
    # elevation variation.
    libration_amplitude_deg = np.degrees(
        2.0 * eccentricity
    )

    print(
        "Longitude libration amplitude from e [deg]: "
        f"{libration_amplitude_deg:.3f}"
    )

    # Control case: a circular orbit isolates the inclination
    # contribution and recovers the textbook elevation.
    circular_elevations_deg = np.array([
        np.degrees(geometric_elevation(hours * 3600.0, 0.0))
        for hours in evaluation_hours
    ])

    print()
    print("Control case with e = 0:")
    print(
        "  Elevation at the epoch [deg]: "
        f"{circular_elevations_deg[0]:.3f}  "
        f"(textbook {book_elevation_epoch_deg:.3f})"
    )
    print(
        "  Peak-to-peak [deg]: "
        f"{np.ptp(circular_elevations_deg):.4f}  "
        f"(textbook {book_circular_spread_deg:.3f})"
    )

    # Verification of the orbital geometry. The epoch elevation is
    # checked against the textbook; the daily variation is checked
    # against the libration the stated eccentricity must produce,
    # not against the textbook table.
    if abs(elevation_epoch_deg - book_elevation_epoch_deg) > 0.05:
        raise AssertionError("Epoch elevation check failed.")

    if not 0.20 < np.ptp(variation_deg) < 0.35:
        raise AssertionError(
            "Daily elevation variation outside the expected range."
        )

    if abs(
        circular_elevations_deg[0] - book_elevation_epoch_deg
    ) > 0.01:
        raise AssertionError("Circular control case failed.")

    # Verification of the refraction model. The textbook prints
    # 0.051 and 0.042 deg, which this implementation reproduces,
    # but with the two temperatures interchanged relative to the
    # printed column headings. See the README.
    np.testing.assert_allclose(
        correction_1_deg,
        0.0425,
        rtol=0.0,
        atol=0.001,
    )

    np.testing.assert_allclose(
        correction_2_deg,
        0.0505,
        rtol=0.0,
        atol=0.001,
    )

    print()
    print("Epoch elevation comparison: passed")
    print("Daily variation consistent with the eccentricity: passed")
    print("Circular control case comparison: passed")
    print("Refraction correction magnitudes: passed")
    print()
    print(
        "The refraction correction is comparable in size to the "
        "daily elevation variation but stays essentially constant "
        "for fixed weather, varying by only 0.0005 deg over a day."
    )


if __name__ == "__main__":
    main()
