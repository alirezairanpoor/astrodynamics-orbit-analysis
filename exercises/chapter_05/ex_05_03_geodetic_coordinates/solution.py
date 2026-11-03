from math import cos, degrees, sin, sqrt

import numpy as np

from satellite_orbits.coordinates.ecef_to_geodetic import (
    ecef_to_geodetic,
)


def main() -> None:
    # Input: Diego Garcia receiver position in WGS84(G873) [m].
    # Coordinate epoch: 1997.0.
    # Source: Montenbruck and Gill, Exercise 5.3, page 192.
    r_ecef_m = np.array([
        1917032.190,
        6029782.349,
        -801376.113,
    ])

    # Inputs: WGS84 reference ellipsoid parameters.
    # Source: Table 5.4, page 189.
    equatorial_radius_km = 6378.137
    flattening = 1.0 / 298.257223563

    # Input conversion: the reusable function expects kilometres.
    r_ecef_km = r_ecef_m / 1000.0

    # Outputs: longitude [rad], geodetic latitude [rad],
    # and ellipsoidal height [km], in that order.
    longitude_rad, latitude_rad, height_km = ecef_to_geodetic(
        r_ecef=r_ecef_km,
        earth_radius=equatorial_radius_km,
        flattening=flattening,
    )

    # Output conversion: degrees and metres for the book comparison.
    longitude_deg = degrees(longitude_rad)
    latitude_deg = degrees(latitude_rad)
    height_m = height_km * 1000.0

    # Reference outputs: rounded values printed in Exercise 5.3.
    longitude_book_deg = 72.36312094
    latitude_book_deg = -7.26654999
    height_book_m = -63.667

    np.testing.assert_allclose(
        longitude_deg,
        longitude_book_deg,
        rtol=0.0,
        atol=1.0e-8,
    )

    np.testing.assert_allclose(
        latitude_deg,
        latitude_book_deg,
        rtol=0.0,
        atol=1.0e-8,
    )

    np.testing.assert_allclose(
        height_m,
        height_book_m,
        rtol=0.0,
        atol=0.001,
    )

    # Verification: reconstruct the Cartesian coordinates using
    # equations 5.83-5.84, independently of the inverse iteration.
    eccentricity_squared = flattening * (2.0 - flattening)

    prime_vertical_radius_km = (
        equatorial_radius_km
        / sqrt(
            1.0
            - eccentricity_squared * sin(latitude_rad)**2
        )
    )

    r_reconstructed_km = np.array([
        (
            (prime_vertical_radius_km + height_km)
            * cos(latitude_rad)
            * cos(longitude_rad)
        ),
        (
            (prime_vertical_radius_km + height_km)
            * cos(latitude_rad)
            * sin(longitude_rad)
        ),
        (
            (
                (1.0 - eccentricity_squared)
                * prime_vertical_radius_km
                + height_km
            )
            * sin(latitude_rad)
        ),
    ])

    reconstruction_error_m = (
        np.linalg.norm(r_reconstructed_km - r_ecef_km)
        * 1000.0
    )

    if reconstruction_error_m > 1.0e-6:
        raise AssertionError("Cartesian reconstruction check failed.")

    print("Exercise 5.3 — Geodetic Coordinates")
    print("Station: Diego Garcia")
    print("Coordinate frame: WGS84(G873)")
    print("Coordinate epoch: 1997.0")
    print("Reference ellipsoid: WGS84")

    print("\nInput Cartesian position [m]:")
    print(np.array2string(r_ecef_m, precision=3))

    print("\nGeodetic coordinates:")
    print(f"East longitude [deg]: {longitude_deg:+.8f}")
    print(f"Latitude       [deg]: {latitude_deg:+.8f}")
    print(f"Height           [m]: {height_m:+.3f}")

    print("\nDifferences from the rounded textbook values:")
    print(
        "Longitude [deg]: "
        f"{longitude_deg - longitude_book_deg:+.3e}"
    )
    print(
        "Latitude  [deg]: "
        f"{latitude_deg - latitude_book_deg:+.3e}"
    )
    print(
        "Height      [m]: "
        f"{height_m - height_book_m:+.3e}"
    )

    print(
        "\nCartesian reconstruction error [m]: "
        f"{reconstruction_error_m:.3e}"
    )

    print("\nTextbook coordinate comparison: passed")
    print("Cartesian reconstruction check: passed")


if __name__ == "__main__":
    main()
