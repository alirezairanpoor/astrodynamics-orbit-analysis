from datetime import datetime, timedelta
from math import pi, sqrt

import numpy as np

from satellite_orbits.constants import MU_EARTH
from satellite_orbits.coordinates.icrs_itrs_state import (
    icrs_to_itrs_matrix_and_rate,
    icrs_to_itrs_state,
    itrs_to_icrs_state,
)


def main() -> None:
    # Inputs: GPS epoch and historical time offsets [s].
    # Source: Montenbruck and Gill, Exercise 5.2, pages 191-192.
    epoch_gps = datetime(1999, 3, 4, 0, 0, 0)

    tai_minus_utc = 32.0
    tai_minus_gps = 19.0

    gps_minus_utc = tai_minus_utc - tai_minus_gps
    tt_minus_utc = tai_minus_utc + 32.184
    ut1_minus_utc = 0.649232

    # Inputs: Pole coordinates supplied for this exercise [arcsec].
    xp_arcsec = 0.06740
    yp_arcsec = 0.24173

    # Input: PRN 15 position in WGS84(G873) [km].
    # WGS84(G873) is treated as ITRS for this exercise.
    r_itrs = np.array([
        19440.953805,
        16881.609273,
        -6777.115092,
    ])

    # Input: Time derivative of the Earth-fixed coordinates [km/s].
    v_itrs = np.array([
        -0.8111827456,
        -0.2573799137,
        -3.0689508125,
    ])

    epoch_utc = epoch_gps - timedelta(seconds=gps_minus_utc)
    mjd_origin = datetime(1858, 11, 17)

    mjd_utc = (
        epoch_utc - mjd_origin
    ).total_seconds() / 86400.0

    # Outputs: Coordinate transformation matrix and rate [1, 1/s].
    U, U_dot = icrs_to_itrs_matrix_and_rate(
        mjd_utc=mjd_utc,
        tt_minus_utc=tt_minus_utc,
        ut1_minus_utc=ut1_minus_utc,
        xp_arcsec=xp_arcsec,
        yp_arcsec=yp_arcsec,
    )

    # Outputs: J2000 position [km] and velocity [km/s].
    r_icrs, v_icrs = itrs_to_icrs_state(
        r_itrs,
        v_itrs,
        U,
        U_dot,
    )

    # Reference outputs: Rounded textbook position and velocity.
    r_book = np.array([
        -23830.593,
        -9747.074,
        -6779.829,
    ])

    v_book = np.array([
        1.561964,
        -1.754346,
        -3.068851,
    ])

    np.testing.assert_allclose(
        r_icrs,
        r_book,
        rtol=0.0,
        atol=0.001,
    )

    np.testing.assert_allclose(
        v_icrs,
        v_book,
        rtol=0.0,
        atol=1.0e-6,
    )

    # Outputs: Osculating-orbit diagnostics from the inertial state.
    # These calculations check the orbit characteristics requested
    # by the exercise; they are not a full element-conversion API.
    radius = np.linalg.norm(r_icrs)
    speed_squared = np.dot(v_icrs, v_icrs)

    h_vector = np.cross(r_icrs, v_icrs)
    h_norm = np.linalg.norm(h_vector)

    eccentricity_vector = (
        np.cross(v_icrs, h_vector) / MU_EARTH
        - r_icrs / radius
    )

    eccentricity = np.linalg.norm(eccentricity_vector)

    semimajor_axis = 1.0 / (
        2.0 / radius - speed_squared / MU_EARTH
    )

    inclination_deg = np.degrees(
        np.arccos(np.clip(h_vector[2] / h_norm, -1.0, 1.0))
    )

    period_hours = (
        2.0 * pi * sqrt(semimajor_axis**3 / MU_EARTH)
        / 3600.0
    )

    np.testing.assert_allclose(
        semimajor_axis, 26561.013, rtol=0.0, atol=0.01
    )

    np.testing.assert_allclose(
        eccentricity, 0.0070606, rtol=0.0, atol=1.0e-6
    )

    np.testing.assert_allclose(
        inclination_deg, 56.338, rtol=0.0, atol=0.001
    )

    if abs(period_hours - 12.0) > 0.1:
        raise AssertionError("The orbital period is not near 12 hours.")

    # Outputs: Recovered Earth-fixed state for a round-trip check.
    r_recovered, v_recovered = icrs_to_itrs_state(
        r_icrs,
        v_icrs,
        U,
        U_dot,
    )

    np.testing.assert_allclose(
        r_recovered, r_itrs, rtol=0.0, atol=1.0e-8
    )

    np.testing.assert_allclose(
        v_recovered, v_itrs, rtol=0.0, atol=1.0e-11
    )

    print("Exercise 5.2 — Velocity in the Earth-fixed Frame")
    print(f"Epoch GPS: {epoch_gps:%Y-%m-%d %H:%M:%S}")
    print(f"Epoch UTC: {epoch_utc:%Y-%m-%d %H:%M:%S}")
    print(f"MJD UTC: {mjd_utc:.12f}")

    print("\nU:")
    print(np.array2string(U, precision=8, suppress_small=True))

    print("\nU_dot [1/s]:")
    print(np.array2string(
        U_dot,
        precision=10,
        suppress_small=False,
    ))

    print("\nPosition J2000 [km]:")
    print(np.array2string(r_icrs, precision=6))

    print("\nVelocity J2000 [km/s]:")
    print(np.array2string(v_icrs, precision=9))

    print("\nOsculating-orbit diagnostics:")
    print(f"Semimajor axis [km]: {semimajor_axis:.3f}")
    print(f"Eccentricity:       {eccentricity:.7f}")
    print(f"Inclination [deg]:  {inclination_deg:.3f}")
    print(f"Keplerian period [h]: {period_hours:.6f}")

    print("\nTextbook state comparison: passed")
    print("Orbital-characteristic checks: passed")
    print("State transformation round-trip: passed")


if __name__ == "__main__":
    main()
