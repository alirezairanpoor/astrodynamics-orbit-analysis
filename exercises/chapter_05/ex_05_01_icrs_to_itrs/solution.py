from datetime import datetime

import numpy as np

from satellite_orbits.coordinates.icrs_to_itrs import (
    earth_rotation_matrix,
    icrs_to_itrs_matrix,
    nutation_matrix_iau1980,
    polar_motion_matrix,
    precession_matrix_iau1976,
)


def main() -> None:
    # Inputs: Exercise 5.1 epoch and historical Earth orientation data.
    # Source: Montenbruck and Gill, Exercise 5.1, page 190.
    epoch_utc = datetime(1999, 3, 4)
    mjd_origin = datetime(1858, 11, 17)

    mjd_utc = (
        epoch_utc - mjd_origin
    ).total_seconds() / 86400.0

    tt_minus_utc = 64.184
    ut1_minus_utc = 0.649232
    xp_arcsec = 0.06740
    yp_arcsec = 0.24173

    mjd_tt = mjd_utc + tt_minus_utc / 86400.0
    mjd_ut1 = mjd_utc + ut1_minus_utc / 86400.0

    P = precession_matrix_iau1976(mjd_tt)

    N, delta_psi, _, epsilon = nutation_matrix_iau1980(mjd_tt)

    Theta = earth_rotation_matrix(
        mjd_ut1,
        delta_psi,
        epsilon,
    )

    Pi = polar_motion_matrix(xp_arcsec, yp_arcsec)

    # Output: Complete coordinate transformation at the exercise epoch.
    U = icrs_to_itrs_matrix(
        mjd_utc=mjd_utc,
        tt_minus_utc=tt_minus_utc,
        ut1_minus_utc=ut1_minus_utc,
        xp_arcsec=xp_arcsec,
        yp_arcsec=yp_arcsec,
    )

    # Reference output: Rounded matrix printed in the textbook.
    U_book = np.array([
        [-0.94737803, 0.32011696, -0.00008431],
        [-0.32011696, -0.94737803, -0.00006363],
        [-0.00010024, -0.00003330, 0.99999999],
    ])

    print("Exercise 5.1 — ICRS to ITRS Transformation")
    print("Epoch: 1999-03-04 00:00:00 UTC")
    print(f"MJD UTC: {mjd_utc:.12f}")
    print(f"MJD TT:  {mjd_tt:.12f}")
    print(f"MJD UT1: {mjd_ut1:.12f}")

    for name, matrix in (
        ("P", P),
        ("N", N),
        ("Theta", Theta),
        ("Pi", Pi),
        ("U = Pi @ Theta @ N @ P", U),
    ):
        print(f"\n{name}:")
        print(np.array2string(
            matrix,
            precision=8,
            suppress_small=True,
        ))

    max_difference = np.max(np.abs(U - U_book))

    print(
        "\nMaximum absolute difference from the rounded book matrix:"
        f" {max_difference:.3e}"
    )

    if not np.allclose(U, U_book, rtol=0.0, atol=1.0e-8):
        raise AssertionError("Textbook matrix comparison failed.")

    print("Textbook matrix comparison: passed")

    # Input: Example J2000 position vector [km].
    r_icrs = np.array([7000.0, 0.0, 0.0])

    # Output: The same position expressed in ITRS [km].
    r_itrs = U @ r_icrs

    print("\nExample position in J2000 [km]:")
    print(r_icrs)

    print("\nExample position in ITRS [km]:")
    print(np.array2string(r_itrs, precision=6))


if __name__ == "__main__":
    main()
