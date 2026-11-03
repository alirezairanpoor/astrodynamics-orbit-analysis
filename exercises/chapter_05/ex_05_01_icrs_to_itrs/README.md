# Exercise 5.1 — ICRS to ITRS Transformation

## Objective

Compute the celestial-to-terrestrial coordinate transformation for
1999-03-04 at 00:00:00 UTC, following Exercise 5.1 of Montenbruck
and Gill, Satellite Orbits, page 190.

The exercise uses the mean equator and equinox of J2000 as an
approximation to ICRS.

## Inputs

The historical Earth orientation data supplied by the textbook are:

- TT − UTC: 64.184 s
- UT1 − UTC: 0.649232 s
- xp: 0.06740 arcsec
- yp: 0.24173 arcsec

These values apply to the exercise epoch. Other epochs require
their corresponding time offsets and Earth orientation parameters.

## Calculation

The complete transformation is:

    U = Pi @ Theta @ N @ P

The matrices are applied to column vectors from right to left:

1. P: J2000 coordinates to mean-of-date coordinates.
2. N: Mean-of-date coordinates to true-of-date coordinates.
3. Theta: Earth rotation about the Celestial Ephemeris Pole.
4. Pi: Polar motion, completing the transformation to ITRS.

TT is used for precession and nutation.
UT1 is used for Earth rotation.

## Models

- IAU 1976 precession.
- Complete IAU 1980 nutation series, evaluated through PyERFA.
- GMST82, evaluated through PyERFA.
- First-order equation of the equinoxes: delta_psi * cos(epsilon).
- First-order polar-motion matrix.

The implementation follows the textbook model rather than a
modern IAU 2006/2000A transformation.

## Files

- `solution.py`: Exercise inputs, intermediate matrices, textbook
  comparison, and an example position transformation.
- `src/satellite_orbits/coordinates/icrs_to_itrs.py`:
  Reusable coordinate-transformation functions.

Dependencies: NumPy and PyERFA.
The PyERFA package is imported as `erfa`.

## Reusable Interface

    U = icrs_to_itrs_matrix(
        mjd_utc,
        tt_minus_utc,
        ut1_minus_utc,
        xp_arcsec,
        yp_arcsec,
    )

Inputs:

- mjd_utc: Modified Julian Date in UTC, in days.
- tt_minus_utc: TT − UTC, in seconds.
- ut1_minus_utc: UT1 − UTC, in seconds.
- xp_arcsec and yp_arcsec: IERS pole coordinates, in arcseconds.

Output:

- U: A dimensionless 3 × 3 transformation matrix.

To transform a position:

    r_itrs = U @ r_icrs

Both vectors describe the same physical position at the same epoch,
expressed in different reference frames. Their length units are
unchanged.

## Textbook Comparison

The rounded reference matrix is:

    [-0.94737803,  0.32011696, -0.00008431]
    [-0.32011696, -0.94737803, -0.00006363]
    [-0.00010024, -0.00003330,  0.99999999]

The script compares the calculated matrix with this reference
using an absolute tolerance of 1e-8 and zero relative tolerance.

This comparison checks agreement with the printed textbook result;
it does not establish the physical accuracy of the model at other
epochs.

## Limitations

- No frame-bias matrix is included.
- Observed corrections to the nutation angles are omitted.
- Polar motion is linearized.
- Earth orientation parameters are supplied by the caller;
  the module does not download or interpolate them.
- Leap-second instants require dedicated time handling.
- Transforming velocity between rotating frames requires an
  additional term involving the time derivative of the
  transformation matrix. Multiplying velocity by U alone is
  insufficient.

## References

Montenbruck and Gill, Satellite Orbits:

- Section 5.3.2: Precession.
- Section 5.3.3: Nutation.
- Section 5.4: Earth rotation and polar motion.
- Exercise 5.1, page 190: Input data and reference matrices.
