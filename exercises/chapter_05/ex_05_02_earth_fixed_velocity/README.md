# Exercise 5.2 — Velocity in the Earth-fixed Frame

## Objective

Transform the position and velocity of GPS satellite PRN 15 from
WGS84(G873), treated as ITRS, to the mean equator and equinox of
J2000, used as the textbook approximation to ICRS.

Unlike position, velocity transformation must account for the
rotation of the terrestrial reference frame.

## Inputs

Epoch:

- GPS: 1999-03-04 00:00:00
- GPS − UTC: 13 s
- UTC: 1999-03-03 23:59:47

Earth orientation parameters supplied by the textbook:

- TT − UTC: 64.184 s
- UT1 − UTC: 0.649232 s
- xp: 0.06740 arcsec
- yp: 0.24173 arcsec

Earth-fixed position [km]:

    [19440.953805, 16881.609273, -6777.115092]

Earth-fixed coordinate velocity [km/s]:

    [-0.8111827456, -0.2573799137, -3.0689508125]

The gravitational parameter used for the orbital diagnostics is
imported from the project's MU_EARTH constant [km^3/s^2].

## Calculation

The forward position and velocity transformations are:

    r_itrs = U @ r_icrs
    v_itrs = U @ v_icrs + U_dot @ r_icrs

The transformation matrix and its approximate derivative are:

    U = Pi @ Theta @ N @ P
    U_dot = Pi @ Theta_dot @ N @ P

As specified in the exercise, precession, nutation, and polar motion
are held constant when computing the derivative. Only Earth rotation
contributes to U_dot.

For the inverse transformation, the implementation solves:

    U @ r_icrs = r_itrs
    U @ v_icrs = v_itrs - U_dot @ r_icrs

Linear solves invert the implemented matrix directly, including
the first-order polar-motion approximation.

## Reusable Functions

Module:

    src/satellite_orbits/coordinates/icrs_itrs_state.py

Functions:

- icrs_to_itrs_matrix_and_rate:
  Computes U and U_dot from the epoch and supplied EOP values.
- icrs_to_itrs_state:
  Transforms J2000 position and velocity to ITRS.
- itrs_to_icrs_state:
  Transforms ITRS position and velocity to J2000.

The matrix calculation reuses the precession, nutation, Earth
rotation, and polar-motion functions from Exercise 5.1.

Units:

- Position: km
- Velocity: km/s
- U: dimensionless
- U_dot: 1/s
- Time offsets: seconds
- Pole coordinates: arcseconds
- Modified Julian Date: days

All state vectors and transformation matrices must refer to the
same physical epoch.

Dependencies: NumPy and PyERFA.

## Results

Position in J2000 [km]:

    [-23830.593391, -9747.073876, -6779.828533]

Velocity in J2000 [km/s]:

    [1.561964406, -1.754345710, -3.068850601]

Osculating-orbit diagnostics:

- Semimajor axis: 26561.013 km
- Eccentricity: 0.0070606
- Inclination: 56.338 deg
- Keplerian period: 11.966728 h

These values describe the expected near-circular GPS orbit,
with an inclination near 56 degrees and a period near 12 hours.

## Verification

The executed solution passed:

- Position comparison with the rounded textbook result:
  absolute tolerance 0.001 km per component.
- Velocity comparison with the rounded textbook result:
  absolute tolerance 1e-6 km/s per component.
- Semimajor-axis, eccentricity, inclination, and period checks.
- J2000-to-ITRS round-trip recovery:
  absolute tolerances of 1e-8 km and 1e-11 km/s per component.

All numerical comparisons use zero relative tolerance.

The round-trip check establishes internal consistency.
The textbook comparison checks agreement with the exercise model;
neither establishes accuracy beyond that model's assumptions.

## Assumptions and Limitations

- WGS84(G873) is treated as ITRS for this exercise.
- J2000 mean-equatorial coordinates approximate ICRS.
- The transformation inherits the IAU 1976/1980 models and
  first-order polar motion from Exercise 5.1.
- U_dot neglects the rates of precession, nutation, and polar motion.
- The supplied EOP values are used as instructed by the exercise.
  Other epochs require appropriate EOP and time-offset values.
- Leap-second instants require dedicated time handling.
- The reported Keplerian period is an osculating two-body diagnostic,
  not a prediction from a perturbed orbit propagation.

## References

Montenbruck and Gill, Satellite Orbits:

- Exercise 5.2, pages 191–192.
- Equation 5.90: Position and velocity transformation.
- Equations 5.91–5.92: Approximate transformation derivative.
