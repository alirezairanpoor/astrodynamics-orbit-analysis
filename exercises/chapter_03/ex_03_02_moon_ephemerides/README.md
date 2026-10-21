# Exercise 3.2 — Moon Ephemerides

## Objective

Compute the Moon's geocentric position using:
1. The low-precision analytical equations in Section 3.3.2.
2. The DE405 Chebyshev coefficients supplied in Exercise 3.2.

Compare the resulting position vectors.

## Setup

- Epochs: 2006-03-14 through 2006-03-18, at midnight TT
- Sampling interval: one day
- Position units: km
- Origin: Earth center
- Reference frame: equatorial J2000 / EME2000

## Reusable Functions

`moon_position_low_precision(mjd_tt)`

Computes the approximate Moon position using equations 3.47–3.51.

`chebyshev_position(mjd, mjd_start, mjd_end, coefficients)`

Evaluates a Cartesian position from supplied Chebyshev coefficients.
The function preserves the origin, reference frame, and position
units of the coefficient dataset.

## Chebyshev Evaluation

The exercise provides 13 coefficients per coordinate:
39 coefficients in total.

The four-day interval is mapped to [-1, +1].
The Chebyshev polynomials are evaluated recursively, and their
weighted sums give the three position components.

No coefficient fitting or numerical orbit integration is performed.

These coefficients must not be used outside their validity interval.

## Validation

All five Chebyshev positions agree with the book's printed
coordinates within 0.001 km per component.

The analytical implementation reproduces our independently executed
results, but differs from the book's analytical table by approximately
2 km in position. The cause remains unresolved.

## Position Differences

Euclidean norm of analytical position minus Chebyshev position:

- 2006-03-14: 170.461 km
- 2006-03-15: 127.107 km
- 2006-03-16: 71.539 km
- 2006-03-17: 39.561 km
- 2006-03-18: 71.163 km

Maximum over the five sampled epochs: 170.461 km.

The book reports approximately 169 km using its tabulated
analytical results.

## Conclusion

The supplied DE405 Chebyshev representation reproduces the book's
reference positions. The low-precision analytical model differs
from it by tens to hundreds of kilometers over the sampled epochs.

Implementation and comparison are complete.
Exact agreement with the book's analytical table remains unresolved.
