# Exercise 3.1 — Gravity Field

## Objective

Measure the CPU time required by the Cunningham algorithm for
gravity models with equal maximum degree and order from 1 to 20.
Compare the measured growth with quadratic scaling.

## Implementation

Reusable function:
`gravity_acceleration_cunningham`

Location:
`src/satellite_orbits/astrodynamics/perturbations/gravity_acceleration_cunningham.py`

The function computes the V/W recurrence terms and sums the
gravitational acceleration contributions.

Inputs:
- ECEF position [km]
- Unnormalized gravity coefficients C and S
- Maximum degree and order
- Gravitational parameter [km^3/s^2]
- Model reference radius [km]

Output:
- ECEF gravitational acceleration [km/s^2]

Central gravity is included when C[0, 0] = 1.

## Timing Setup

- Fixed ECEF position: [7000, 1200, 2300] km
- Maximum degree = maximum order
- Model sizes: 1 through 20
- 10,000 evaluations per model size
- CPU time measured with time.process_time()
- Printing excluded from the timed interval

Coefficients:
- C[0, 0] = 1
- C[2, 0] = -J2_EARTH
- All other coefficients are zero

This is a central-gravity + J2 dataset, not a complete degree-20
Earth gravity model. The implementation evaluates all requested
terms without skipping zero coefficients, making this setup
suitable for measuring its computational scaling.

## Validation

The degree-zero result was compared with:

a = -mu * r / |r|^3

The central-gravity check passed.

This check alone does not validate all higher-degree terms.

## Results

Recorded CPU times for 10,000 evaluations:

- 4x4 model: 1.1406 s
- 10x10 model: 3.5156 s
- 20x20 model: 12.8281 s

Measured ratio:

t(10) / t(4) = 3.082

Pure quadratic ratio:

(10 / 4)^2 = 6.25

A quadratic fit to the measured times gave:

t(N) = 0.03239455*N^2 - 0.03204780*N + 0.5410088

Here, t is in seconds for 10,000 evaluations.

## Interpretation

The timing trend is consistent with approximately quadratic
computational growth. The number of evaluated harmonic terms
increases approximately as N^2.

Fixed overhead and lower-order computational costs mean that
timing ratios need not equal the corresponding pure N^2 ratios.

The measured 10x10-to-4x4 ratio is close to the approximate
factor of three reported in the book. Absolute timings depend
on hardware, software, and implementation.
