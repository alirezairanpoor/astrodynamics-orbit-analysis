# Exercise 4.2 — Fourth-Order Gauss–Jackson

## Objective

Integrate the normalized planar two-body problem using the
fourth-order Gauss–Jackson method.

Compare the final state with a Keplerian reference solution
and count acceleration evaluations.

## Test problem

- Gravitational parameter: mu = 1.
- Semi-major axis: a = 1.
- Eccentricity: e = 0.1.
- Initial time: 0.
- Final time: 20.
- Initial position: [1 - e, 0].
- Initial velocity: [0, sqrt((1 + e) / (1 - e))].

All quantities use normalized units.

State:

    [x, y, vx, vy]

Acceleration:

    acceleration = -position / norm(position)^3

## Implementation

Reusable integrator:

    satellite_orbits.astrodynamics.propagation.integrate_gauss_jackson4

Initialization reuses:

    satellite_orbits.astrodynamics.propagation.integrate_rk4

Three backward RK4 steps generate states at -h, -2h, and -3h.
Accelerations at these times and at t = 0 initialize the
backward differences and the first and second sums.

The explicit Gauss–Jackson formulas then calculate position
and velocity using the stored acceleration history.

The initialization requires:
- 12 acceleration evaluations for the three RK4 steps.
- 4 evaluations for the initial acceleration history.

Each subsequent step requires one new acceleration evaluation:

    total evaluations = number of steps + 16

The implemented interface accepts acceleration depending on
time and position. It does not directly support acceleration
depending on velocity, such as atmospheric Drag.

## Error definition

The reference state is obtained from Kepler's equation.

    error = norm(numerical_state - reference_state)
    digits = -log10(error)

The norm includes both position and velocity components in
normalized units. It is not a position error in meters.

## Results

| Steps | Step size | Evaluations | State error | Digits |
|---:|---:|---:|---:|---:|
| 100 | 0.200000 | 116 | 2.455686e-03 | 2.61 |
| 300 | 0.066667 | 316 | 3.363006e-06 | 5.47 |
| 600 | 0.033333 | 616 | 7.705286e-08 | 7.11 |
| 1000 | 0.020000 | 1016 | 4.539193e-09 | 8.34 |
| 1500 | 0.013333 | 1516 | 5.174494e-10 | 9.29 |
| 2000 | 0.010000 | 2016 | 1.166413e-10 | 9.93 |
| 3000 | 0.006667 | 3016 | 1.499746e-11 | 10.82 |
| 4000 | 0.005000 | 4016 | 3.577889e-12 | 11.45 |

## Verification

- All evaluation counts equal N + 16.
- State errors agree with the textbook within the selected
  comparison tolerances: rtol = 0.002 and atol = 5e-13.

Both checks passed.

Small differences in the final digits can occur because of
floating-point arithmetic.

## Comparison with RK4

For 1000 steps in this test:

- RK4: 4000 evaluations and state error 1.996416e-7.
- Gauss–Jackson: 1016 evaluations and state error 4.539193e-9.

Using acceleration history reduces the number of evaluations
and produces a smaller endpoint error for this particular test.

This does not establish that Gauss–Jackson is preferable for
every orbit or force model.
