# Exercise 4.1 — Fourth-Order Runge–Kutta

## Objective

Apply classical RK4 to the normalized planar two-body problem
and compare the final state with a Keplerian reference solution.

Evaluate accuracy and computational cost for different fixed
step sizes.

## Test problem

- Gravitational parameter: mu = 1.
- Semi-major axis: a = 1.
- Eccentricity: e = 0.1.
- Initial time: 0.
- Final time: 20.
- Initial position: [1 - e, 0].
- Initial velocity: [0, sqrt((1 + e) / (1 - e))].

All quantities use normalized units.

The state is [x, y, vx, vy], and acceleration is:

    acceleration = -position / norm(position)^3

## Implementation

The exercise reuses:

    satellite_orbits.astrodynamics.propagation.integrate_rk4

The reference state is computed by solving Kepler's equation
at the final time and evaluating position and velocity.

The reported error is the Euclidean norm of the complete
normalized state difference:

    error = norm(numerical_state - reference_state)

The accuracy indicator is:

    digits = -log10(error)

Derivative evaluations are counted during execution.
Classical RK4 requires four evaluations per step.

## Results

| Steps | Step size | Evaluations | State error | Digits |
|---:|---:|---:|---:|---:|
| 50 | 0.400000 | 200 | 1.953483e-01 | 0.71 |
| 100 | 0.200000 | 400 | 6.663267e-03 | 2.18 |
| 250 | 0.080000 | 1000 | 9.051025e-05 | 4.04 |
| 500 | 0.040000 | 2000 | 4.012045e-06 | 5.40 |
| 750 | 0.026667 | 3000 | 6.847233e-07 | 6.16 |
| 1000 | 0.020000 | 4000 | 1.996416e-07 | 6.70 |
| 1500 | 0.013333 | 6000 | 3.608269e-08 | 7.44 |
| 2000 | 0.010000 | 8000 | 1.088728e-08 | 7.96 |

## Verification

- All evaluation counts equal four times the number of steps.
- All final-state errors agree with the rounded textbook values
  within the selected relative tolerance of 0.001.

Both checks passed.

## Interpretation

Smaller steps reduce numerical error at increased computational cost.

For sufficiently small steps, the global error of classical RK4
scales approximately as h^4. Halving the step therefore approaches
a sixteenfold error reduction.

The observed reduction from 1000 to 2000 steps is approximately 18.3.

This exercise measures numerical integration error against a
Keplerian reference. It does not measure force-model error or
uncertainty in a real satellite orbit.
