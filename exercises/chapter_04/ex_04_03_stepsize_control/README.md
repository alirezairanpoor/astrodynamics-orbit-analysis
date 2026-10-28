# Exercise 4.3 — Stepsize Control for Eccentric Orbits

## Objective

Integrate the normalized two-body problem with eccentricity
e = 0.9 using the Shampine–Gordon DE multistep method.

Monitor the accepted step size as a function of time and
distance from the central body.

## Test problem

- Gravitational parameter: mu = 1.
- Semi-major axis: a = 1.
- Eccentricity: e = 0.9.
- Initial time: 0.
- Final time: 20.
- Initial position: [1 - e, 0].
- Initial velocity: [0, sqrt((1 + e) / (1 - e))].

All quantities use normalized units.

The orbit starts at pericenter:

- Pericenter distance: 0.1.
- Apocenter distance: 1.9.
- Orbital period: 2*pi.

State:

    [x, y, vx, vy]

Acceleration:

    acceleration = -position / norm(position)^3

## Implementation

Reusable integrator:

    satellite_orbits.astrodynamics.propagation.integrate_de

The implementation adapts the Shampine–Gordon STEP algorithm,
using variable-step, variable-order Adams predictor–corrector
formulas.

The method:
- Predicts the next state from stored derivative information.
- Estimates local error.
- Rejects and retries an unsuitable step.
- Corrects an accepted prediction.
- Selects the next step size and order.

Orders range from 1 to 12.

The driver supports forward integration with scalar absolute
tolerance and an unweighted Euclidean error norm. It returns
accepted step endpoints without dense interpolation and limits
the last step to the requested final time.

It does not reproduce the complete historical DE/ODE driver API.

Dependencies: NumPy and Matplotlib.

## Reference solution and checks

Kepler's equation provides independent reference states at
the accepted integration times.

The reported state error is:

    error = norm(numerical_state - reference_state)

This includes normalized position and velocity components;
it is not a position error in meters.

The exercise uses an absolute local tolerance of 1e-8.
An additional run with tolerance 1e-10 checks sensitivity
to the tolerance setting.

## Results

| Quantity | Absolute tolerance 1e-8 | Absolute tolerance 1e-10 |
|---|---:|---:|
| Accepted steps | 792 | 1116 |
| Derivative evaluations | 1606 | 2254 |
| Orders used | 1–12 | 1–12 |
| Smallest accepted step | 1.249406954e-6 | 1.249406954e-7 |
| Largest accepted step | 1.966387035e-1 | 1.353334858e-1 |
| Final state error | 2.395785422e-6 | 3.330893851e-9 |
| Maximum sampled state error | 2.447190431e-4 | 3.994382225e-7 |

The tighter tolerance reduced both the endpoint error and
the maximum error sampled on each run's accepted grid.

For the 1e-8 run, excluding the first orbit and the final
endpoint-limited step:

- Median step for r < 0.2: 1.896592433e-3.
- Median step for r > 1.8: 1.268897161e-1.
- Ratio of these medians: 66.90.

The smallest step reported for the full run includes startup
and should not be interpreted as a typical pericenter step.

## Plots and output files

The program produces:

    outputs/accepted_steps.csv
    outputs/stepsize_control.png

The CSV contains:

    t_start, t_end, accepted_step, r_start, order

The two plots show:
1. Accepted step size against the starting time of each step.
2. Accepted step size against the starting distance of each step.

These are actual integration steps, not uniformly requested
output intervals.

## Interpretation

Near pericenter, rapid motion requires smaller integration steps.
Near apocenter, slower motion permits larger steps.

The step size is selected from estimated numerical error.
No explicit step-size-versus-distance law is imposed.

The step size can remain constant over several accepted steps,
which explains the plateaus in the time plot.

Different step sizes can occur at the same distance because
the current order and integration history also affect selection.

This behavior is consistent with the qualitative result in
Figure 4.10 of the textbook. Exact reproduction of its step
sequence is not claimed.

## Limitations

The local tolerance is not a bound on the accumulated global
state error.

Adaptive step and order decisions may vary with floating-point
roundoff and implementation details.

This test evaluates numerical integration of an idealized
two-body model. It does not quantify real-orbit uncertainty
or force-model error.

The absolute-only interface requires appropriately scaled
state components before use in other problems.

## Algorithm reference

Lawrence F. Shampine and Marilyn K. Gordon,
Computer Solution of Ordinary Differential Equations (1975).

Python adaptation of the STEP algorithm distributed in the
ODE library:

https://people.math.sc.edu/Burkardt/f77_src/ode/ode.f

The integrate_de.py adaptation is identified as
LGPL-3.0-or-later in its source header.
