# Exercise 3.3 — Accelerations

## Objective

Find the altitudes where simplified perturbation acceleration
magnitudes are equal:

- J20 versus Moon and Sun
- J22 versus Moon and Sun
- Atmospheric drag versus solar radiation pressure

Equality of magnitudes does not imply cancellation of the vectors.

## Implementation

Reusable functions:

- gravity_third_body_balance_radius
- drag_srp_balance_altitude
- interpolate_density

The gravitational balance radius is calculated algebraically
using the simplified equations supplied in the exercise.

The drag/SRP balance altitude is found using bisection.
Density is interpolated exponentially between tabulated heights.

## Units and Coefficients

- Distances and altitudes: km
- Gravitational parameters: km^3/s^2
- Gravity coefficients: normalized, from Table 3.3
- Density: kg/m^3
- Solar radiation pressure: N/m^2

The drag/SRP calculation converts distances and the gravitational
parameter to SI internally.

Table 3.8 densities are converted using:
1 g/km^3 = 1e-12 kg/m^3.

## Gravitational Results

Balance altitudes:

- J20 / Moon: 35984.527 km
- J20 / Sun: 43118.422 km
- J22 / Moon: 8750.724 km
- J22 / Sun: 11298.432 km

These are close to the book's reported values:
35983, 43117, 8750, and 11298 km, respectively.

Exact numerical agreement was not obtained with the explicitly
listed constants and coefficients used here.

## Drag / SRP Assumptions

- Circular orbit
- Full sunlight
- Solar pressure at 1 AU: 4.560e-6 N/m^2
- C_D = 2.3
- C_R = 1.3
- Atmospheric rotation and winds neglected
- Same effective area for drag and SRP

Under these assumptions, the area-to-mass ratio cancels from
the balance equation.

## Atmospheric Model

The minimum and maximum Harris–Priester density envelopes
from Table 3.8 are evaluated separately.

Both envelopes correspond to mean solar activity.
They do not represent minimum and maximum solar activity.

Only the 500–1000 km table segment is used.
The full Sun-direction-dependent Harris–Priester model
is not implemented in this exercise.

## Drag / SRP Results

- Minimum-density envelope: 592.679 km
- Maximum-density envelope: 785.568 km

Higher density allows drag to remain dominant to a greater altitude.

The book reports approximately 752 km. That value lies between
these results, but this does not constitute reproduction of
the book's specific atmospheric calculation.


