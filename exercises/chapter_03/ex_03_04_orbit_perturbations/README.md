# Exercise 3.4 — Orbit Perturbations

Based on Chapter 3 of Montenbruck & Gill, *Satellite Orbits*.

## Objective

Evaluate how simplifying the force model changes the propagated
satellite position.

Each restricted model starts from the same initial state as the
reference model. The reported quantity is the maximum sampled
position difference over the specified propagation interval.

The reference trajectory is a numerical model, not an observed
true orbit.

## Initial conditions

Epoch: 1999-03-01 00:00:00 UTC.

| Parameter | LEO | GEO |
|---|---:|---:|
| Semi-major axis [km] | 7178.0 | 42166.0 |
| Eccentricity | 0.001 | 0.0004 |
| Inclination [deg] | 98.57 | 0.02 |
| RAAN [deg] | 0.0 | 0.0 |
| AOP [deg] | 0.0 | 0.0 |
| Mean anomaly [deg] | 0.0 | 0.0 |
| Effective area [m²] | 5.0 | 10.0 |
| Mass [kg] | 1000.0 | 1000.0 |
| Radiation pressure coefficient | 1.3 | 1.3 |
| Drag coefficient | 2.3 | 2.3 |

LEO propagation intervals:
- One initial Keplerian period: approximately 6052.240 s.
- One day.

GEO propagation intervals:
- One day.
- Two days.

## Force model

The reference model includes:

- JGM-3 Earth gravity through degree/order 20/20.
- Differential third-body acceleration from the Sun and Moon.
- Direct solar radiation pressure with Earth shadowing.
- Atmospheric drag for LEO.

Earth gravity coefficients are loaded from the implemented JGM-3
table and converted to the unnormalized convention required by
the Cunningham V/W evaluation.

The gravity function receives an ECEF position and returns ECEF
acceleration, including central gravity. This acceleration is
rotated to ECI before combining it with the other accelerations.

Sun and Moon positions are obtained from the implemented
low-precision analytical ephemerides.

The shadow function uses the apparent overlap of the Earth and
Sun disks, including partial illumination in penumbra.

LEO density uses the Harris–Priester mean-solar-activity table,
with exponent 6 and an apex offset of 30 degrees. Altitude is
measured above the reference ellipsoid.

The atmosphere is assumed to co-rotate with Earth; winds are
neglected.

GEO drag is omitted by assumption. Its reported zero is not the
result of evaluating atmospheric density at GEO altitude.

## Restricted models

Gravity truncation cases:

- 2/0
- 2/2
- 4/4
- 10/10

These retain all available coefficients within the specified
degree/order limits, including central gravity. Other forces
remain enabled.

Individual force omission cases:

- No Sun gravity
- No Moon gravity
- No SRP
- No Drag

These retain Earth gravity at degree/order 20/20.

Disabling Sun gravity does not disable SRP.

## Numerical propagation

The state contains ECI position and velocity:

    state = [x, y, z, vx, vy, vz]

Its derivative contains velocity and total acceleration:

    derivative = [vx, vy, vz, ax, ay, az]

Propagation uses classical fixed-step RK4:

- LEO step: 5 s.
- GEO step: 30 s.

For each satellite, the reference orbit is also propagated with
half the step size. Positions are compared at common epochs.

The selected acceptance threshold for this comparison is 0.1 m.

| Reference trajectory | Maximum sampled step-halving difference [m] |
|---|---:|
| LEO, one day | 0.014618 |
| GEO, two days | 0.000214 |

Both reference runs satisfy the selected threshold.

This checks sensitivity to the integration step. It does not
establish the physical accuracy of the force model or constitute
an independent validation of every restricted trajectory.

## Results

Maximum sampled position differences relative to the reference
trajectory, in meters:

| Restricted model | LEO 1 rev | LEO 1 day | GEO 1 day | GEO 2 days |
|---|---:|---:|---:|---:|
| 2/0 | 599.921 | 5028.362 | 671.153 | 2534.752 |
| 2/2 | 223.804 | 3040.385 | 2.295 | 9.925 |
| 4/4 | 147.704 | 1927.199 | 0.137 | 0.556 |
| 10/10 | 22.915 | 458.885 | 0.000 | 0.000 |
| No Sun | 2.865 | 33.787 | 3142.823 | 4834.574 |
| No Moon | 6.004 | 65.917 | 5080.416 | 5439.380 |
| No SRP | 0.697 | 9.761 | 408.505 | 816.080 |
| No Drag | 0.463 | 104.542 | 0.000 | 0.000 |

Values displayed as 0.000 may be below the printed precision.
The GEO No Drag entries are zero by construction.

Maxima are evaluated on the integration output grid. For the LEO
one-revolution endpoint, which falls between grid samples, the
position-difference vector is linearly interpolated to include
that endpoint approximately.

## Comparison with the textbook

Textbook reference values, in meters:

| Restricted model | LEO 1 rev | LEO 1 day | GEO 1 day | GEO 2 days |
|---|---:|---:|---:|---:|
| 2/0 | 600 | 5028 | 671 | 2534 |
| 2/2 | 224 | 3038 | 2 | 10 |
| 4/4 | 148 | 1925 | 0 | 0 |
| 10/10 | 23 | 459 | 0 | 0 |
| No Sun | 3 | 34 | 3143 | 4834 |
| No Moon | 6 | 66 | 5080 | 5438 |
| No SRP | 1 | 14 | 415 | 830 |
| No Drag | 1 | 105 | 0 | 0 |

The four LEO one-revolution gravity-truncation results reproduce
the textbook values when rounded to meters.

Other gravity and third-body results are close, but agreement
is not exact. For example, the LEO one-day 4/4 result is
1927.199 m versus 1925 m in the textbook.

The LEO one-day Drag result is 104.542 m versus 105 m.

SRP shows a more noticeable discrepancy:

- LEO one day: 9.761 m versus 14 m.
- GEO one day: 408.505 m versus 415 m.
- GEO two days: 816.080 m versus 830 m.

The SRP and apparent-disk shadow implementations were reviewed
against equations 3.75 and 3.85–3.94. No formula-level error was
identified in that review.

The cause of the remaining SRP discrepancy has not been
established. It is retained as an unresolved validation item;
parameters were not adjusted to force agreement.

## Assumptions and limitations

- ECI/ECEF conversion uses Earth rotation only; precession,
  nutation, and polar motion are omitted.
- UT1 is approximated by UTC.
- TT − UTC is fixed at 64.184 s for the exercise dates.
- Analytical Sun and Moon ephemerides have limited precision.
- JGM-3 coefficients use the precision of the implemented
  textbook table.
- SRP uses constant effective area and radiation coefficient.
- Shadow geometry assumes spherical Earth and Sun and uniform
  solar disk brightness.
- Harris–Priester density represents mean solar activity.
- GEO atmospheric drag is neglected.
- Reported maxima are sampled estimates, not continuously
  optimized extrema.
- Exact reproduction of the complete textbook table has not
  been achieved.

## Physical interpretation

For these specific orbits and propagation intervals:

- Higher-degree Earth gravity terms strongly affect LEO
  position prediction.
- Sun and Moon perturbations become especially important
  for GEO.
- Drag accumulates into a significant LEO position difference
  over one day.
- The force-model detail required depends on orbit altitude,
  propagation duration, and the desired position accuracy.

## Status

All requested propagation cases have been implemented and run.
Reference step-halving checks passed the selected criterion.
Remaining textbook discrepancies are documented, including
the unresolved SRP difference.
