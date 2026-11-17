# Exercise 6.3 — User Clock Error from GPS Pseudorange

## Objective

Determine the clock error of the Rogue GPS receiver at Goldstone
from dual-frequency pseudoranges to PRN 1, using a precise IGS
ephemeris and a rigorous light-time solution.

Reference: Montenbruck and Gill, Satellite Orbits,
Exercise 6.3, pages 230-232.

## Physical Background

A pseudorange is not a range. It is the difference between the
reception time read on the receiver clock and the transmission time
read on the satellite clock, scaled by the speed of light. Both
clocks are offset from GPS system time, so

    rho = |r_r(t) - r_t(t - tau)| + c * (dt_r - dt_t)

Everything in this equation is known except the receiver clock
offset dt_r, which is what the exercise asks for:

    c * dt_r = rho - d + c * dt_t

The satellite clock offset dt_t is about 40 microseconds, worth
roughly 12000 m, so it is the dominant term in the answer and must
be taken from the precise ephemeris rather than neglected.

Three modelling steps stand between the raw observation and d.

### Ionospheric delay

The ionospheric group delay scales as 1 / frequency^2, so two
frequencies give two equations in one unknown electron content and
the delay can be eliminated rather than modelled:

    rho = [rho_2 - (f1/f2)^2 * rho_1] / [1 - (f1/f2)^2]

For the GPS L1 and L2 frequencies the weights are +2.5457 on L1 and
-1.5457 on L2, matching equation 6.72. Being large and of opposite
sign, they amplify the measurement noise by about a factor of three.
That is an acceptable trade against a systematic error of tens of
metres.

### Light time

The signal takes about 70 ms to reach the receiver, during which
the satellite moves roughly 270 m. Only the along-line-of-sight
component matters, and it amounts to 24 m here.

### Frame consistency

This is the step that distinguishes this exercise from Exercise 6.1.
The ephemeris and the receiver coordinates are both Earth-fixed, so
no inertial transformation is wanted. But the satellite position at
emission refers to the Earth-fixed frame of that epoch, while the
receiver refers to the frame at reception, and the Earth has turned
by omega * tau between them. Equation 6.77 removes the mismatch:

    c * tau = |Rz(omega * tau) * r_WGS(t - tau) - R_WGS|

Rz here is the coordinate rotation, which expresses a fixed vector
in a rotated frame. It is the transpose of the matrix that rotates
a vector within a fixed frame. Using the wrong one of the two
produces an error of twice the correction, about 9 m at the first
epoch, which was how the sign was identified during development.

## Inputs

PRN 1 precise ephemeris from the IGS, 1998/02/19, thirteen epochs
at 15 min intervals from 08:00 to 11:00 GPS time, with satellite
clock offsets.

Pseudoranges at six epochs at 30 min intervals from 08:30 to 11:00:
P code at L2 and C/A code at L1. Observations are restricted to
epochs with precise ephemeris data to avoid Selective Availability.

Receiver position at Goldstone, WGS84 [km]:

    X = -2353.614128
    Y = -4641.385447
    Z = +3676.976501

All lengths are converted to km, the unit used throughout this
project.

## Reusable Functions

New modules introduced by this exercise:

    src/satellite_orbits/ephemerides/
        interpolate_tabulated_position.py
    src/satellite_orbits/observations/
        ionosphere_free_pseudorange.py
        light_time_earth_fixed.py

interpolate_tabulated_position evaluates a Lagrange polynomial by
Neville's algorithm over a node window centred on the requested
epoch. It sits beside chebyshev_position, since both interpolate a
tabulated ephemeris; it is not a propagator and does not belong
with the integrators.

light_time_earth_fixed is unit-agnostic. The position vectors, the
speed of light and the tolerance must all use the same length unit.
An early version of this exercise passed metre positions with the
project's km speed of light, which inflated the travel time by a
factor of 1000 and was caught by a sanity check on tau.

## Results

    GPS       Iono-free      Distance   SatClock   ClockError
    hh:mm:ss        [m]           [m]        [m]          [m]
    08:30:00.0   21096582.633  21090908.389  12000.301     17674.55
    09:00:00.0   20519970.005  20514343.594  12049.746     17676.16
    09:30:00.0   20282712.756  20277032.664  11995.499     17675.59
    10:00:00.0   20375843.885  20370215.515  12046.068     17674.44
    10:30:00.0   20751684.441  20746014.482  12003.550     17673.51
    11:00:00.0   21340060.814  21334385.767  11999.341     17674.39

    Mean clock bias: 17674.77 m
    Maximum spread:  2.65 m

The clock drift is negligible over the 2.5 h interval, so the mean
of the six residuals gives the clock bias. The spread of 2.65 m
reflects the accuracy of the ionosphere-free pseudoranges, not a
variation of the receiver clock.

## Two Findings Worth Recording

### Interpolation order is not the limiting factor

Repeating the computation at orders 5, 7, 9, 11 and 12 gives clock
errors identical to the second decimal in metres. With a 15 min
grid, even a moderate order already interpolates the GPS orbit far
below the accuracy of the pseudoranges. Order 9 is kept as a
conventional choice, not because the result depends on it.

### The relativistic clock term must not be applied

Equation 6.44 gives a periodic relativistic correction

    dt_rel = -2 * (r . v) / c^2

worth -2.6 to +0.1 m over this interval. Applying it degrades the
agreement with the textbook rather than improving it, because the
IGS precise clock values already include it. The correction is
needed when the satellite clock offset is reconstructed from the
broadcast ephemeris coefficients af0, af1 and af2, which is not the
case here. It is therefore deliberately omitted.

## Verification

All six clock errors match the textbook with an absolute tolerance
of 0.5 m, the mean bias matches to better than 0.01 m, and the
residual spread stays below the 3 m bound implied by the textbook
discussion.

Unlike Exercises 6.1 and 6.2, there is no systematic reference
frame residual here, because no inertial transformation is involved
and the sidereal-time model is never used.

## Scope and Limitations

- Only the first-order ionospheric term is removed. Higher-order
  terms remain, at the millimetre to centimetre level.
- No tropospheric delay. For a satellite well above the horizon at
  Goldstone this is of the order of a few metres and is absorbed
  into the clock bias estimate.
- Antenna phase centre offsets, site displacements and multipath
  are not modelled.
- The receiver clock is treated as constant over the interval. No
  drift term is estimated.
