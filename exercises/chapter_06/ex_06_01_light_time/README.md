# Exercise 6.1 — Light-Time Iteration

## Objective

Solve the downleg and upleg light-time equations for two-way range
measurements of a polar satellite tracked from a ground station in
central Europe, and show how many fixed-point iterations each leg
requires to reach millimetre accuracy.

Reference: Montenbruck and Gill, Satellite Orbits,
Exercise 6.1, page 229.

## Physical Background

A range measurement is not an instantaneous geometric distance. The
signal needs a finite time to travel, so the satellite is no longer
where it was when the signal left it. The downleg light-time
equation states this self-referentially:

    c * tau_d = |r(t - tau_d) - R(t)|

The unknown travel time appears on both sides, which is why a
fixed-point iteration is used. Starting from tau_d = 0 gives the
geometric distance, and each further pass corrects the satellite
position for the motion during the signal travel.

For the upleg the roles are reversed. The transponding epoch is
already fixed by the downleg solution, so the satellite position is
held constant and the station is traced back along the Earth
rotation:

    c * tau_u = |r(t - tau_d) - R(t - tau_d - tau_u)|

The upleg correction is much smaller than the downleg correction
because the inertial station velocity, about 0.31 km/s at this
latitude, is far below the relative satellite velocity.

## Inputs

Satellite, with the ascending node passage as the element epoch:

    Semi-major axis   = 6378.137 + 960 km
    Eccentricity      = 0
    Inclination       = 97 deg
    RAAN              = 130.7 deg
    Argument of perigee = 0
    Mean anomaly at epoch = 0

Epoch of the ascending node passage: 1997/01/01 0h UTC, MJD 50449.0.

Ground station:

    East longitude = 11 deg
    Latitude       = 48 deg
    Height         = 0 (not specified by the exercise)

Evaluation times: 6, 15, and 24 minutes past the ascending node,
which fall near the start, the culmination, and the end of the pass.

## Reusable Functions

New modules introduced by this exercise:

    src/satellite_orbits/coordinates/geodetic_to_ecef.py
    src/satellite_orbits/observations/ground_station_state.py
    src/satellite_orbits/observations/light_time_downleg.py
    src/satellite_orbits/observations/light_time_upleg.py
    src/satellite_orbits/observations/two_way_range.py
    src/satellite_orbits/astrodynamics/keplerian/
        propagate_keplerian_state.py

New constants: C_LIGHT, OMEGA_EARTH, SECONDS_PER_DAY.

The two light-time functions take the moving endpoint as a callable
of time and the fixed endpoint as a vector. This keeps them
independent of how the trajectory is produced, so the same functions
serve a Keplerian orbit here and an interpolated GPS ephemeris in
Exercise 6.3. Both return the travel time together with the signal
path length after each iteration, which is what the exercise asks to
tabulate.

OMEGA_EARTH is the time derivative of Greenwich Mean Sidereal Time,
not 2*pi/86400. The sidereal day is about 236 s shorter than the
solar day, and using the solar value would introduce a 0.27 percent
error in the station velocity, worth about 0.8 m/s.

## Results

    UTC        Distance     Down I-1  Down I-2   Up I-1      Range
    hh:mm:ss        [m]          [m]      [mm]      [m]        [m]
    00:06:00.0    3644876.6       75.8       1.6      1.8    3644954.3
    00:15:00.0    1167048.1       -6.4       0.0      0.5    1167042.2
    00:24:00.0    4110979.3      -87.0       1.8     -0.7    4110891.6

The downleg needs two iterations to reach the millimetre level; the
upleg is already converged after one.

The downleg correction changes sign between the first and the last
evaluation because it is proportional to the range rate: the
satellite approaches the station before culmination and recedes
afterwards. Near culmination the range rate vanishes and so does the
correction, which is why the 15-minute row shows a second-iteration
correction of zero at the printed resolution.

The upleg first-iteration corrections are small and also change
sign, following the component of the station velocity along the line
of sight.

## Verification

The light-time corrections are compared against the textbook with an
absolute tolerance of 0.05 in the printed units, and all three
columns pass at all three epochs.

The absolute distances and two-way ranges carry a systematic
residual with respect to the textbook:

    Geometric distance: [-1.97, -1.78, +0.60] m
    Two-way range:      [-1.94, -1.75, +0.64] m

This residual was investigated rather than tuned away. Applying a
single constant offset to the Greenwich Hour Angle and scanning it
reduces the RMS of the three distance residuals from 1.57 m to
0.047 m at an offset of about -0.195 arcsec. A single rotation
therefore explains all three residuals simultaneously, including the
sign reversal of the third one, which a scale or position error
could not do.

An offset of 0.2 arcsec corresponds to roughly 14 ms of UT1 and is
the expected size of the effects the exercise instructs to omit: the
equation of the equinoxes between GMST and GAST, polar motion at the
1997 epoch, and the particular GMST formulation adopted. No offset
is applied in the code, since a correction fitted to one exercise
would corrupt the reusable modules.

The light-time corrections are immune to this residual because they
are differences of two evaluations a few milliseconds apart, so a
constant rotation cancels. This is why they match the textbook
exactly while the absolute distances do not. The absolute values are
only checked against a 2.5 m bound.

## Scope and Limitations

- Keplerian motion only. No precession, nutation, polar motion, or
  UT1 corrections, as specified by the exercise.
- The station is placed on the WGS84 ellipsoid because no height is
  given.
- Aberration is not modelled. For this geometry it amounts to about
  0.6 arcsec and does not affect the range.
- No tropospheric or ionospheric path delay; those enter in
  Exercise 6.4.
