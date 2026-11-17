# Exercise 6.4 — Tropospheric Refraction

## Objective

Compute the daily variation of the elevation of a near-
geostationary satellite seen from the ground station of
Exercise 6.1, and compare that variation with the amplitude of the
tropospheric refraction correction for two air temperatures.

Reference: Montenbruck and Gill, Satellite Orbits,
Exercise 6.4, page 232.

## Physical Background

### The bending angle

A signal entering the atmosphere at zenith distance z0 travels
inside it at a smaller angle z, by Snell's law:

    n sin z = sin z0

Substituting the elevation E = 90 deg - z and expanding in
dE = E - E0 gives equation 6.50:

    dE = (n - 1) / tan E0

The correction vanishes at the zenith, where the signal crosses the
atmospheric layers at normal incidence, and diverges towards the
horizon, where the plane-atmosphere assumption fails. It should not
be trusted below a few degrees of elevation.

### Refractivity

The refractivity N = (n - 1) * 1e6 splits into a dry and a wet
term, equations 6.53 to 6.55:

    N1 = 77.624 * p / T

    N2 = 371900 * e / T^2 - 12.92 * e / T

    e = 6.10 * fh * exp[17.15 * Tc / (234.7 + Tc)]

The two terms respond to temperature in opposite directions. The
dry term falls as 1/T, while the wet term rises steeply because the
saturation vapour pressure grows exponentially: about 12.3 hPa at
10 C against 42.4 hPa at 30 C. The net effect is therefore not
obvious in advance and has to be evaluated.

## Inputs

Satellite:

    Semi-major axis = 42164 km
    Eccentricity    = 0.00296
    Inclination     = 0.05 deg
    RAAN            = 150.7 deg

The argument of perigee and the mean anomaly are not specified by
the exercise and are taken as zero at the epoch.

Epoch: 1997/01/01 0h UTC. Elevation evaluated every 3 h over a day.

Station: 11 deg East, 48 deg North, on the WGS84 ellipsoid.

Meteorology: p = 1024 hPa, fh = 0.7, T1 = 283 K, T2 = 303 K.

## Reusable Functions

New module introduced by this exercise:

    src/satellite_orbits/observations/tropospheric_refraction.py

providing water_vapour_pressure, tropospheric_refractivity and
elevation_refraction_correction.

Also added:

    src/satellite_orbits/coordinates/ecef_to_enz_rotation.py

This is a pure rotation from ECEF to the local ENZ frame, taking
the station-to-satellite vector as input. It is the exact inverse
of enz_to_ecef.

The existing ecef_to_enz performs both the translation and the
rotation, and builds the station position on a sphere of radius
R_EARTH from the geodetic latitude. For a mid-latitude station that
places the station 24.3 km away from its position on the reference
ellipsoid, which at geostationary distance is comparable to the
entire effect this exercise measures. Separating the translation
from the rotation lets the caller build the station position with
geodetic_to_ecef instead. ecef_to_enz is left unchanged, since the
spherical approximation is adequate where it is already used.

## Results

    Elevation at the epoch: 23.340 deg

    T = 283.0 K  ->  N1 = 280.873, N2 =  39.173, N = 320.046
    T = 303.0 K  ->  N1 = 262.333, N2 = 118.514, N = 380.848

    UTC        E - E0     dE(T1)     dE(T2)
    hh:mm:ss    [deg]      [deg]      [deg]
    00:00:00.0     0.000      0.042      0.051
    03:00:00.0    -0.087      0.043      0.051
    06:00:00.0    -0.107      0.043      0.051
    09:00:00.0    -0.050      0.043      0.051
    12:00:00.0     0.051      0.042      0.050
    15:00:00.0     0.136      0.042      0.050
    18:00:00.0     0.155      0.042      0.050
    21:00:00.0     0.098      0.042      0.050
    24:00:00.0    -0.003      0.043      0.051

## Engineering Conclusion

The refraction correction, 0.042 to 0.051 deg, is of the same order
as the daily elevation variation, but it behaves very differently.
It varies by only 0.0005 deg over the day for fixed weather, so it
acts as an almost constant bias. The orbital motion, by contrast,
is a genuine signal.

The 20 K temperature difference changes the correction by
0.0081 deg, matching the 0.008 deg the textbook quotes. At
geostationary distance 0.001 deg corresponds to about 700 m of
along-track or cross-track position, so unmodelled weather
variations map directly into metres of apparent position error.

## Two Discrepancies with the Textbook

Both were investigated rather than tuned away.

### The two temperature columns appear interchanged

The textbook prints dE(T1) = 0.051 deg and dE(T2) = 0.042 deg for
T1 = 283 K and T2 = 303 K. This implementation produces exactly
those two values but assigned to the opposite temperatures.

The physics decides the direction. Between 283 K and 303 K the wet
refractivity rises from 39.2 to 118.5, a factor of three, while the
dry term only falls from 280.9 to 262.3. The exponential growth of
the saturation vapour pressure dominates, so warmer and more humid
air bends the signal more. The 303 K case must give the larger
correction.

Supporting evidence: the difference between the two values is
0.0081 deg, and the textbook itself states in the surrounding text
that the 20 K difference changes the elevation by about 0.008 deg.
The magnitudes are therefore correct and only the column headings
disagree.

### The daily elevation variation is larger than reported

The textbook reports a peak-to-peak variation of about 0.058 deg
with the maximum 6 h past the epoch. This implementation gives
0.263 deg with the maximum at 18 h.

The cause was isolated by switching the orbital parameters off one
at a time:

    as given      p2p = 0.2627   sub-lon spread 0.6773
    e = 0         p2p = 0.0932   sub-lon spread 0.0022
    i = 0         p2p = 0.3553   sub-lon spread 0.6773
    e = 0, i = 0  p2p = 0.0011   sub-lon spread 0.0022

The e = 0, i = 0 case gives a flat elevation, which confirms that
the underlying model carries no spurious drift: the semi-major axis
of 42164 km is geostationary to within 0.17 km, worth 0.0022 deg of
longitude drift per day.

The dominant term is the eccentric libration. A satellite with
eccentricity e oscillates in longitude with an amplitude of 2e
radians, here 0.339 deg, producing the observed 0.68 deg
peak-to-peak excursion of the sub-satellite longitude. With the
stated eccentricity this motion is unavoidable.

A control case with e = 0 recovers the textbook elevation at the
epoch almost exactly, 23.366 deg against 23.363 deg, which is
strong evidence that the textbook modelled the daily motion without
the eccentric libration. The textbook also attributes its 0.060 deg
variation to the inclination, not to the eccentricity.

This explanation is incomplete. The circular control case still
gives a peak-to-peak variation of 0.0932 deg against the textbook's
0.058 deg, a factor of 1.6 that remains unexplained. Sweeping the
unspecified argument of perigee and mean anomaly over a full grid
produced no combination reproducing both the textbook elevation and
its variation, so the residual is not a phase convention. The
discrepancy is recorded here as open.

The engineering conclusion of the exercise is unaffected, and in
fact strengthened: the refraction correction is a near-constant
bias while the orbital motion is a real signal, whichever of the
two variation amplitudes is adopted.

## Scope and Limitations

- A plane atmosphere of constant refractivity is assumed. The
  more accurate Hopfield-Goad model of equations 6.59 to 6.65 is
  not implemented. Against Table 6.4 this simple model agrees to
  within a few percent above 20 deg elevation but is 8 percent
  high at 10 deg.
- Ray-path curvature is neglected.
- Azimuth is not affected by refraction and is not corrected.
- No ionospheric correction. The troposphere is non-dispersive for
  radio waves, so the result does not depend on carrier frequency.
