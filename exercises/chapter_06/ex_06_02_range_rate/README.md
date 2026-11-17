# Exercise 6.2 — Range Rate Modeling

## Objective

Compute average two-way range rate measurements for the tracking
scenario of Exercise 6.1 with a Doppler count interval of 1 s, and
show that the measurement is well approximated by the instantaneous
range rate at the middle of the count interval.

Reference: Montenbruck and Gill, Satellite Orbits,
Exercise 6.2, page 229.

## Physical Background

A Doppler system does not measure an instantaneous frequency shift.
It counts accumulated cycles over a count interval tc, and that
count is equivalent to the change of the two-way range across the
interval:

    average_range_rate(t) = [rho(t) - rho(t - tc)] / tc

Each of the two range values requires a full downleg and upleg
light-time solution, so a single Doppler observation involves four
light-time iterations.

Two separate approximations are at work when this is replaced by a
simple geometric derivative.

The first is the Taylor truncation. Expanding the range difference
around the centre of the interval gives

    [d0(t) - d0(t - tc)] / tc = d0_dot(t - tc/2) + O(tc^2)

so the average is centred on the mid-epoch, not on the time tag.
This is why the instantaneous value must be evaluated at
t - tc/2 rather than at t.

The second is the omission of the light-time correction itself,
whose leading term is

    -(rho_dot^2 + rho * rho_ddot) / c

The results below show that the second effect dominates the first
by more than an order of magnitude.

## Inputs

Identical to Exercise 6.1: circular polar orbit at 960 km altitude,
i = 97 deg, RAAN = 130.7 deg, ascending node passage at
1997/01/01 0h UTC, station at 11 deg East and 48 deg North on the
WGS84 ellipsoid.

Doppler count interval: tc = 1 s.
Evaluation times: every 3 minutes from 6 to 24 minutes past the
ascending node.

## Reusable Functions

New modules introduced by this exercise:

    src/satellite_orbits/observations/
        instantaneous_range_rate.py
        average_two_way_range_rate.py

Everything else is reused unchanged from Exercise 6.1, including
the light-time functions and the station and satellite state
providers.

A deliberate asymmetry: average_two_way_range_rate uses the
rigorous two_way_range model, while instantaneous_range_rate is
purely geometric and takes no light time into account. This is the
point of the exercise. Making both rigorous, or both geometric,
would remove the quantity under study.

The time tag convention follows the textbook: it refers to the end
of the count interval.

## Results

    UTC        Average  Instantaneous      Delta
    hh:mm:ss     [m/s]          [m/s]      [m/s]
    ----------------------------------------------
    00:06:00.0    -6385.695      -6385.555     -0.140
    00:09:00.0    -6113.313      -6113.162     -0.151
    00:12:00.0    -4575.874      -4575.721     -0.153
    00:15:00.0     1488.575       1488.742     -0.167
    00:18:00.0     5500.479       5500.640     -0.161
    00:21:00.0     6286.692       6286.841     -0.149
    00:24:00.0     6392.680       6392.814     -0.134

The range rate changes sign between 12 and 15 minutes, which marks
the culmination of the pass.

The difference between the two models stays between -0.134 and
-0.167 m/s and does not vanish anywhere, not even at culmination
where the range rate itself passes through zero. This confirms that
the difference is not driven by the Taylor truncation, which the
textbook bounds at 9 mm/s, but by the neglected light-time
correction.

## Engineering Conclusion

Real Doppler tracking systems reach a measurement accuracy of about
1 mm/s. The instantaneous model is in error by 0.13 to 0.17 m/s
here, two orders of magnitude worse. The simplified model is
therefore unusable for precise range rate modelling, even though it
is adequate for a quick visualisation of a pass.

## Verification

The model difference, which is the quantity the exercise asks for,
matches the textbook at every epoch with an absolute tolerance of
0.002 m/s.

The absolute range rates carry the same systematic reference-frame
residual as Exercise 6.1:

    Average:       [-0.0008, -0.0015, -0.0016, 0.0106,
                     0.0053, 0.0021, 0.0007] m/s
    Instantaneous: [-0.0010, -0.0014, -0.0010, 0.0100,
                     0.0052, 0.0020, 0.0010] m/s

This residual peaks at culmination, at about 11 mm/s, and falls to
roughly 1 mm/s at both ends of the pass. The shape is consistent
with the sidereal-time explanation documented in the Exercise 6.1
README: near culmination the range rate is most sensitive to the
geometry, so a small error in the Earth rotation angle has its
largest effect there.

Both the average and the instantaneous residuals follow the same
pattern and nearly cancel in their difference, which is why the
difference column reproduces the textbook exactly while the
absolute values do not. The absolute values are only checked
against a 0.02 m/s bound.

## Scope and Limitations

- The transponder turn-around ratio and the cycle counting
  electronics are not modelled. The measurement is treated at the
  level of equation 6.37, as the exercise intends.
- The rotational Doppler bias of a spinning spacecraft, equation
  6.40, is not included.
- No media corrections. Tropospheric refraction enters in
  Exercise 6.4.
- 