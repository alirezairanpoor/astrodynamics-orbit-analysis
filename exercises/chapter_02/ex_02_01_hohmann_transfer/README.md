# Exercise 2.1 — Hohmann Orbit Raising

## Objective

Compute the total velocity increment required for a Hohmann transfer between two circular Earth orbits.

The example considered here raises a spacecraft from an altitude of:

- Initial altitude: 750 km
- Final altitude: 775 km

The implementation also verifies the first-order approximation for a small change in orbital radius.

## Method

The orbital radii are obtained from

\[
r_1 = R_E + h_1
\]

\[
r_2 = R_E + h_2
\]

The circular velocities are

\[
v_1 = \sqrt{\frac{\mu}{r_1}}
\]

\[
v_2 = \sqrt{\frac{\mu}{r_2}}
\]

The semi-major axis of the Hohmann transfer ellipse is

\[
a_t = \frac{r_1+r_2}{2}
\]

Using the vis-viva equation, the transfer-orbit velocities at perigee and apogee are calculated.

The two required burns are then

\[
\Delta v_1 = v_p-v_1
\]

\[
\Delta v_2 = v_2-v_a
\]

and the total velocity increment is

\[
\Delta v_{\text{total}}=\Delta v_1+\Delta v_2
\]

For a small change in orbital radius, the exercise also considers the first-order approximation

\[
\Delta v_{\text{total}} \approx |v_1-v_2|
\]

## Result

The Python implementation gives:

```text
Initial orbital radius : 7128.137 km
Final orbital radius   : 7153.137 km

Initial circular velocity : 7.477921 km/s
Final circular velocity   : 7.464842 km/s

First burn  Δv1 : 6.542 m/s
Second burn Δv2 : 6.537 m/s
Total Δv         : 13.079 m/s

First-order approximation |v1 - v2| : 13.079 m/s
```

The result agrees with the reference solution in Montenbruck & Gill.

## Source

Oliver Montenbruck and Eberhard Gill, *Satellite Orbits: Models, Methods and Applications*, Exercise 2.1.