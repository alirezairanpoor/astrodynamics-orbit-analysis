# Exercise 2.4 — Topocentric Satellite Motion

## Objective

Compute the apparent motion of a satellite relative to a ground station.

The required observables are:

- Azimuth
- Elevation
- Range
- Visibility duration
- Maximum elevation

The satellite is propagated for 30 minutes from the reference epoch.

## Given Data

### Satellite

```text
Altitude  = 960 km
e         = 0
i         = 97 deg
RAAN      = 130.7 deg
```

At the reference epoch, the satellite crosses the equator at the ascending node.

Because the orbit is circular, `AOP` is undefined. For the numerical representation used here,

```text
AOP = 0
M0  = 0
```

is chosen so that the satellite is located at the ascending node at the reference epoch.

### Ground Station

```text
Longitude = +11 deg East
Latitude  = +48 deg
```

### Reference Epoch

```text
MJD = 50449.0
```

corresponding to 1 January 1997.

The exercise data follow Montenbruck and Gill, Exercise 2.4.

## Method

The computation is performed through a sequence of independent transformations.

### 1. Propagate the Mean Anomaly

For an unperturbed Keplerian orbit,

\[
M(t)=M_0+n\Delta t
\]

with

\[
n=\sqrt{\frac{\mu}{a^3}}
\]

This is implemented by:

```python
propagate_mean_anomaly()
```

### 2. Orbital Elements to ECI

The orbital elements

```text
a / e / i / RAAN / AOP / M
```

are converted to the satellite position and velocity in the Earth-Centered Inertial frame:

```text
r_eci
v_eci
```

using:

```python
elements_to_states()
```

### 3. Compute Greenwich Hour Angle

The current Modified Julian Date is

\[
MJD =
MJD_0
+
\frac{\Delta t}{86400}
\]

and is converted to the Greenwich Hour Angle.

### 4. ECI to ECEF

The satellite position is rotated from the inertial frame into the Earth-fixed frame:

```text
r_eci
  ↓
r_ecef
```

using:

```python
eci_to_ecef()
```

### 5. ECEF to Local ENZ

The satellite position relative to the ground station is expressed in the local East-North-Zenith frame:

```text
r_ecef
   ↓
s_enz
```

using:

```python
ecef_to_enz()
```

### 6. Local Observables

The local ENZ vector is converted to:

```text
Azimuth
Elevation
Range
```

using:

```python
enz_to_az_el_range()
```

The complete transformation chain is therefore

```text
a / e / i / RAAN / AOP / M0
              ↓
     propagate_mean_anomaly()
              ↓
              M
              ↓
      elements_to_states()
              ↓
            r_eci
              ↓
         eci_to_ecef()
              ↓
            r_ecef
              ↓
         ecef_to_enz()
              ↓
            s_enz
              ↓
      enz_to_az_el_range()
              ↓
 Azimuth / Elevation / Range
```

## Visibility

The satellite is considered visible when

\[
Elevation > 0
\]

The rise and set times are estimated by linear interpolation between consecutive simulation samples surrounding the horizon crossing.

## Results

The simulation gives a maximum elevation of approximately

```text
54.6 deg
```

The satellite rises shortly after minute 6 and sets shortly before minute 23.

The reference solution reports a visibility of approximately 18 minutes and a peak elevation of roughly 55 degrees.

The computed trajectory closely reproduces the tabulated reference values. For example, the book gives an elevation of 54.9 degrees at minute 14 and 53.7 degrees at minute 15.

## Conclusion

This exercise connects orbital motion with observations made from
