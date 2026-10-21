# Exercise 2.3 — Osculating Elements

## Objective

Determine the six classical Keplerian orbital elements from a given ECI position and velocity state.

The given state vectors are

\[
\mathbf r_{ECI} =
\begin{bmatrix}
10000 \\
40000 \\
-5000
\end{bmatrix}
\text{ km}
\]

and

\[
\mathbf v_{ECI} =
\begin{bmatrix}
-1.5 \\
1.0 \\
-0.1
\end{bmatrix}
\text{ km/s}
\]

The required orbital elements are

\[
a,\ e,\ i,\ RAAN,\ AOP,\ M
\]

where:

- \(a\): semi-major axis
- \(e\): eccentricity
- \(i\): inclination
- \(RAAN\): right ascension of the ascending node
- \(AOP\): argument of perigee
- \(M\): mean anomaly

## Method

The orbital-plane orientation is first obtained from the angular momentum vector

\[
\mathbf h =
\mathbf r_{ECI}
\times
\mathbf v_{ECI}
\]

The direction of \(\mathbf h\) determines \(i\) and \(RAAN\).

The semi-latus rectum is calculated from

\[
p = \frac{h^2}{\mu}
\]

The semi-major axis is obtained from the vis-viva equation:

\[
a =
\left(
\frac{2}{r}
-
\frac{v^2}{\mu}
\right)^{-1}
\]

The eccentricity is then

\[
e =
\sqrt{
1-\frac{p}{a}
}
\]

The eccentric anomaly \(E\) is determined from the ECI position and velocity state.

The mean anomaly is obtained from Kepler's equation:

\[
M =
E
-
e \times \sin(E)
\]

The argument of latitude is calculated from the geometry of the position vector.

The true anomaly \(TA\) is calculated from \(E\) and \(e\).

Finally,

\[
AOP =
argument\_of\_latitude
-
TA
\]

## Results

The Python implementation gives approximately

```text
a    : 25015.181 km
e    : 0.7079772
i    : 6.971 deg
RAAN : 173.290 deg
AOP  : 91.553 deg
M    : 144.225 deg
```

These values agree with the reference solution given in Montenbruck and Gill.

## Conclusion

This exercise demonstrates the transformation

\[
\mathbf r_{ECI},\mathbf v_{ECI}
\longrightarrow
a,e,i,RAAN,AOP,M
\]

A single ECI position and velocity state contains enough information to determine the size, shape, orientation, and current phase of an unperturbed Keplerian orbit.

The reverse transformation is implemented separately through `elements_to_states()`, allowing consistent conversion in both directions:

\[
\mathbf r_{ECI},\mathbf v_{ECI}
\leftrightarrow
a,e,i,RAAN,AOP,M
\]

## Source

Oliver Montenbruck and Eberhard Gill, *Satellite Orbits: Models, Methods and Applications*, Exercise 2.3.
