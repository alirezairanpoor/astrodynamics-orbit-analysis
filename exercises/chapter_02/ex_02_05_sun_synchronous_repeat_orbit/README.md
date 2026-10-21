# Exercise 2.5 — Sun-Synchronous Repeat Orbit

## Objective

Determine the altitude and inclination of a circular Earth orbit that satisfies two requirements simultaneously:

1. The orbit must be **sun-synchronous**.
2. The ground track must repeat after

```text
K = 3 days
N = 43 draconic orbits
```

The required outputs are therefore:

```text
a
h
i
```

where:

- `a` = semi-major axis
- `h` = altitude
- `i` = inclination

The exercise also asks us to verify that the longitude of the ascending node repeats after the specified number of crossings.

---

## 1. Repeat-Ground-Track Requirement

The orbit must complete 43 draconic orbits in 3 days.

Therefore the required draconic period is

\[
T_N
=
\frac{K}{N}
\text{ day}
\]

or

\[
T_N
=
\frac{3}{43}
\text{ day}
\]

which gives approximately

\[
T_N
=
6027.907
\text{ s}
\]

This value is the **target draconic period** and remains fixed throughout the iteration.

A draconic period is the time between two consecutive crossings of the same orbital node.

---

## 2. Initial Estimate of `a`

At first, the difference between the draconic period and the Keplerian orbital period is neglected:

\[
T_N \approx T_0
\]

with

\[
T_0
=
2\pi
\sqrt{
\frac{a^3}{\mu}
}
\]

Solving for `a` gives an initial estimate:

\[
a
=
\left[
\mu
\left(
\frac{T_N}{2\pi}
\right)^2
\right]^{1/3}
\]

The initial result is approximately

```text
a = 7158.747 km
```

The book explicitly recommends using this approximation as the first guess before applying the secular perturbation corrections.

---

## 3. Sun-Synchronous Inclination

For a sun-synchronous orbit, the secular RAAN rate caused by the Earth's \(J_2\) perturbation must match the apparent annual motion of the Sun:

\[
\dot{RAAN}
=
0.985647240^\circ/day
\]

The nodal precession caused by \(J_2\) is

\[
\dot{RAAN}
=
-\frac{3}{2}
n_0
J_2
\left(
\frac{R_E}{a}
\right)^2
\cos i
\]

where

\[
n_0
=
\sqrt{
\frac{\mu}{a^3}
}
\]

For each estimate of `a`, this equation is solved for `i`.

---

## 4. J2 Corrections

The draconic period is not exactly equal to the Keplerian period.

The exercise includes two additional secular effects:

### Mean-motion correction

\[
\Delta n
=
-\frac{3}{4}
n_0
J_2
\left(
\frac{R_E}{a}
\right)^2
\left(
1-3\cos^2 i
\right)
\]

The effective mean motion is therefore

\[
n
=
n_0+\Delta n
\]

### Argument-of-perigee rate

\[
\dot{AOP}
=
-\frac{3}{4}
n_0
J_2
\left(
\frac{R_E}{a}
\right)^2
\left(
1-5\cos^2 i
\right)
\]

These perturbations modify the node-to-node orbital period.

---

## 5. Draconic Period

The corrected draconic period is

\[
T_N
=
\frac{2\pi}
{n+\dot{AOP}}
\]

and since

\[
n=n_0+\Delta n
\]

we obtain

\[
T_N
=
\frac{2\pi}
{n_0+\Delta n+\dot{AOP}}
\]

The required \(T_N\) is already known from the repeat-ground-track condition.

Therefore we solve for a corrected Keplerian mean motion:

\[
n_0
=
\frac{2\pi}{T_N}
-
\Delta n
-
\dot{AOP}
\]

and then obtain a new semi-major axis from

\[
a
=
\left(
\frac{\mu}{n_0^2}
\right)^{1/3}
\]

This new value of `a` produces a new value of `i`, which produces new \(J_2\) corrections.

The process is repeated until the result converges.

---

## 6. Iteration Process

The complete iteration is:

```text
K / N
  ↓
TN
  ↓
assume TN ≈ T0
  ↓
initial a
  ↓
SSO condition
  ↓
i
  ↓
J2 corrections
  ↓
delta_n + AOP_rate
  ↓
corrected n0
  ↓
new a
  ↓
new i
  ↓
repeat
```

The computed values are approximately:

```text
Iteration 0
a = 7158.747 km
h = 780.610 km
i = 98.518 deg

Iteration 1
a = 7153.139 km
h = 775.002 km
i = 98.495 deg

Iteration 2
a = 7153.120 km
h = 774.983 km
i = 98.495 deg
```

The reference solution in the book gives approximately:

```text
a = 7153.123 km
h = 774.986 km
i = 98.498 deg
```

which closely agrees with the implementation.

---

## 7. Keplerian Period vs Draconic Period

The two periods are different.

The Keplerian period is

\[
T_0
=
\frac{2\pi}{n_0}
\]

while the draconic period is

\[
T_N
=
\frac{2\pi}
{n_0+\Delta n+\dot{AOP}}
\]

Therefore

\[
T_N \neq T_0
\]

In the secular model used in this exercise, both are constant for the final orbit
