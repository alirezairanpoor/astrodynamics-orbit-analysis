# Exercise 2.2 — Solving Kepler's Equation

## Objective

Solve Kepler's equation

\[
M = E - e\sin E
\]

using two iterative numerical methods:

- Newton's method
- Fixed-point iteration

The exercise compares the convergence behavior and computational effort of both methods for an eccentricity of

\[
e = 0.72
\]

at two different mean anomalies:

\[
M = 4^\circ
\]

and

\[
M = 50^\circ
\]

## Theory

Kepler's equation cannot generally be solved explicitly for the eccentric anomaly \(E\), so an iterative numerical method is required.

### Newton's Method

Kepler's equation is written as

\[
f(E)=E-e\sin E-M=0
\]

and Newton's method updates the estimate according to

\[
E_{i+1}
=
E_i-
\frac{E_i-e\sin E_i-M}
{1-e\cos E_i}
\]

Each Newton iteration evaluates both a sine and a cosine function.

### Fixed-Point Iteration

Rearranging Kepler's equation gives

\[
E=M+e\sin E
\]

which leads to the iteration

\[
E_{i+1}=M+e\sin E_i
\]

Each fixed-point iteration requires only one sine evaluation.

## Results

### Case 1 — \(M=4^\circ\)

```text
Newton method
E                 : 13.9335999836 deg
Iterations        : 4
Trig evaluations  : 8

Fixed-point method
E                 : 13.9335999737 deg
Iterations        : 58
Trig evaluations  : 58
```

Both methods converge to approximately

\[
E \approx 13.9336^\circ
\]

Newton's method converges much faster in this case.

### Case 2 — \(M=50^\circ\)

```text
Newton method
E                 : 91.2432498975 deg
Iterations        : 6
Trig evaluations  : 12

Fixed-point method
E                 : 91.2432498975 deg
Iterations        : 8
Trig evaluations  : 8
```

Both methods converge to approximately

\[
E \approx 91.24325^\circ
\]

In this case, fixed-point iteration requires slightly more iterations but fewer total trigonometric evaluations.

## Conclusion

The exercise demonstrates that the most efficient numerical method depends on both the convergence rate and the computational cost of each iteration.

Newton's method generally converges faster because it uses information about the slope of the equation. However, each Newton iteration is computationally more expensive because both sine and cosine must be evaluated.

Fixed-point iteration is simpler and cheaper per iteration. Near the pericenter of a relatively eccentric orbit, it may converge slowly, while farther from the pericenter it can become competitive with Newton's method.

The exercise therefore illustrates an important numerical-analysis principle:

> A method requiring fewer iterations is not necessarily computationally cheaper.

The exact iteration counts depend on the selected convergence tolerance and stopping criterion.

## Source

Oliver Montenbruck and Eberhard Gill, *Satellite Orbits: Models, Methods and Applications*, Exercise 2.2.