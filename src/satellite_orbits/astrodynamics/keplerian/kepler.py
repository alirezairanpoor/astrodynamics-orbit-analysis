from math import sin, cos


def solve_kepler_newton(
    M: float,
    e: float,
    tolerance: float = 1e-10,
) -> tuple[float, int, int]:
    """
    Solve Kepler's equation using Newton's method.

    Parameters
    ----------
    M : float
        Mean anomaly [rad].

    e : float
        Orbital eccentricity [-].

    tolerance : float
        Convergence tolerance [rad].

    Returns
    -------
    tuple
        E : Eccentric anomaly [rad].
        iterations : Number of iterations.
        trig_evaluations : Number of trigonometric evaluations.
    """

    E = M
    iterations = 0
    trig_evaluations = 0

    while True:
        correction = (
            E - e * sin(E) - M
        ) / (
            1 - e * cos(E)
        )

        trig_evaluations += 2

        E_new = E - correction
        iterations += 1

        if abs(E_new - E) < tolerance:
            return E_new, iterations, trig_evaluations

        E = E_new


def solve_kepler_fixed_point(
    M: float,
    e: float,
    tolerance: float = 1e-10,
) -> tuple[float, int, int]:
    """
    Solve Kepler's equation using fixed-point iteration.
    """

    E = M
    iterations = 0
    trig_evaluations = 0

    while True:
        E_new = M + e * sin(E)

        trig_evaluations += 1
        iterations += 1

        if abs(E_new - E) < tolerance:
            return E_new, iterations, trig_evaluations

        E = E_new
