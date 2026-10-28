import numpy as np

from satellite_orbits.astrodynamics.propagation.integrate_rk4 import (
    integrate_rk4,
)


def reference_state(t, e):
    """Keplerian state for the normalized planar test problem.

    mu = 1, a = 1, initial mean anomaly = 0.
    All quantities use normalized units.
    """
    M = t

    # Solve Kepler's equation: E - e*sin(E) = M.
    E = M

    for _ in range(30):
        correction = (
            (E - e * np.sin(E) - M)
            / (1.0 - e * np.cos(E))
        )
        E -= correction

        if abs(correction) < 1e-14:
            break
    else:
        raise RuntimeError("Kepler iteration did not converge.")

    denominator = 1.0 - e * np.cos(E)
    beta = np.sqrt(1.0 - e**2)

    return np.array([
        np.cos(E) - e,
        beta * np.sin(E),
        -np.sin(E) / denominator,
        beta * np.cos(E) / denominator,
    ])


def main():
    e = 0.1
    duration = 20.0

    state0 = np.array([
        1.0 - e,
        0.0,
        0.0,
        np.sqrt((1.0 + e) / (1.0 - e)),
    ])

    state_reference = reference_state(duration, e)

    step_counts = [
        50, 100, 250, 500,
        750, 1000, 1500, 2000,
    ]

    book_errors = np.array([
        1.953e-1,
        6.663e-3,
        9.051e-5,
        4.012e-6,
        6.847e-7,
        1.996e-7,
        3.608e-8,
        1.089e-8,
    ])

    measured_errors = []

    print("Exercise 4.1 — Fourth-order Runge-Kutta")
    print("Normalized planar two-body problem")
    print("e = 0.1, initial time = 0, final time = 20")
    print()
    print(
        f"{'Steps':>6} "
        f"{'Step size':>10} "
        f"{'Calls':>7} "
        f"{'State error':>14} "
        f"{'Digits':>8}"
    )
    print("-" * 49)

    for n_steps in step_counts:
        function_calls = 0

        def derivative_function(t, state):
            nonlocal function_calls
            function_calls += 1

            r = state[:2]
            v = state[2:]

            acceleration = -r / np.linalg.norm(r)**3

            return np.concatenate((v, acceleration))

        _, states = integrate_rk4(
            derivative_function,
            state0,
            duration,
            n_steps,
        )

        error = float(
            np.linalg.norm(states[-1] - state_reference)
        )
        digits = -np.log10(error)

        measured_errors.append(error)

        if function_calls != 4 * n_steps:
            raise RuntimeError("Unexpected RK4 evaluation count.")

        print(
            f"{n_steps:6d} "
            f"{duration / n_steps:10.6f} "
            f"{function_calls:7d} "
            f"{error:14.6e} "
            f"{digits:8.2f}"
        )

    # Allow for rounding of the textbook's printed values.
    np.testing.assert_allclose(
        measured_errors,
        book_errors,
        rtol=1e-3,
        atol=0.0,
    )

    print("\nFunction evaluation count: passed")
    print("Agreement with rounded textbook errors: passed")


if __name__ == "__main__":
    main()