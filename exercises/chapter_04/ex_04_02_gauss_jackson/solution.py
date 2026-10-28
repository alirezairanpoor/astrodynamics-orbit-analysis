import numpy as np

from satellite_orbits.astrodynamics.propagation.integrate_gauss_jackson4 import (
    integrate_gauss_jackson4,
)


def reference_state(t, e):
    """Keplerian reference for the normalized planar problem."""
    M = t
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


def acceleration_function(t, r):
    """Normalized central gravity: mu = 1."""
    return -r / np.linalg.norm(r)**3


def main():
    e = 0.1
    duration = 20.0

    r0 = np.array([1.0 - e, 0.0])
    v0 = np.array([0.0, np.sqrt((1.0 + e) / (1.0 - e))])

    state_reference = reference_state(duration, e)

    step_counts = [
        100, 300, 600, 1000,
        1500, 2000, 3000, 4000,
    ]

    book_errors = np.array([
        2.456e-3,
        3.363e-6,
        7.705e-8,
        4.539e-9,
        5.174e-10,
        1.166e-10,
        1.499e-11,
        3.577e-12,
    ])

    measured_errors = []

    print("Exercise 4.2 — Fourth-order Gauss-Jackson")
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
        _, states, calls = integrate_gauss_jackson4(
            acceleration_function,
            r0,
            v0,
            duration,
            n_steps,
        )

        error = float(
            np.linalg.norm(states[-1] - state_reference)
        )
        measured_errors.append(error)

        if calls != n_steps + 16:
            raise RuntimeError("Unexpected acceleration evaluation count.")

        print(
            f"{n_steps:6d} "
            f"{duration / n_steps:10.6f} "
            f"{calls:7d} "
            f"{error:14.6e} "
            f"{-np.log10(error):8.2f}"
        )

    # Printed reference rounding and floating-point differences
    # matter particularly for the smallest errors.
    np.testing.assert_allclose(
        measured_errors,
        book_errors,
        rtol=2e-3,
        atol=5e-13,
    )

    print("\nAcceleration evaluation count: passed")
    print("Book comparison within stated tolerances: passed")


if __name__ == "__main__":
    main()
