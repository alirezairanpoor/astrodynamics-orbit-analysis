from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from satellite_orbits.astrodynamics.propagation.integrate_de import (
    integrate_de,
)


def derivative(t, state):
    """Normalized planar two-body equations: mu = 1."""
    r = state[:2]
    v = state[2:]

    acceleration = -r / np.linalg.norm(r)**3

    return np.concatenate((v, acceleration))


def reference_states(times, e):
    """Keplerian reference: mu = a = 1, pericenter at t = 0."""
    M = (np.asarray(times) + np.pi) % (2.0 * np.pi) - np.pi

    E = np.where(M >= 0.0, np.pi, -np.pi)

    for _ in range(60):
        correction = (
            (E - e * np.sin(E) - M)
            / (1.0 - e * np.cos(E))
        )
        E -= correction

        if np.max(np.abs(correction)) < 1e-14:
            break
    else:
        raise RuntimeError("Kepler iteration failed.")

    beta = np.sqrt(1.0 - e**2)
    denominator = 1.0 - e * np.cos(E)

    return np.column_stack((
        np.cos(E) - e,
        beta * np.sin(E),
        -np.sin(E) / denominator,
        beta * np.cos(E) / denominator,
    ))


def summarize(times, states, orders, calls, e, tolerance):
    reference = reference_states(times, e)

    errors = np.linalg.norm(
        states - reference,
        axis=1,
    )

    steps = np.diff(times)

    print(f"\nAbsolute local tolerance: {tolerance:.1e}")
    print(f"Accepted steps: {len(orders)}")
    print(f"Derivative calls: {calls}")
    print(f"Orders used: {orders.min()} to {orders.max()}")
    print(f"Smallest accepted step: {steps.min():.9e}")
    print(f"Largest accepted step:  {steps.max():.9e}")
    print(f"Final state error:       {errors[-1]:.9e}")
    print(f"Maximum sampled error:  {errors.max():.9e}")

    return errors


def main():
    e = 0.9
    duration = 20.0
    tolerance = 1e-8

    state0 = np.array([
        1.0 - e,
        0.0,
        0.0,
        np.sqrt((1.0 + e) / (1.0 - e)),
    ])

    print("Exercise 4.3 — Stepsize Control for Eccentric Orbits")
    print("Shampine-Gordon DE step algorithm")
    print("Normalized units: mu = a = 1")
    print("e = 0.9, initial time = 0, final time = 20")

    times, states, orders, calls = integrate_de(
        derivative,
        state0,
        duration,
        atol=tolerance,
    )

    errors = summarize(
        times,
        states,
        orders,
        calls,
        e,
        tolerance,
    )

    # Repeat with a tighter tolerance to check numerical sensitivity.
    times_fine, states_fine, orders_fine, calls_fine = integrate_de(
        derivative,
        state0,
        duration,
        atol=1e-10,
    )

    errors_fine = summarize(
        times_fine,
        states_fine,
        orders_fine,
        calls_fine,
        e,
        1e-10,
    )

    if errors_fine[-1] >= errors[-1]:
        raise RuntimeError(
            "Tighter tolerance did not improve final state error."
        )

    if errors_fine.max() >= errors.max():
        raise RuntimeError(
            "Tighter tolerance did not improve maximum sampled error."
        )

    # Actual accepted step lengths, not requested output intervals.
    steps = np.diff(times)

    # Associate each step with the distance at its starting point.
    radii = np.linalg.norm(
        states[:-1, :2],
        axis=1,
    )

    # Exclude startup during the first orbit and the final capped step.
    settled = times[:-1] >= 2.0 * np.pi
    settled[-1] = False

    near_pericenter = settled & (radii < 0.2)
    near_apocenter = settled & (radii > 1.8)

    if not (
        near_pericenter.any()
        and near_apocenter.any()
    ):
        raise RuntimeError(
            "Insufficient samples for step-size comparison."
        )

    h_pericenter = np.median(steps[near_pericenter])
    h_apocenter = np.median(steps[near_apocenter])

    if h_apocenter <= h_pericenter:
        raise RuntimeError("Unexpected step-size behavior.")

    print(
        "\nMedian step near pericenter: "
        f"{h_pericenter:.9e}"
    )
    print(
        "Median step near apocenter:  "
        f"{h_apocenter:.9e}"
    )
    print(
        "Median step ratio, apocenter/pericenter: "
        f"{h_apocenter / h_pericenter:.2f}"
    )

    print("\nTolerance-refinement checks: passed")
    print("Pericenter/apocenter step-size behavior: passed")
    print("Local tolerance is not a bound on global state error.")

    # Store the numerical data and plot beside this solution.
    output_dir = Path(__file__).resolve().parent / "outputs"
    output_dir.mkdir(exist_ok=True)

    csv_path = output_dir / "accepted_steps.csv"

    np.savetxt(
        csv_path,
        np.column_stack((
            times[:-1],
            times[1:],
            steps,
            radii,
            orders,
        )),
        delimiter=",",
        header="t_start,t_end,accepted_step,r_start,order",
        comments="",
    )

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(11, 4),
        constrained_layout=True,
    )

    axes[0].plot(
        times[:-1],
        steps,
        linewidth=0.9,
    )
    axes[0].set_xlabel("Time at step start (normalized)")
    axes[0].set_ylabel("Accepted step size (normalized)")
    axes[0].set_xlim(0.0, duration)

    axes[1].scatter(
        radii,
        steps,
        s=7,
        alpha=0.65,
    )
    axes[1].set_xlabel("Distance at step start (normalized)")
    axes[1].set_ylabel("Accepted step size (normalized)")
    axes[1].set_xlim(0.0, 2.0)

    for ax in axes:
        ax.set_ylim(bottom=0.0)
        ax.grid(alpha=0.25)

    fig.suptitle(
        "Exercise 4.3 — Shampine-Gordon DE, e = 0.9, atol = 1e-8"
    )

    plot_path = output_dir / "stepsize_control.png"
    fig.savefig(plot_path, dpi=180)

    print(f"\nCSV: {csv_path}")
    print(f"Plot: {plot_path}")

    plt.show()


if __name__ == "__main__":
    main()
