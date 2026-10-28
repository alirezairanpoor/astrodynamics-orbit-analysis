import numpy as np

from .integrate_rk4 import integrate_rk4


def integrate_gauss_jackson4(
    acceleration_function,
    r0,
    v0,
    duration,
    n_steps,
):
    """Explicit fourth-order Gauss-Jackson integration.

    Solves:
        r'' = acceleration_function(t, r)

    Acceleration must not depend on velocity.

    Starts at t = 0 and uses a fixed positive step.
    RK4 generates the three required backward starting states.

    Returns
    -------
    times : ndarray
    states : ndarray
        Position components followed by velocity components.
    acceleration_calls : int
        Includes initialization evaluations.
    """
    r0 = np.asarray(r0, dtype=float)
    v0 = np.asarray(v0, dtype=float)

    if (
        r0.ndim != 1
        or r0.size == 0
        or v0.shape != r0.shape
        or not np.all(np.isfinite(r0))
        or not np.all(np.isfinite(v0))
    ):
        raise ValueError("Provide matching finite position and velocity arrays.")

    if not np.isfinite(duration) or duration <= 0:
        raise ValueError("Duration must be finite and positive.")

    if not isinstance(n_steps, (int, np.integer)) or n_steps < 1:
        raise ValueError("n_steps must be a positive integer.")

    dimension = r0.size
    h = duration / n_steps
    acceleration_calls = 0

    def evaluate_acceleration(t, r):
        nonlocal acceleration_calls
        acceleration_calls += 1

        value = np.asarray(
            acceleration_function(t, r),
            dtype=float,
        )

        if value.shape != r0.shape or not np.all(np.isfinite(value)):
            raise ValueError("Acceleration has an invalid shape or value.")

        return value

    # Integrate backward using s = -t.
    # Therefore d(state)/ds = -d(state)/dt.
    def backward_derivative(s, state):
        r = state[:dimension]
        v = state[dimension:]

        return -np.concatenate((
            v,
            evaluate_acceleration(-s, r),
        ))

    _, backward_states = integrate_rk4(
        backward_derivative,
        np.concatenate((r0, v0)),
        3.0 * h,
        3,
    )

    # Accelerations at t = 0, -h, -2h, -3h.
    history = np.array([
        evaluate_acceleration(
            -j * h,
            backward_states[j, :dimension],
        )
        for j in range(4)
    ])

    # Backward differences at t = 0.
    differences = np.array([
        history[0],
        history[0] - history[1],
        history[0] - 2.0 * history[1] + history[2],
        history[0] - 3.0 * history[1]
        + 3.0 * history[2] - history[3],
    ])

    d0, d1, d2, d3 = differences

    # Initialize the first and second sums.
    first_sum = (
        v0 / h
        + d0 / 2.0
        + d1 / 12.0
        + d2 / 24.0
        + 19.0 * d3 / 720.0
    )

    second_sum = (
        r0 / h**2
        + first_sum
        - d0 / 12.0
        + d2 / 240.0
        + d3 / 240.0
    )

    times = np.linspace(0.0, duration, n_steps + 1)
    states = np.empty((n_steps + 1, 2 * dimension))
    states[0] = np.concatenate((r0, v0))

    for j in range(n_steps):
        d0, d1, d2, d3 = differences

        # Explicit position and velocity formulas.
        r_next = h**2 * (
            second_sum
            + d0 / 12.0
            + d1 / 12.0
            + 19.0 * d2 / 240.0
            + 3.0 * d3 / 40.0
        )

        v_next = h * (
            first_sum
            + d0 / 2.0
            + 5.0 * d1 / 12.0
            + 3.0 * d2 / 8.0
            + 251.0 * d3 / 720.0
        )

        states[j + 1] = np.concatenate((r_next, v_next))

        # One new acceleration evaluation per step.
        updated = np.empty_like(differences)
        updated[0] = evaluate_acceleration(
            times[j + 1],
            r_next,
        )

        for k in range(1, 4):
            updated[k] = updated[k - 1] - differences[k - 1]

        differences = updated

        first_sum = first_sum + differences[0]
        second_sum = second_sum + first_sum

    return times, states, acceleration_calls
