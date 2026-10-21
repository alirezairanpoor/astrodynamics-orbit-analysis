import numpy as np


def integrate_rk4(derivative_function, state0, duration, n_steps):
    """Integrate an ODE using fixed-step classical RK4.

    derivative_function is called as f(t, state).
    Integration starts at t = 0.

    Returns
    -------
    times : ndarray, shape (n_steps + 1,)
    states : ndarray, shape (n_steps + 1, len(state0))
    """
    state0 = np.asarray(state0, dtype=float)

    if state0.ndim != 1 or not np.all(np.isfinite(state0)):
        raise ValueError("Initial state must be a finite 1D array.")

    if not np.isfinite(duration) or duration <= 0:
        raise ValueError("Duration must be finite and positive.")

    if not isinstance(n_steps, (int, np.integer)) or n_steps < 1:
        raise ValueError("n_steps must be a positive integer.")

    times = np.linspace(0.0, duration, n_steps + 1)
    states = np.empty((n_steps + 1, state0.size))
    states[0] = state0

    def evaluate(t, state):
        result = np.asarray(
            derivative_function(t, state),
            dtype=float,
        )
        if result.shape != state0.shape:
            raise ValueError("Derivative has an incorrect shape.")
        if not np.all(np.isfinite(result)):
            raise ValueError("Derivative contains nonfinite values.")
        return result

    for j in range(n_steps):
        t = times[j]
        h = times[j + 1] - t
        state = states[j]

        k1 = evaluate(t, state)

        k2 = evaluate(
            t + h / 2.0,
            state + h * k1 / 2.0,
        )

        k3 = evaluate(
            t + h / 2.0,
            state + h * k2 / 2.0,
        )

        k4 = evaluate(
            t + h,
            state + h * k3,
        )

        states[j + 1] = state + h * (
            k1 + 2.0 * k2 + 2.0 * k3 + k4
        ) / 6.0

        if not np.all(np.isfinite(states[j + 1])):
            raise RuntimeError("Integration produced a nonfinite state.")

    return times, states
