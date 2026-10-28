"""Shampine-Gordon variable-order Adams PECE step algorithm.

Python adaptation of STEP from Shampine and Gordon's ODE package.

Source:
https://people.math.sc.edu/Burkardt/f77_src/ode/ode.f

SPDX-License-Identifier: LGPL-3.0-or-later

Scope:
- Forward integration.
- Scalar absolute tolerance.
- Accepted-step output without dense interpolation.
- Final step limited to the requested endpoint.

The complete legacy DE/ODE driver interface is not reproduced.
"""

import numpy as np


def integrate_de(
    function,
    state0,
    duration,
    atol=1e-8,
    max_steps=100000,
):
    """Integrate a first-order system using the DE step algorithm.

    Parameters
    ----------
    function : callable
        Called as function(t, state).
        Returns the state derivative.
    state0 : array_like
        Initial state at t = 0.
    duration : float
        Positive propagation duration.
    atol : float
        Absolute local error-control tolerance.
    max_steps : int
        Maximum number of accepted steps.

    Returns
    -------
    times : ndarray
        Initial time followed by accepted step endpoints.
    states : ndarray
        State at each returned time.
    orders : ndarray
        orders[j] is the order used between times[j]
        and times[j + 1].
    calls : int
        Total derivative evaluations, including rejected attempts.

    Notes
    -----
    Error control uses an unweighted Euclidean norm.
    No relative tolerance is used.

    State components must be scaled appropriately.
    This interface is intended for normalized problems.
    """
    y = np.asarray(state0, dtype=float).copy()

    if (
        y.ndim != 1
        or not y.size
        or not np.all(np.isfinite(y))
    ):
        raise ValueError(
            "state0 must be a finite nonempty 1D array."
        )

    if not np.isfinite(duration) or duration <= 0:
        raise ValueError(
            "duration must be finite and positive."
        )

    if not np.isfinite(atol) or atol <= 0:
        raise ValueError(
            "atol must be finite and positive."
        )

    if (
        not isinstance(max_steps, (int, np.integer))
        or max_steps < 1
    ):
        raise ValueError(
            "max_steps must be a positive integer."
        )

    calls = 0

    def evaluate(t, state):
        nonlocal calls
        calls += 1

        result = np.asarray(
            function(t, state),
            dtype=float,
        )

        if (
            result.shape != y.shape
            or not np.all(np.isfinite(result))
        ):
            raise ValueError(
                "Invalid derivative shape or value."
            )

        return result

    # Preserve one-based coefficient indexing from the algorithm.
    phi = np.zeros((y.size, 17))
    alpha = np.zeros(14)
    beta = np.ones(14)
    psi = np.zeros(14)
    sigma = np.zeros(14)
    v = np.zeros(14)
    w = np.zeros(14)
    g = np.zeros(14)

    g[1:3] = [1.0, 0.5]
    sigma[1] = 1.0

    gstr = np.array([
        0.0,
        0.500,
        0.0833,
        0.0417,
        0.0264,
        0.0188,
        0.0143,
        0.0114,
        0.00936,
        0.00789,
        0.00679,
        0.00592,
        0.00524,
        0.00468,
    ])

    eps = np.finfo(float).eps
    half_tol = 0.5 * atol
    t = 0.0

    derivative = evaluate(t, y)
    phi[:, 1] = derivative

    derivative_norm = np.linalg.norm(derivative)
    h = duration

    if (
        derivative_norm > 0
        and atol < 16 * derivative_norm * h**2
    ):
        h = 0.25 * np.sqrt(atol / derivative_norm)

    roundoff = 2 * eps * np.linalg.norm(y)
    compensated = half_tol <= 100 * roundoff

    phase1 = True
    k = 1
    kold = 0
    ns = 0
    hold = 0.0

    times = [t]
    states = [y.copy()]
    orders = []

    for _ in range(max_steps):
        h = min(h, duration - t)

        if half_tol < 2 * eps * np.linalg.norm(y):
            raise RuntimeError(
                "Tolerance is below the roundoff limit."
            )

        failures = 0

        while True:
            if (
                h <= 0
                or t + h == t
                or h < 4 * eps * abs(t)
            ):
                raise RuntimeError(
                    "Step size is below machine resolution."
                )

            kp1 = k + 1
            kp2 = k + 2

            if h != hold:
                ns = 0

            if ns <= kold:
                ns += 1

            # Coefficients for the nonuniform step history.
            if k >= ns:
                beta[ns] = 1.0
                alpha[ns] = 1.0 / ns
                temp = h * ns
                sigma[ns + 1] = 1.0

                for i in range(ns + 1, k + 1):
                    old_psi = psi[i - 1]
                    psi[i - 1] = temp

                    beta[i] = (
                        beta[i - 1] * temp / old_psi
                    )

                    temp = old_psi + h
                    alpha[i] = h / temp

                    sigma[i + 1] = (
                        i * alpha[i] * sigma[i]
                    )

                psi[k] = temp

                if ns == 1:
                    for i in range(1, k + 1):
                        v[i] = 1.0 / (i * (i + 1))
                        w[i] = v[i]

                else:
                    if k > kold:
                        v[k] = 1.0 / (k * kp1)

                        for j in range(1, ns - 1):
                            i = k - j
                            v[i] -= alpha[j + 1] * v[i + 1]

                    for i in range(1, kp1 - ns + 1):
                        v[i] -= alpha[ns] * v[i + 1]
                        w[i] = v[i]

                    g[ns + 1] = w[1]

                for i in range(ns + 2, kp1 + 1):
                    for j in range(1, kp2 - i + 1):
                        w[j] -= alpha[i - 1] * w[j + 1]

                    g[i] = w[1]

            # Predict using modified divided differences.
            for i in range(ns + 1, k + 1):
                phi[:, i] *= beta[i]

            phi[:, kp2] = phi[:, kp1]
            phi[:, kp1] = 0.0

            increment = np.zeros_like(y)

            for i in range(k, 0, -1):
                increment += g[i] * phi[:, i]
                phi[:, i] += phi[:, i + 1]

            if compensated:
                tau = h * increment - phi[:, 15]
                predicted = y + tau
                phi[:, 16] = (predicted - y) - tau
            else:
                predicted = y + h * increment

            yp = evaluate(t + h, predicted)
            delta = yp - phi[:, 1]

            scaled = h * np.linalg.norm(delta)
            error = scaled * (g[k] - g[kp1])
            erk = scaled * sigma[kp1] * gstr[k]

            erkm1 = 0.0
            erkm2 = 0.0

            if k >= 2:
                erkm1 = (
                    h
                    * sigma[k]
                    * gstr[k - 1]
                    * np.linalg.norm(phi[:, k] + delta)
                )

            if k >= 3:
                erkm2 = (
                    h
                    * sigma[k - 1]
                    * gstr[k - 2]
                    * np.linalg.norm(phi[:, k - 1] + delta)
                )

            knew = k

            if k >= 3 and max(erkm1, erkm2) <= erk:
                knew = k - 1
            elif k == 2 and erkm1 <= 0.5 * erk:
                knew = k - 1

            if error <= atol:
                break

            # Reject the step and restore the history.
            phase1 = False

            for i in range(1, k + 1):
                phi[:, i] = (
                    phi[:, i] - phi[:, i + 1]
                ) / beta[i]

            for i in range(2, k + 1):
                psi[i - 1] = psi[i] - h

            failures += 1
            factor = 0.5

            if failures >= 3:
                knew = 1

            if failures > 3 and half_tol < 0.25 * erk:
                factor = np.sqrt(half_tol / erk)

            h *= factor
            k = knew

        # Accept and correct the predicted state.
        kold = k
        hold = h

        correction = h * g[kp1] * delta

        if compensated:
            rho = correction - phi[:, 16]
            y = predicted + rho
            phi[:, 15] = (y - predicted) - rho
        else:
            y = predicted + correction

        t += h
        yp = evaluate(t, y)

        # Update differences for the next step.
        phi[:, kp1] = yp - phi[:, 1]
        phi[:, kp2] = phi[:, kp1] - phi[:, kp2]

        for i in range(1, k + 1):
            phi[:, i] += phi[:, kp1]

        times.append(t)
        states.append(y.copy())
        orders.append(kold)

        if t >= duration:
            return (
                np.array(times),
                np.array(states),
                np.array(orders),
                calls,
            )

        # Select the next order.
        if knew == k - 1 or k == 12:
            phase1 = False

        if phase1:
            k += 1
            erk = 0.0

        elif knew == k - 1:
            k -= 1
            erk = erkm1

        elif kp1 <= ns:
            erkp1 = (
                h
                * gstr[kp1]
                * np.linalg.norm(phi[:, kp2])
            )

            if k == 1:
                if erkp1 < 0.5 * erk:
                    k += 1
                    erk = erkp1

            elif erkm1 <= min(erk, erkp1):
                k -= 1
                erk = erkm1

            elif erkp1 < erk and k < 12:
                k += 1
                erk = erkp1

        # Select the next step size.
        hnew = 2.0 * h

        if (
            not phase1
            and half_tol < erk * 2.0**(k + 1)
        ):
            hnew = h

            if half_tol < erk:
                factor = (
                    half_tol / erk
                ) ** (1.0 / (k + 1))

                hnew = h * max(
                    0.5,
                    min(0.9, factor),
                )

        h = hnew

    raise RuntimeError(
        "Maximum accepted-step count exceeded."
    )
