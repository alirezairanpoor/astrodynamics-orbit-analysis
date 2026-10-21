import numpy as np


def gravity_acceleration_cunningham(
    r_ecef: np.ndarray,
    C: np.ndarray,
    S: np.ndarray,
    n_max: int,
    m_max: int,
    mu: float,
    earth_radius: float,
) -> np.ndarray:
    """
    Compute gravitational acceleration using Cunningham recurrences.

    Inputs
    ------
    r_ecef : position in ECEF [km]
    C, S : unnormalized gravity coefficients, indexed by [n, m]
    n_max : maximum degree
    m_max : maximum order
    mu : gravitational parameter [km^3/s^2]
    earth_radius : reference radius of the gravity model [km]

    Returns
    -------
    Acceleration in ECEF [km/s^2].

    Central gravity is included when C[0, 0] = 1.
    """

    r_ecef = np.asarray(r_ecef, dtype=float)
    C = np.asarray(C, dtype=float)
    S = np.asarray(S, dtype=float)

    if r_ecef.shape != (3,):
        raise ValueError("r_ecef must have shape (3,).")

    if not 0 <= m_max <= n_max:
        raise ValueError("Require 0 <= m_max <= n_max.")

    for coefficients in (C, S):
        if (
            coefficients.ndim != 2
            or coefficients.shape[0] <= n_max
            or coefficients.shape[1] <= m_max
        ):
            raise ValueError("Coefficient arrays are too small.")

    # --------------------------------------------------------
    # POSITION-DEPENDENT FACTORS
    # --------------------------------------------------------

    r2 = float(r_ecef @ r_ecef)

    if r2 == 0.0:
        raise ValueError("Position must be nonzero.")

    x0, y0, z0 = earth_radius * r_ecef / r2
    rho = earth_radius**2 / r2

    # Acceleration requires V and W up to degree n_max + 1
    # and order m_max + 1.
    V = np.zeros((n_max + 2, m_max + 2))
    W = np.zeros_like(V)

    V[0, 0] = earth_radius / np.sqrt(r2)

    # --------------------------------------------------------
    # CUNNINGHAM RECURRENCES — EQUATIONS 3.29–3.31
    # --------------------------------------------------------

    for m in range(m_max + 2):

        # Diagonal terms: n = m
        if m > 0:
            V[m, m] = (2 * m - 1) * (
                x0 * V[m - 1, m - 1]
                - y0 * W[m - 1, m - 1]
            )

            W[m, m] = (2 * m - 1) * (
                x0 * W[m - 1, m - 1]
                + y0 * V[m - 1, m - 1]
            )

        # Remaining terms in column m: n > m
        for n in range(m + 1, n_max + 2):
            V[n, m] = (
                (2 * n - 1) * z0 * V[n - 1, m]
                / (n - m)
            )

            W[n, m] = (
                (2 * n - 1) * z0 * W[n - 1, m]
                / (n - m)
            )

            if n >= m + 2:
                factor = (
                    (n + m - 1) * rho
                    / (n - m)
                )

                V[n, m] -= factor * V[n - 2, m]
                W[n, m] -= factor * W[n - 2, m]

    # --------------------------------------------------------
    # ACCELERATION SUM — EQUATION 3.33
    # --------------------------------------------------------

    ax = 0.0
    ay = 0.0
    az = 0.0

    for n in range(n_max + 1):
        for m in range(min(n, m_max) + 1):

            c = C[n, m]
            s = S[n, m]

            if m == 0:
                ax -= c * V[n + 1, 1]
                ay -= c * W[n + 1, 1]

            else:
                factor = (n - m + 1) * (n - m + 2)

                ax += 0.5 * (
                    -c * V[n + 1, m + 1]
                    - s * W[n + 1, m + 1]
                    + factor * (
                        c * V[n + 1, m - 1]
                        + s * W[n + 1, m - 1]
                    )
                )

                ay += 0.5 * (
                    -c * W[n + 1, m + 1]
                    + s * V[n + 1, m + 1]
                    + factor * (
                        -c * W[n + 1, m - 1]
                        + s * V[n + 1, m - 1]
                    )
                )

            az -= (n - m + 1) * (
                c * V[n + 1, m]
                + s * W[n + 1, m]
            )

    a_ecef = (
        mu / earth_radius**2
        * np.array([ax, ay, az])
    )

    return a_ecef
