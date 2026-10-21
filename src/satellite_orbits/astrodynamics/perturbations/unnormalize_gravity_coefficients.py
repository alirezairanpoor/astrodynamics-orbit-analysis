from math import factorial, sqrt

import numpy as np


def unnormalize_gravity_coefficients(
    C_normalized: np.ndarray,
    S_normalized: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Convert fully normalized gravity coefficients to the
    unnormalized convention used by the Cunningham function.

    Uses the normalization convention of Montenbruck & Gill.

    Inputs:
        C_normalized, S_normalized:
            Matching two-dimensional coefficient arrays.
            Coefficients must already be dimensionless;
            apply any table scale factor before calling.

    Returns:
        C, S:
            New arrays containing unnormalized coefficients.

    Intended for the degree-20 model used in this project.
    """

    C_normalized = np.asarray(C_normalized, dtype=float)
    S_normalized = np.asarray(S_normalized, dtype=float)

    if (
        C_normalized.ndim != 2
        or S_normalized.ndim != 2
        or C_normalized.shape != S_normalized.shape
    ):
        raise ValueError(
            "C and S must be matching two-dimensional arrays."
        )

    C = np.zeros_like(C_normalized)
    S = np.zeros_like(S_normalized)

    n_max = C.shape[0] - 1
    m_max = C.shape[1] - 1

    for n in range(n_max + 1):
        for m in range(min(n, m_max) + 1):

            order_factor = 1.0 if m == 0 else 2.0

            normalization_factor = sqrt(
                order_factor
                * (2 * n + 1)
                * factorial(n - m)
                / factorial(n + m)
            )

            C[n, m] = (
                normalization_factor
                * C_normalized[n, m]
            )

            S[n, m] = (
                normalization_factor
                * S_normalized[n, m]
            )

    return C, S
