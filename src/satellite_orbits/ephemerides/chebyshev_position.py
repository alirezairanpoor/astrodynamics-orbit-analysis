import numpy as np


def chebyshev_position(
    mjd: float,
    mjd_start: float,
    mjd_end: float,
    coefficients: np.ndarray,
) -> np.ndarray:
    """
    Evaluate Cartesian position from Chebyshev coefficients.

    Inputs
    ------
    mjd : evaluation epoch
    mjd_start, mjd_end : coefficient validity interval
    coefficients : array with shape (N, 3)
        Each row contains [ax_k, ay_k, az_k] in km.

    All epochs must use the time scale of the coefficient set.

    Output
    ------
    Position [km], in the origin and frame of the coefficients.

    Convention: the coefficient of T0 is used in full.
    """

    coefficients = np.asarray(coefficients, dtype=float)

    if (
        coefficients.ndim != 2
        or coefficients.shape[1] != 3
        or coefficients.shape[0] == 0
    ):
        raise ValueError("coefficients must have shape (N, 3).")

    if mjd_end <= mjd_start:
        raise ValueError("mjd_end must be greater than mjd_start.")

    if not mjd_start <= mjd <= mjd_end:
        raise ValueError("Epoch is outside the coefficient interval.")

    # --------------------------------------------------------
    # SCALE TIME TO [-1, +1]
    # --------------------------------------------------------

    tau = (
        2.0 * (mjd - mjd_start)
        / (mjd_end - mjd_start)
        - 1.0
    )

    # --------------------------------------------------------
    # INITIAL TERMS
    # --------------------------------------------------------

    T_previous = 1.0
    position = coefficients[0].copy()

    if len(coefficients) == 1:
        return position

    T_current = tau
    position += coefficients[1] * T_current

    # --------------------------------------------------------
    # RECURRENCE AND POSITION SUM
    # --------------------------------------------------------

    for k in range(2, len(coefficients)):
        T_next = (
            2.0 * tau * T_current
            - T_previous
        )

        position += coefficients[k] * T_next

        T_previous = T_current
        T_current = T_next

    return position
