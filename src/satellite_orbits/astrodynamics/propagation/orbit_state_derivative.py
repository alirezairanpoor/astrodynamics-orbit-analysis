import numpy as np


def orbit_state_derivative(t, state, acceleration_function):
    """Return the time derivative of the orbital state.

    Parameters
    ----------
    t : float
        Elapsed time since the initial epoch [s].
    state : array_like, shape (6,)
        Position [km] followed by velocity [km/s].
    acceleration_function : callable
        Called as acceleration_function(t, r_eci, v_eci).
        Must return acceleration in ECI [km/s^2].

    Returns
    -------
    ndarray, shape (6,)
        Velocity [km/s] followed by acceleration [km/s^2].
    """
    state = np.asarray(state, dtype=float)

    if state.shape != (6,):
        raise ValueError("State must have shape (6,).")

    if not np.all(np.isfinite(state)):
        raise ValueError("State must contain finite values.")

    r_eci = state[:3]
    v_eci = state[3:]

    a_eci = np.asarray(
        acceleration_function(t, r_eci, v_eci),
        dtype=float,
    )

    if a_eci.shape != (3,) or not np.all(np.isfinite(a_eci)):
        raise ValueError(
            "Acceleration must be a finite vector with shape (3,)."
        )

    return np.concatenate((v_eci, a_eci))
