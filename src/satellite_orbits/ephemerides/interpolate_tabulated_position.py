import numpy as np


def interpolate_tabulated_position(
    times: np.ndarray,
    positions: np.ndarray,
    time: float,
    order: int = 9,
) -> np.ndarray | float:
    """
    Interpolate a tabulated ephemeris with a Lagrange polynomial.

    The polynomial is evaluated by Neville's algorithm, which
    builds it by repeated weighted averaging instead of forming
    the Lagrange basis explicitly.

    The nodes used are selected as a window centred on the
    requested time. Centring matters: evaluating near the edge of
    the node set raises the interpolation error sharply.

    Precise GPS ephemerides are tabulated at 15 min intervals, for
    which an order near 9 keeps the interpolation error well below
    the centimetre level. Linear interpolation over the same
    interval would be in error by tens of kilometres.

    Parameters
    ----------
    times : np.ndarray
        Tabulated epochs, shape (n,), strictly increasing. Any
        consistent time argument may be used.

    positions : np.ndarray
        Tabulated values, shape (n,) for scalar data or (n, m) for
        vector data such as Cartesian positions.

    time : float
        Epoch at which the interpolated value is requested.

    order : int
        Polynomial order. The window holds order + 1 nodes.

    Returns
    -------
    np.ndarray or float
        Interpolated value, matching the shape of a single row of
        positions.
    """

    times = np.asarray(times, dtype=float)
    positions = np.asarray(positions, dtype=float)

    number_of_nodes = order + 1

    if number_of_nodes > len(times):
        raise ValueError(
            "The requested order needs more nodes than available."
        )

    # Centre the node window on the requested epoch, then clamp it
    # to the available range.
    insertion_index = int(
        np.searchsorted(times, time)
    )

    start = insertion_index - (number_of_nodes + 1) // 2

    start = max(
        0,
        min(
            start,
            len(times) - number_of_nodes,
        ),
    )

    stop = start + number_of_nodes

    time_window = times[start:stop]
    table = positions[start:stop].copy()

    # Neville's algorithm. Each pass raises the polynomial order by
    # one and shortens the table by one entry.
    for degree in range(1, number_of_nodes):
        for index in range(number_of_nodes - degree):
            table[index] = (
                (time - time_window[index + degree]) * table[index]
                + (time_window[index] - time) * table[index + 1]
            ) / (
                time_window[index] - time_window[index + degree]
            )

    result = table[0]

    if positions.ndim == 1:
        return float(result)

    return result
