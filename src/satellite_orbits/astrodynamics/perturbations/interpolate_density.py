import numpy as np


def interpolate_density(
    h: float,
    heights: np.ndarray,
    densities: np.ndarray,
) -> float:
    """
    Exponentially interpolate an atmospheric density table.

    Inputs:
        h: altitude [km]
        heights: strictly increasing table altitudes [km]
        densities: positive table densities [kg/m^3]

    Output:
        Density at h [kg/m^3]

    Extrapolation outside the table is not allowed.
    """

    heights = np.asarray(heights, dtype=float)
    densities = np.asarray(densities, dtype=float)

    if (
        heights.ndim != 1
        or densities.ndim != 1
        or len(heights) != len(densities)
        or len(heights) < 2
    ):
        raise ValueError("Provide matching one-dimensional tables.")

    if np.any(np.diff(heights) <= 0.0):
        raise ValueError("Heights must be strictly increasing.")

    if np.any(densities <= 0.0):
        raise ValueError("Densities must be positive.")

    if not heights[0] <= h <= heights[-1]:
        raise ValueError("Altitude is outside the density table.")

    log_density = np.interp(
        h,
        heights,
        np.log(densities),
    )

    return float(np.exp(log_density))