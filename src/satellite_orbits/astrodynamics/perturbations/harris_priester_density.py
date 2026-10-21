import numpy as np

from .interpolate_density import interpolate_density


# Columns: altitude [km], minimum density, maximum density.
# Density values in the source table are in g/km^3.
_HP_TABLE = np.array([
    [100, 497400, 497400],
    [120, 24900, 24900],
    [130, 8377, 8710],
    [140, 3899, 4059],
    [150, 2122, 2215],
    [160, 1263, 1344],
    [170, 800.8, 875.8],
    [180, 528.3, 601.0],
    [190, 361.7, 429.7],
    [200, 255.7, 316.2],
    [210, 183.9, 239.6],
    [220, 134.1, 185.3],
    [230, 99.49, 145.5],
    [240, 74.88, 115.7],
    [250, 57.09, 93.08],
    [260, 44.03, 75.55],
    [270, 34.30, 61.82],
    [280, 26.97, 50.95],
    [290, 21.39, 42.26],
    [300, 17.08, 35.26],
    [320, 10.99, 25.11],
    [340, 7.214, 18.19],
    [360, 4.824, 13.37],
    [380, 3.274, 9.955],
    [400, 2.249, 7.492],
    [420, 1.558, 5.684],
    [440, 1.091, 4.355],
    [460, 0.7701, 3.362],
    [480, 0.5474, 2.612],
    [500, 0.3916, 2.042],
    [520, 0.2819, 1.605],
    [540, 0.2042, 1.267],
    [560, 0.1488, 1.005],
    [580, 0.1092, 0.7997],
    [600, 0.08070, 0.6390],
    [620, 0.06012, 0.5123],
    [640, 0.04519, 0.4121],
    [660, 0.03430, 0.3325],
    [680, 0.02632, 0.2691],
    [700, 0.02043, 0.2185],
    [720, 0.01607, 0.1779],
    [740, 0.01281, 0.1452],
    [760, 0.01036, 0.1190],
    [780, 0.008496, 0.09776],
    [800, 0.007069, 0.08059],
    [840, 0.004680, 0.05741],
    [880, 0.003200, 0.04210],
    [920, 0.002210, 0.03130],
    [960, 0.001560, 0.02360],
    [1000, 0.001150, 0.01810],
], dtype=float)

_HEIGHTS = _HP_TABLE[:, 0]

# Conversion: g/km^3 -> kg/m^3.
_MIN_DENSITIES = _HP_TABLE[:, 1] * 1e-12
_MAX_DENSITIES = _HP_TABLE[:, 2] * 1e-12


def harris_priester_density(
    h: float,
    r_eci: np.ndarray,
    r_sun_eci: np.ndarray,
    exponent: float,
    lag_angle: float = np.deg2rad(30.0),
) -> float:
    """Return atmospheric density [kg/m^3].

    Parameters
    ----------
    h : float
        Altitude above the reference ellipsoid [km].
        Supported range: 100 to 1000 km.
    r_eci : array_like, shape (3,)
        Geocentric satellite position in ECI [km].
    r_sun_eci : array_like, shape (3,)
        Geocentric Sun position in the same ECI frame [km].
    exponent : float
        Harris-Priester exponent, usually between 2 and 6.
    lag_angle : float
        Eastward angular offset of the density maximum
        from the Sun direction [rad].

    Notes
    -----
    Uses the mean-solar-activity table.
    The ECI z-axis approximates Earth's rotation axis.
    """
    r_eci = np.asarray(r_eci, dtype=float)
    r_sun_eci = np.asarray(r_sun_eci, dtype=float)

    if r_eci.shape != (3,) or r_sun_eci.shape != (3,):
        raise ValueError("Position vectors must have shape (3,).")

    if not (
        np.all(np.isfinite(r_eci))
        and np.all(np.isfinite(r_sun_eci))
        and np.isfinite(h)
        and np.isfinite(exponent)
        and np.isfinite(lag_angle)
    ):
        raise ValueError("Inputs must be finite.")

    if not _HEIGHTS[0] <= h <= _HEIGHTS[-1]:
        raise ValueError("Altitude must be between 100 and 1000 km.")

    if exponent <= 0:
        raise ValueError("Exponent must be positive.")

    r_norm = np.linalg.norm(r_eci)
    r_sun_norm = np.linalg.norm(r_sun_eci)

    if r_norm == 0 or r_sun_norm == 0:
        raise ValueError("Position vectors must be nonzero.")

    # Interpolate both density envelopes at the given altitude.
    rho_min = interpolate_density(h, _HEIGHTS, _MIN_DENSITIES)
    rho_max = interpolate_density(h, _HEIGHTS, _MAX_DENSITIES)

    satellite_direction = r_eci / r_norm
    sun_direction = r_sun_eci / r_sun_norm

    # Rotate the Sun direction eastward about the ECI z-axis.
    c = np.cos(lag_angle)
    s = np.sin(lag_angle)

    apex_direction = np.array([
        c * sun_direction[0] - s * sun_direction[1],
        s * sun_direction[0] + c * sun_direction[1],
        sun_direction[2],
    ])

    cos_psi = np.clip(
        np.dot(satellite_direction, apex_direction),
        -1.0,
        1.0,
    )

    weight = (0.5 * (1.0 + cos_psi)) ** (0.5 * exponent)

    return float(rho_min + (rho_max - rho_min) * weight)
