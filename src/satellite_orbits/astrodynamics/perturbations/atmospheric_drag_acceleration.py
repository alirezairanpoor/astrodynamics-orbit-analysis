import numpy as np


def atmospheric_drag_acceleration(
    r_eci: np.ndarray,
    v_eci: np.ndarray,
    density: float,
    area: float,
    mass: float,
    C_D: float,
    earth_rotation_rate: float = 7.292115e-5,
) -> np.ndarray:
    """
    Compute atmospheric drag acceleration in ECI.

    Inputs:
        r_eci: Earth-to-satellite position [km]
        v_eci: satellite inertial velocity [km/s]
        density: local atmospheric density [kg/m^3]
        area: effective cross-sectional area for drag [m^2]
        mass: satellite mass [kg]
        C_D: drag coefficient
        earth_rotation_rate: Earth rotation rate [rad/s]

    Output:
        Drag acceleration in ECI [km/s^2].

    Assumptions:
        Atmosphere rotates rigidly with Earth.
        Winds are neglected.
        Earth's rotation axis is approximated by the ECI z-axis.
    """

    r_eci = np.asarray(r_eci, dtype=float)
    v_eci = np.asarray(v_eci, dtype=float)

    if r_eci.shape != (3,) or v_eci.shape != (3,):
        raise ValueError("State vectors must have shape (3,).")

    if mass <= 0.0:
        raise ValueError("Mass must be positive.")

    if density < 0.0 or area < 0.0 or C_D < 0.0:
        raise ValueError(
            "Density, area, and C_D must be nonnegative."
        )

    # --------------------------------------------------------
    # ATMOSPHERIC VELOCITY IN ECI
    # --------------------------------------------------------

    omega_earth_eci = np.array([
        0.0,
        0.0,
        earth_rotation_rate,
    ])

    v_atmosphere_eci = np.cross(
        omega_earth_eci,
        r_eci,
    )  # km/s

    # --------------------------------------------------------
    # VELOCITY RELATIVE TO THE ATMOSPHERE
    # --------------------------------------------------------

    v_relative_eci = v_eci - v_atmosphere_eci  # km/s

    # Convert to SI for the drag calculation.
    v_relative_m_s = v_relative_eci * 1000.0
    relative_speed_m_s = np.linalg.norm(v_relative_m_s)

    # --------------------------------------------------------
    # DRAG ACCELERATION
    # --------------------------------------------------------

    a_drag_m_s2 = (
        -0.5
        * C_D
        * density
        * area
        / mass
        * relative_speed_m_s
        * v_relative_m_s
    )

    return a_drag_m_s2 / 1000.0
