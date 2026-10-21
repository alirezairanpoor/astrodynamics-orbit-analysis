from math import radians


def mjd_to_greenwich_hour_angle(
    mjd: float,
) -> float:
    """
    Convert Modified Julian Date (MJD)
    to Greenwich Hour Angle.

    Parameters
    ----------
    mjd : float
        Modified Julian Date [days].

    Returns
    -------
    float
        Greenwich Hour Angle [rad].
    """

    J2000_MJD = 51544.5

    days_since_j2000 = (
        mjd - J2000_MJD
    )

    greenwich_hour_angle_deg = (
        280.4606
        + 360.9856473 * days_since_j2000
    ) % 360.0

    greenwich_hour_angle = radians(
        greenwich_hour_angle_deg
    )

    return greenwich_hour_angle
