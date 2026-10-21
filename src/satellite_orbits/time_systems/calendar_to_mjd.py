from math import floor


def calendar_to_mjd(
    year: int,
    month: int,
    day: int,
    hour: int = 0,
    minute: int = 0,
    second: float = 0.0,
) -> float:
    """
    Convert a Gregorian calendar date and time
    to Modified Julian Date (MJD).

    Parameters
    ----------
    year : int
        Calendar year.

    month : int
        Calendar month.

    day : int
        Calendar day.

    hour : int
        UTC hour.

    minute : int
        UTC minute.

    second : float
        UTC second.

    Returns
    -------
    float
        Modified Julian Date.
    """

    if month <= 2:
        y = year - 1
        m = month + 12
    else:
        y = year
        m = month

    B = (
        floor(y / 400)
        - floor(y / 100)
        + floor(y / 4)
    )

    fraction_of_day = (
        hour
        + minute / 60
        + second / 3600
    ) / 24

    D = (
        day
        + fraction_of_day
    )

    mjd = (
        365 * y
        - 679004
        + B
        + floor(
            30.6001
            * (m + 1)
        )
        + D
    )

    return mjd
