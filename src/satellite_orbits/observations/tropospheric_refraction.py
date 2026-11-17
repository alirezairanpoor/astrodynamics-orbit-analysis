from math import exp, tan


def water_vapour_pressure(
    temperature: float,
    relative_humidity: float,
) -> float:
    """
    Compute the partial pressure of water vapour.

    The saturation pressure grows exponentially with temperature,
    which is why the wet refractivity term rises steeply in warm
    air even at constant relative humidity.

    Source: equation 6.55, after Grosskopf (1970).

    Parameters
    ----------
    temperature : float
        Air temperature [K].

    relative_humidity : float
        Relative humidity, between 0 and 1 [-].

    Returns
    -------
    float
        Partial pressure of water vapour [hPa].
    """

    temperature_celsius = temperature - 273.15

    pressure = (
        6.10
        * relative_humidity
        * exp(
            17.15 * temperature_celsius
            / (234.7 + temperature_celsius)
        )
    )

    return pressure


def tropospheric_refractivity(
    pressure: float,
    temperature: float,
    relative_humidity: float,
) -> tuple[float, float]:
    """
    Compute the dry and wet tropospheric refractivity.

    The troposphere is non-dispersive for radio waves, so the
    result does not depend on the carrier frequency.

    The dry term falls with temperature while the wet term rises
    steeply, so the net temperature dependence is not obvious in
    advance and has to be evaluated.

    Source: equations 6.53 and 6.54, after Goad and Goodman (1974).

    Parameters
    ----------
    pressure : float
        Atmospheric pressure [hPa].

    temperature : float
        Air temperature [K].

    relative_humidity : float
        Relative humidity, between 0 and 1 [-].

    Returns
    -------
    tuple
        dry_refractivity : Dry component N1 [-].
        wet_refractivity : Wet component N2 [-].
    """

    vapour_pressure = water_vapour_pressure(
        temperature,
        relative_humidity,
    )

    dry_refractivity = (
        77.624
        * pressure
        / temperature
    )

    wet_refractivity = (
        371900.0
        * vapour_pressure
        / temperature ** 2
        - 12.92
        * vapour_pressure
        / temperature
    )

    return (
        dry_refractivity,
        wet_refractivity,
    )


def elevation_refraction_correction(
    elevation: float,
    pressure: float,
    temperature: float,
    relative_humidity: float,
) -> float:
    """
    Compute the tropospheric bending of the elevation angle.

    A plane atmosphere of constant refractivity is assumed. Snell's
    law then gives the bending angle directly in terms of the
    refractive index and the elevation.

    The correction vanishes at the zenith, where the signal crosses
    the layers at normal incidence, and diverges towards the
    horizon, where the plane-atmosphere assumption breaks down. The
    model should not be trusted below a few degrees of elevation.

    Azimuth measurements are not affected by refraction.

    Source: equation 6.50.

    Parameters
    ----------
    elevation : float
        Geometric elevation of the satellite [rad].

    pressure : float
        Atmospheric pressure [hPa].

    temperature : float
        Air temperature [K].

    relative_humidity : float
        Relative humidity, between 0 and 1 [-].

    Returns
    -------
    float
        Elevation correction to be added to the geometric
        elevation [rad].
    """

    dry_refractivity, wet_refractivity = (
        tropospheric_refractivity(
            pressure,
            temperature,
            relative_humidity,
        )
    )

    refractive_index_excess = (
        dry_refractivity
        + wet_refractivity
    ) * 1.0e-6

    correction = (
        refractive_index_excess
        / tan(elevation)
    )

    return correction
