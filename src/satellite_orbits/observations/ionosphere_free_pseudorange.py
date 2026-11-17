def ionosphere_free_pseudorange(
    pseudorange_f1: float,
    pseudorange_f2: float,
    frequency_f1: float,
    frequency_f2: float,
) -> float:
    """
    Remove the first-order ionospheric delay from dual-frequency
    pseudoranges.

    The ionospheric group delay scales as 1 / frequency^2, so two
    measurements at different frequencies allow the total electron
    content to be eliminated rather than modelled. The weights are
    fixed by requiring that the true range keeps a unit coefficient
    while the electron-content terms cancel.

    The weights are large and of opposite sign, so the measurement
    noise of the combination is about three times that of a single
    observable for the GPS L1 and L2 frequencies. This is the price
    paid for removing a systematic error of tens of metres.

    Only the first-order term is removed. Higher-order ionospheric
    terms remain, at the millimetre to centimetre level.

    Parameters
    ----------
    pseudorange_f1 : float
        Pseudorange measured at the first frequency [m].

    pseudorange_f2 : float
        Pseudorange measured at the second frequency [m].

    frequency_f1 : float
        First carrier frequency [Hz].

    frequency_f2 : float
        Second carrier frequency [Hz].

    Returns
    -------
    float
        Ionosphere-free pseudorange [m].
    """

    frequency_ratio_squared = (
        frequency_f1 / frequency_f2
    ) ** 2

    pseudorange = (
        pseudorange_f2
        - frequency_ratio_squared * pseudorange_f1
    ) / (
        1.0 - frequency_ratio_squared
    )

    return float(pseudorange)