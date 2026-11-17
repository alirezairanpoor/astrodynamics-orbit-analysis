import numpy as np

from satellite_orbits.constants import (
    C_LIGHT,
    GPS_L1_FREQUENCY,
    GPS_L2_FREQUENCY,
    OMEGA_EARTH,
)
from satellite_orbits.ephemerides import (
    interpolate_tabulated_position,
)
from satellite_orbits.observations import (
    ionosphere_free_pseudorange,
    light_time_earth_fixed,
)


def main() -> None:
    # Inputs: PRN 1 precise ephemeris from the IGS, 1998/02/19.
    # Source: Montenbruck and Gill, Exercise 6.3, pages 230-231.
    #
    # Positions are converted to km, the length unit used
    # throughout this project. Satellite clock offsets are in
    # microseconds.
    ephemeris_hours = np.array([
        8.00, 8.25, 8.50, 8.75,
        9.00, 9.25, 9.50, 9.75,
        10.00, 10.25, 10.50, 10.75,
        11.00,
    ])

    ephemeris_positions_km = np.array([
        [-15504.291797, -21530.763883, -1271.498273],
        [-15284.290679, -21684.703684, 1573.435406],
        [-14871.711829, -21600.510259, 4391.350089],
        [-14242.843546, -21306.712708, 7133.948741],
        [-13380.818523, -20837.175663, 9754.366309],
        [-12276.418004, -20229.688085, 12207.953668],
        [-10928.585710, -19524.421024, 14453.015617],
        [-9344.633744, -18762.314034, 16451.492281],
        [-7540.134384, -17983.451817, 18169.574686],
        [-5538.503062, -17225.491970, 19578.246580],
        [-3370.289205, -16522.202377, 20653.745961],
        [-1072.201838, -15902.162018, 21377.940941],
        [1314.093678, -15387.672739, 21738.615794],
    ])

    ephemeris_clock_microseconds = np.array([
        40.018233, 40.097295, 40.028697, 40.154941,
        40.193626, 40.039288, 40.012677, 39.883106,
        40.181357, 40.328261, 40.039533, 40.052642,
        40.025493,
    ])

    # Pseudorange measurements: P code at L2 and C/A code at L1.
    # Converted to km to match the ephemeris.
    measurement_hours = np.array([
        8.50, 9.00, 9.50, 10.00, 10.50, 11.00,
    ])

    pseudorange_p2_km = np.array([
        21096577.475, 20519964.850, 20282706.954,
        20375838.496, 20751678.769, 21340055.129,
    ]) / 1000.0

    pseudorange_c1_km = np.array([
        21096579.501, 20519966.875, 20282709.233,
        20375840.613, 20751680.997, 21340057.362,
    ]) / 1000.0

    # Receiver position at Goldstone in WGS84 [km].
    receiver_position_km = np.array([
        -2353.614128,
        -4641.385447,
        3676.976501,
    ])

    # Interpolation order for the 15 min ephemeris grid.
    interpolation_order = 9

    ephemeris_seconds = ephemeris_hours * 3600.0
    measurement_seconds = measurement_hours * 3600.0

    def satellite_position(t):
        return interpolate_tabulated_position(
            ephemeris_seconds,
            ephemeris_positions_km,
            t,
            interpolation_order,
        )

    # Reference output: mean clock bias printed in the textbook.
    book_clock_bias_m = 17674.77
    book_clock_errors_m = (
        17674.52, 17676.14, 17675.59,
        17674.45, 17673.53, 17674.41,
    )

    clock_errors_m = []

    print("Exercise 6.3 - User Clock Error from GPS Pseudorange")
    print("Receiver: Rogue at Goldstone, WGS84")
    print("Satellite: PRN 1, IGS precise ephemeris")
    print(f"Interpolation order: {interpolation_order}")
    print()
    print(
        "GPS       Iono-free      Distance   SatClock"
        "   ClockError"
    )
    print(
        "hh:mm:ss        [m]           [m]        [m]"
        "          [m]"
    )
    print("-" * 66)

    for index, reception_time in enumerate(measurement_seconds):

        # Remove the first-order ionospheric delay before any
        # geometric modelling.
        pseudorange_km = ionosphere_free_pseudorange(
            pseudorange_c1_km[index],
            pseudorange_p2_km[index],
            GPS_L1_FREQUENCY,
            GPS_L2_FREQUENCY,
        )

        # Rigorous light-time solution in the Earth-fixed frame.
        _, distance_km = light_time_earth_fixed(
            satellite_position,
            receiver_position_km,
            reception_time,
            C_LIGHT,
            OMEGA_EARTH,
        )

        # The measurement epochs coincide with ephemeris epochs,
        # so the tabulated satellite clock offset is used without
        # interpolation, as the exercise prescribes.
        ephemeris_index = int(
            np.argmin(
                np.abs(ephemeris_seconds - reception_time)
            )
        )

        satellite_clock_km = (
            C_LIGHT
            * ephemeris_clock_microseconds[ephemeris_index]
            * 1.0e-6
        )

        # Solve the pseudorange equation for the receiver clock
        # offset, expressed as a range.
        clock_error_km = (
            pseudorange_km
            - distance_km
            + satellite_clock_km
        )

        clock_errors_m.append(clock_error_km * 1000.0)

        hours = int(measurement_hours[index])
        minutes = int(
            round((measurement_hours[index] - hours) * 60.0)
        )

        print(
            f"{hours:02d}:{minutes:02d}:00.0 "
            f"{pseudorange_km * 1000.0:14.3f} "
            f"{distance_km * 1000.0:13.3f} "
            f"{satellite_clock_km * 1000.0:10.3f} "
            f"{clock_error_km * 1000.0:12.2f}"
        )

    clock_errors_m = np.array(clock_errors_m)
    mean_clock_bias_m = float(np.mean(clock_errors_m))
    deviations_m = clock_errors_m - mean_clock_bias_m

    print()
    print("GPS       ClockError      Delta")
    print("hh:mm:ss         [m]        [m]")
    print("-" * 33)

    for index in range(len(measurement_hours)):
        hours = int(measurement_hours[index])
        minutes = int(
            round((measurement_hours[index] - hours) * 60.0)
        )
        print(
            f"{hours:02d}:{minutes:02d}:00.0 "
            f"{clock_errors_m[index]:12.2f} "
            f"{deviations_m[index]:10.2f}"
        )

    print()
    print(f"Mean clock bias [m]: {mean_clock_bias_m:.2f}")
    print(f"Textbook        [m]: {book_clock_bias_m:.2f}")
    print(
        f"Difference      [m]: "
        f"{mean_clock_bias_m - book_clock_bias_m:+.2f}"
    )
    print(
        f"Maximum spread  [m]: "
        f"{np.max(clock_errors_m) - np.min(clock_errors_m):.2f}"
    )

    # Verification against the textbook values.
    np.testing.assert_allclose(
        clock_errors_m,
        book_clock_errors_m,
        rtol=0.0,
        atol=0.5,
    )

    if abs(mean_clock_bias_m - book_clock_bias_m) > 0.5:
        raise AssertionError("Mean clock bias check failed.")

    # The spread of the residuals reflects the accuracy of the
    # ionosphere-free pseudoranges, not a drift of the receiver
    # clock, which is negligible over 2.5 h.
    if np.max(np.abs(deviations_m)) > 3.0:
        raise AssertionError("Residual spread check failed.")

    print()
    print("Textbook clock error comparison: passed")
    print("Mean clock bias comparison: passed")
    print("Residual spread bound: passed")


if __name__ == "__main__":
    main()
