from math import radians, degrees, pi

from satellite_orbits.constants import (
    MU_EARTH,
    R_EARTH,
    J2_EARTH,
)
from satellite_orbits.astrodynamics import (
    semi_major_axis_from_period,
    semi_major_axis_from_draconic_period,
    sun_synchronous_inclination,
    j2_secular_rates,
)


# ============================================================
# USER INPUTS
# ============================================================

# Repeat-ground-track requirement
K = 3       # repeat cycle [days]
N = 43      # number of draconic orbits

# Required sun-synchronous RAAN rate
RAAN_rate_target = radians(
    0.985647240
) / 86400.0  # rad/s


# ============================================================
# DRACONIC PERIOD
# ============================================================

TN = (
    K * 86400.0
    / N
)


# ============================================================
# ITERATION 0 — INITIAL GUESS
# ============================================================

# Initially ignore J2 corrections:
# TN ≈ Keplerian orbital period
a = semi_major_axis_from_period(
    TN,
    MU_EARTH,
)

# SSO inclination corresponding to this a
i = sun_synchronous_inclination(
    a,
    RAAN_rate_target,
    MU_EARTH,
    R_EARTH,
    J2_EARTH,
)


print("Exercise 2.5 — Sun-Synchronous Repeat Orbit")
print("------------------------------------------------")

print(f"Draconic period TN : {TN:.3f} s")
print()

print("Iteration 0")
print(f"a : {a:.3f} km")
print(f"h : {a - R_EARTH:.3f} km")
print(f"i : {degrees(i):.3f} deg")


# ============================================================
# REFINEMENT ITERATIONS
# ============================================================

number_of_iterations = 2

for iteration in range(
    1,
    number_of_iterations + 1,
):

    # --------------------------------------------------------
    # 1. J2 secular rates for current a and i
    # --------------------------------------------------------

    (
        n0,
        RAAN_rate,
        AOP_rate,
        delta_n,
    ) = j2_secular_rates(
        a,
        i,
        MU_EARTH,
        R_EARTH,
        J2_EARTH,
    )

    # --------------------------------------------------------
    # 2. Refine a using the draconic period
    # --------------------------------------------------------

    a = semi_major_axis_from_draconic_period(
        TN,
        delta_n,
        AOP_rate,
        MU_EARTH,
    )

    # --------------------------------------------------------
    # 3. Recalculate SSO inclination using new a
    # --------------------------------------------------------

    i = sun_synchronous_inclination(
        a,
        RAAN_rate_target,
        MU_EARTH,
        R_EARTH,
        J2_EARTH,
    )

    # --------------------------------------------------------
    # OUTPUT
    # --------------------------------------------------------

    print()
    print(f"Iteration {iteration}")

    print(
        f"AOP_rate : "
        f"{degrees(AOP_rate) * 86400:.5f} deg/day"
    )

    print(
        f"delta_n  : "
        f"{degrees(delta_n) * 86400:.5f} deg/day"
    )

    print(
        f"a        : "
        f"{a:.3f} km"
    )

    print(
        f"h        : "
        f"{a - R_EARTH:.3f} km"
    )

    print(
        f"i        : "
        f"{degrees(i):.3f} deg"
    )


# ============================================================
# FINAL RESULT
# ============================================================

print()
print("Final Result")
print("------------")

print(
    f"a = {a:.3f} km"
)

print(
    f"h = {a - R_EARTH:.3f} km"
)

print(
    f"i = {degrees(i):.3f} deg"
)
# ============================================================
# REPEAT GROUND TRACK VERIFICATION
# ============================================================

# Greenwich Hour Angle rate from Chapter 2
greenwich_hour_angle_rate = radians(
    360.9856473
) / 86400.0  # rad/s


print()
print("Repeat Ground Track Verification")
print("--------------------------------")
print("Crossing | Longitude of Ascending Node [deg]")
print("--------------------------------------------")


for crossing_number in range(N + 1):

    elapsed_time = (
        crossing_number * TN
    )

    lambda_RAAN = (
        (
            RAAN_rate_target
            - greenwich_hour_angle_rate
        )
        * elapsed_time
    )

    lambda_RAAN_deg = (
        degrees(lambda_RAAN)
        % 360.0
    )

    # Avoid printing 360 deg instead of 0 deg
    if abs(lambda_RAAN_deg - 360.0) < 1e-6:
        lambda_RAAN_deg = 0.0

    print(
        f"{crossing_number:8d} | "
        f"{lambda_RAAN_deg:10.3f}"
    )
