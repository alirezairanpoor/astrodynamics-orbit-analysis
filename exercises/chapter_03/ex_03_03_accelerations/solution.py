import numpy as np

from satellite_orbits.constants import MU_EARTH, R_EARTH

from satellite_orbits.astrodynamics.perturbations.gravity_third_body_balance_radius import (
    gravity_third_body_balance_radius,
)
from satellite_orbits.astrodynamics.perturbations.drag_srp_balance_altitude import (
    drag_srp_balance_altitude,
)
from satellite_orbits.astrodynamics.perturbations.interpolate_density import (
    interpolate_density,
)


# ============================================================
# GRAVITY COEFFICIENTS — NORMALIZED, BOOK TABLE 3.3
# ============================================================

gravity_terms = [
    ("J20", -484.165368e-6, 0.0),
    ("J22", 2.439261e-6, -1.400266e-6),
]


# ============================================================
# THIRD-BODY DATA
# ============================================================

third_bodies = [
    ("Moon", 4902.801, 384400.0),
    ("Sun", 1.32712440018e11, 149597870.691),
]


# ============================================================
# GRAVITATIONAL BALANCE ALTITUDES
# ============================================================

print("Exercise 3.3 — Accelerations")
print()
print("Gravity term | Third body | Balance altitude [km]")
print("-------------------------------------------------")

for term_name, C_normalized, S_normalized in gravity_terms:
    for body_name, mu_body, distance_body in third_bodies:

        r_balance = gravity_third_body_balance_radius(
            n=2,
            C_normalized=C_normalized,
            S_normalized=S_normalized,
            mu_earth=MU_EARTH,
            earth_radius=R_EARTH,
            mu_third_body=mu_body,
            third_body_distance=distance_body,
        )

        h_balance = r_balance - R_EARTH

        print(
            f"{term_name:12s} | "
            f"{body_name:10s} | "
            f"{h_balance:21.3f}"
        )


# ============================================================
# HARRIS–PRIESTER DENSITY ENVELOPES — BOOK TABLE 3.8
# ============================================================

# Mean solar activity.
# Only the 500–1000 km portion needed here is included.
#
# Columns:
# altitude [km], minimum density [g/km^3], maximum density [g/km^3]

density_table = np.array([
    [500, 0.3916,   2.042],
    [520, 0.2819,   1.605],
    [540, 0.2042,   1.267],
    [560, 0.1488,   1.005],
    [580, 0.1092,   0.7997],
    [600, 0.08070,  0.6390],
    [620, 0.06012,  0.5123],
    [640, 0.04519,  0.4121],
    [660, 0.03430,  0.3325],
    [680, 0.02632,  0.2691],
    [700, 0.02043,  0.2185],
    [720, 0.01607,  0.1779],
    [740, 0.01281,  0.1452],
    [760, 0.01036,  0.1190],
    [780, 0.008496, 0.09776],
    [800, 0.007069, 0.08059],
    [840, 0.004680, 0.05741],
    [880, 0.003200, 0.04210],
    [920, 0.002210, 0.03130],
    [960, 0.001560, 0.02360],
    [1000, 0.001150, 0.01810],
])

heights = density_table[:, 0]

# 1 g/km^3 = 1e-12 kg/m^3
densities_min = density_table[:, 1] * 1e-12
densities_max = density_table[:, 2] * 1e-12


def density_min(h: float) -> float:
    return interpolate_density(
        h,
        heights,
        densities_min,
    )


def density_max(h: float) -> float:
    return interpolate_density(
        h,
        heights,
        densities_max,
    )


# ============================================================
# DRAG / SOLAR RADIATION PRESSURE BALANCE
# ============================================================

C_D = 2.3
C_R = 1.3
solar_pressure = 4.560e-6  # N/m^2, at 1 AU

print()
print("Drag / SRP balance")
print("Harris–Priester envelopes, mean solar activity")
print("-----------------------------------------------")

for label, density_function in (
    ("Minimum-density envelope", density_min),
    ("Maximum-density envelope", density_max),
):

    h_balance = drag_srp_balance_altitude(
        density_function=density_function,
        h_lower=500.0,
        h_upper=1000.0,
        C_D=C_D,
        C_R=C_R,
        solar_pressure=solar_pressure,
        mu=MU_EARTH,
        earth_radius=R_EARTH,
    )

    print(f"{label}: {h_balance:.3f} km")

print()
print("Book reference: approximately 752 km.")
print("The balance altitude depends on the density model and conditions.")
