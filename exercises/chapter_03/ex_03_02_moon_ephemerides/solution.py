import numpy as np

from satellite_orbits.time_systems import calendar_to_mjd

from satellite_orbits.ephemerides.moon_position_low_precision import (
    moon_position_low_precision,
)
from satellite_orbits.ephemerides.chebyshev_position import (
    chebyshev_position,
)


# ============================================================
# COEFFICIENT INTERVAL
# ============================================================

# Exercise epochs are expressed in TT.
mjd_start = calendar_to_mjd(2006, 3, 14)
mjd_end = calendar_to_mjd(2006, 3, 18)


# ============================================================
# DE405 CHEBYSHEV COEFFICIENTS FROM EXERCISE 3.2
# ============================================================

# Each row: [ax_k, ay_k, az_k]
# Units: km
# Origin: Earth center
# Frame: equatorial J2000

coefficients = np.array([
    [-0.383089044877e6, -0.379891721705e5, -0.178496690739e5],
    [ 0.218158411755e5, -0.143611643157e6, -0.788257550332e5],
    [ 0.179067292901e5,  0.187126702787e4,  0.880684692614e3],
    [-0.836928063412e2,  0.112734362473e4,  0.618395886330e3],
    [-0.628266733052e2,  0.932891213817e0,  0.103331218595e1],
    [-0.459274434235e0, -0.191932684131e1, -0.104949867328e1],
    [ 0.491167202820e-1, -0.266517663332e-1, -0.150337371963e-1],
    [ 0.770804039288e-3,  0.104558913449e-2,  0.569056416308e-3],
    [-0.125935992206e-3, -0.359077689124e-4, -0.186297523287e-4],
    [ 0.500271026611e-5, -0.123405162037e-4, -0.680012420654e-5],
    [ 0.107044869186e-5,  0.180479239596e-6,  0.902057208454e-7],
    [ 0.172472464344e-8,  0.525522632334e-7,  0.287891446432e-7],
    [-0.269667589577e-8,  0.543313967009e-9,  0.319822827700e-9],
])


# ============================================================
# COMPUTE POSITIONS
# ============================================================

days = range(14, 19)

positions_analytical = []
positions_chebyshev = []

for day in days:

    mjd_tt = calendar_to_mjd(2006, 3, day)

    r_analytical_eci = moon_position_low_precision(mjd_tt)

    r_chebyshev_eci = chebyshev_position(
        mjd=mjd_tt,
        mjd_start=mjd_start,
        mjd_end=mjd_end,
        coefficients=coefficients,
    )

    positions_analytical.append(r_analytical_eci)
    positions_chebyshev.append(r_chebyshev_eci)

positions_analytical = np.array(positions_analytical)
positions_chebyshev = np.array(positions_chebyshev)


# ============================================================
# VALIDATE CHEBYSHEV RESULTS AGAINST THE BOOK
# ============================================================

book_positions_chebyshev = np.array([
    [-386976.783,  106369.219,  61240.442],
    [-403002.331,   34008.826,  21741.255],
    [-401058.650,  -39859.480, -18729.305],
    [-381019.563, -111859.423, -58322.341],
    [-343513.403, -178603.217, -95176.374],
])

# The book prints positions to 0.001 km.
np.testing.assert_allclose(
    positions_chebyshev,
    book_positions_chebyshev,
    rtol=0.0,
    atol=0.001,
)


# ============================================================
# POSITION DIFFERENCES
# ============================================================

delta_r_eci = (
    positions_analytical
    - positions_chebyshev
)

position_differences = np.linalg.norm(
    delta_r_eci,
    axis=1,
)


# ============================================================
# OUTPUT
# ============================================================

print("Exercise 3.2 — Moon Ephemerides")
print("Chebyshev book comparison: passed")

for label, positions in (
    ("Analytical model", positions_analytical),
    ("DE405 Chebyshev", positions_chebyshev),
):
    print()
    print(label)
    print("Date (TT)  |       x [km] |       y [km] |       z [km]")
    print("----------------------------------------------------------")

    for day, position in zip(days, positions):
        print(
            f"2006-03-{day:02d} | "
            f"{position[0]:12.3f} | "
            f"{position[1]:12.3f} | "
            f"{position[2]:12.3f}"
        )

print()
print("Position difference: analytical minus Chebyshev")
print("Date (TT)  | Norm [km]")
print("------------------------")

for day, difference in zip(days, position_differences):
    print(
        f"2006-03-{day:02d} | "
        f"{difference:9.3f}"
    )

print()
print(
    "Maximum difference at the five sampled epochs: "
    f"{np.max(position_differences):.3f} km"
)
