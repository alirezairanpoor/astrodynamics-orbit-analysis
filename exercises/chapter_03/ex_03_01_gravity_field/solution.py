from time import process_time

import numpy as np

from satellite_orbits.constants import (
    MU_EARTH,
    R_EARTH,
    J2_EARTH,
)
from satellite_orbits.astrodynamics.perturbations.gravity_acceleration_cunningham import (
    gravity_acceleration_cunningham,
)


# ============================================================
# USER INPUTS
# ============================================================

# Fixed evaluation position in ECEF [km]
# Chosen for this timing experiment.
r_ecef = np.array([
    7000.0,
    1200.0,
    2300.0,
])

maximum_degree = 20
number_of_evaluations = 10_000


# ============================================================
# COEFFICIENTS FOR THE TIMING EXPERIMENT
# ============================================================

# Unnormalized coefficients.
# This is a central-gravity + J2 model, not a full degree-20
# Earth gravity dataset.
#
# The acceleration function still evaluates every requested
# term, including those whose coefficients are zero.

C = np.zeros((maximum_degree + 1, maximum_degree + 1))
S = np.zeros_like(C)

C[0, 0] = 1.0
C[2, 0] = -J2_EARTH


# ============================================================
# CENTRAL-GRAVITY CHECK
# ============================================================

a_cunningham = gravity_acceleration_cunningham(
    r_ecef,
    C,
    S,
    0,
    0,
    MU_EARTH,
    R_EARTH,
)

r = np.linalg.norm(r_ecef)

a_reference = -MU_EARTH * r_ecef / r**3

np.testing.assert_allclose(
    a_cunningham,
    a_reference,
    rtol=1e-12,
    atol=1e-15,
)

print("Exercise 3.1 — Gravity Field")
print("--------------------------------------------")
print("Central-gravity check: passed")


# ============================================================
# CPU TIMING
# ============================================================

degrees = np.arange(1, maximum_degree + 1)
cpu_times = []

print()
print(f"Evaluations per model: {number_of_evaluations}")
print()
print("Degree/order | CPU time [s] | Time per call [us]")
print("------------------------------------------------")

for degree in degrees:

    n_max = int(degree)
    m_max = int(degree)

    # Warm-up call outside the timed interval.
    gravity_acceleration_cunningham(
        r_ecef,
        C,
        S,
        n_max,
        m_max,
        MU_EARTH,
        R_EARTH,
    )

    start_time = process_time()

    for _ in range(number_of_evaluations):
        gravity_acceleration_cunningham(
            r_ecef,
            C,
            S,
            n_max,
            m_max,
            MU_EARTH,
            R_EARTH,
        )

    elapsed_time = process_time() - start_time
    cpu_times.append(elapsed_time)

    time_per_call_us = (
        elapsed_time
        / number_of_evaluations
        * 1e6
    )

    print(
        f"{n_max:6d}/{m_max:<5d} | "
        f"{elapsed_time:12.4f} | "
        f"{time_per_call_us:18.2f}",
        flush=True,
    )


# ============================================================
# QUADRATIC FIT
# ============================================================

cpu_times = np.array(cpu_times)

# Fit: CPU time = b2 * N^2 + b1 * N + b0
b2, b1, b0 = np.polyfit(
    degrees,
    cpu_times,
    deg=2,
)

print()
print("Quadratic fit: t(N) = b2*N^2 + b1*N + b0")
print(f"b2 = {b2:.6e} s")
print(f"b1 = {b1:.6e} s")
print(f"b0 = {b0:.6e} s")


# ============================================================
# COMPARE 10x10 AND 4x4 MODELS
# ============================================================

time_4 = cpu_times[3]
time_10 = cpu_times[9]

print()
print(f"Measured t(10)/t(4): {time_10 / time_4:.3f}")
print(f"Pure N^2 prediction: {(10 / 4)**2:.3f}")
