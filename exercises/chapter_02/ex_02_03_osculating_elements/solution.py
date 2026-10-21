from math import degrees

import numpy as np

from satellite_orbits.constants import MU_EARTH
from satellite_orbits.astrodynamics import states_to_elements


r_eci = np.array([
    10000.0,
    40000.0,
    -5000.0,
])

v_eci = np.array([
    -1.5,
    1.0,
    -0.1,
])


a, e, i, RAAN, AOP, M = states_to_elements(
    r_eci,
    v_eci,
    MU_EARTH,
)


print("Exercise 2.3 — Osculating Elements")
print("-----------------------------------")

print(f"a    : {a:.3f} km")
print(f"e    : {e:.7f}")
print(f"i    : {degrees(i):.3f} deg")
print(f"RAAN : {degrees(RAAN):.3f} deg")
print(f"AOP  : {degrees(AOP):.3f} deg")
print(f"M    : {degrees(M):.3f} deg")
