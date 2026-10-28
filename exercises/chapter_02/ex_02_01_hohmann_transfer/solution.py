from satellite_orbits.constants import MU_EARTH, R_EARTH
from satellite_orbits.astrodynamics.keplerian.two_body import circular_velocity
from satellite_orbits.astrodynamics.orbit_design.maneuvers import hohmann_transfer


# Exercise 2.1
# Orbit raising from 750 km to 775 km altitude

h1 = 750.0  # km
h2 = 775.0  # km

# Convert altitude above Earth to orbital radius
r1 = R_EARTH + h1
r2 = R_EARTH + h2

# Exact Hohmann transfer
delta_v1, delta_v2, delta_v_total = hohmann_transfer(
    MU_EARTH,
    r1,
    r2,
)

# Circular velocities before and after the maneuver
v1 = circular_velocity(MU_EARTH, r1)
v2 = circular_velocity(MU_EARTH, r2)

# First-order approximation from the exercise
delta_v_approx = abs(v1 - v2)


print("Exercise 2.1 — Hohmann Orbit Raising")
print("-------------------------------------")

print(f"Initial orbital radius : {r1:.3f} km")
print(f"Final orbital radius   : {r2:.3f} km")

print(f"\nInitial circular velocity : {v1:.6f} km/s")
print(f"Final circular velocity   : {v2:.6f} km/s")

print(f"\nFirst burn  Δv1 : {delta_v1 * 1000:.3f} m/s")
print(f"Second burn Δv2 : {delta_v2 * 1000:.3f} m/s")
print(f"Total Δv         : {delta_v_total * 1000:.3f} m/s")

print(
    f"First-order approximation |v1 - v2| : "
    f"{delta_v_approx * 1000:.3f} m/s"
)