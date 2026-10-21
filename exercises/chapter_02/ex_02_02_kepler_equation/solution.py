from math import radians, degrees

from satellite_orbits.astrodynamics.keplerian.kepler import (
    solve_kepler_newton,
    solve_kepler_fixed_point,
)


e = 0.72

M_values_deg = [
    4.0,
    50.0,
]


for M_deg in M_values_deg:

    M = radians(M_deg)

    E_newton, iterations_newton, trig_newton = solve_kepler_newton(
        M,
        e,
    )

    E_fixed, iterations_fixed, trig_fixed = solve_kepler_fixed_point(
        M,
        e,
    )

    print(f"\nM = {M_deg:.1f} deg")
    print("----------------------------------------")

    print("\nNewton method")
    print(f"E                 : {degrees(E_newton):.10f} deg")
    print(f"Iterations        : {iterations_newton}")
    print(f"Trig evaluations  : {trig_newton}")

    print("\nFixed-point method")
    print(f"E                 : {degrees(E_fixed):.10f} deg")
    print(f"Iterations        : {iterations_fixed}")
    print(f"Trig evaluations  : {trig_fixed}")
