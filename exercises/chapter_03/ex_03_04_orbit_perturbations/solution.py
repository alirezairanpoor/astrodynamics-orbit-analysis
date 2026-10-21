import numpy as np

from satellite_orbits.time_systems.calendar_to_mjd import (
    calendar_to_mjd,
)
from satellite_orbits.time_systems.mjd_to_greenwich_hour_angle import (
    mjd_to_greenwich_hour_angle,
)
from satellite_orbits.astrodynamics.keplerian.elements_to_states import (
    elements_to_states,
)
from satellite_orbits.ephemerides.sun_position_low_precision import (
    sun_position_low_precision,
)
from satellite_orbits.ephemerides.moon_position_low_precision import (
    moon_position_low_precision,
)
from satellite_orbits.astrodynamics.perturbations.jgm3_gravity_model import (
    load_jgm3_gravity_model,
    JGM3_MU,
    JGM3_RADIUS,
)
from satellite_orbits.astrodynamics.perturbations.total_acceleration import (
    total_acceleration,
)
from satellite_orbits.astrodynamics.propagation.orbit_state_derivative import (
    orbit_state_derivative,
)
from satellite_orbits.astrodynamics.propagation.integrate_rk4 import (
    integrate_rk4,
)


SECONDS_PER_DAY = 86400.0

# Valid for the exercise dates in March 1999.
TT_MINUS_UTC = 64.184
UT1_MINUS_UTC = 0.0  # Approximation.

MU_SUN = 1.32712440018e11
MU_MOON = 4902.800066

# Numerical acceptance criterion chosen for this exercise.
CONVERGENCE_LIMIT_M = 0.1


def main():
    mjd_epoch_utc = calendar_to_mjd(
        1999, 3, 1, 0, 0, 0.0,
    )
    C, S = load_jgm3_gravity_model()

    # label, degree, order, omitted effect
    models = [
        ("2/0", 2, 0, None),
        ("2/2", 2, 2, None),
        ("4/4", 4, 4, None),
        ("10/10", 10, 10, None),
        ("No Sun", 20, 20, "sun"),
        ("No Moon", 20, 20, "moon"),
        ("No SRP", 20, 20, "srp"),
        ("No Drag", 20, 20, "drag"),
    ]

    satellites = [
        {
            "name": "LEO",
            "a": 7178.0,
            "e": 0.001,
            "i_deg": 98.57,
            "area": 5.0,
            "duration": SECONDS_PER_DAY,
            "step": 5.0,
            "drag": True,
        },
        {
            "name": "GEO",
            "a": 42166.0,
            "e": 0.0004,
            "i_deg": 0.02,
            "area": 10.0,
            "duration": 2.0 * SECONDS_PER_DAY,
            "step": 30.0,
            "drag": False,
        },
    ]

    results = {}

    for satellite in satellites:
        name = satellite["name"]

        a = satellite["a"]
        e = satellite["e"]
        i = np.deg2rad(satellite["i_deg"])
        RAAN = 0.0
        AOP = 0.0
        M = 0.0

        r0_eci, v0_eci = elements_to_states(
            a, e, i, RAAN, AOP, M, JGM3_MU,
        )
        state0 = np.concatenate((r0_eci, v0_eci))

        period = 2.0 * np.pi * np.sqrt(a**3 / JGM3_MU)
        duration = satellite["duration"]
        n_steps = int(np.ceil(duration / satellite["step"]))

        if name == "LEO":
            intervals = [period, SECONDS_PER_DAY]
        else:
            intervals = [SECONDS_PER_DAY, 2.0 * SECONDS_PER_DAY]

        def propagate(
            n_max=20,
            m_max=20,
            omitted=None,
            steps=n_steps,
        ):
            include_drag = (
                satellite["drag"] and omitted != "drag"
            )

            def acceleration_function(t, r_eci, v_eci):
                mjd_utc = mjd_epoch_utc + t / SECONDS_PER_DAY
                mjd_tt = mjd_utc + TT_MINUS_UTC / SECONDS_PER_DAY
                mjd_ut1 = mjd_utc + UT1_MINUS_UTC / SECONDS_PER_DAY

                r_sun_eci = sun_position_low_precision(mjd_tt)
                r_moon_eci = moon_position_low_precision(mjd_tt)

                greenwich_hour_angle = (
                    mjd_to_greenwich_hour_angle(mjd_ut1)
                )

                return total_acceleration(
                    r_eci,
                    v_eci,
                    greenwich_hour_angle=greenwich_hour_angle,
                    r_sun_eci=r_sun_eci,
                    r_moon_eci=r_moon_eci,
                    C=C,
                    S=S,
                    n_max=n_max,
                    m_max=m_max,
                    mu_earth=JGM3_MU,
                    gravity_reference_radius=JGM3_RADIUS,
                    mu_sun=MU_SUN,
                    mu_moon=MU_MOON,
                    area=satellite["area"],
                    mass=1000.0,
                    C_R=1.3,
                    C_D=2.3,
                    hp_exponent=6.0,
                    include_sun=(omitted != "sun"),
                    include_moon=(omitted != "moon"),
                    include_srp=(omitted != "srp"),
                    include_drag=include_drag,
                )

            def derivative_function(t, state):
                return orbit_state_derivative(
                    t,
                    state,
                    acceleration_function,
                )

            return integrate_rk4(
                derivative_function,
                state0,
                duration,
                steps,
            )

        print(f"\n{name}: reference model", flush=True)
        print(
            f"Duration: {duration:.3f} s | "
            f"Step: {duration / n_steps:.3f} s",
            flush=True,
        )

        times, reference = propagate()

        print(f"{name}: step-halving check...", flush=True)
        times_fine, reference_fine = propagate(steps=2 * n_steps)

        np.testing.assert_allclose(
            times,
            times_fine[::2],
            rtol=0.0,
            atol=1e-8,
        )

        convergence_difference_m = 1000.0 * np.linalg.norm(
            reference[:, :3] - reference_fine[::2, :3],
            axis=1,
        )

        maximum_difference = float(
            np.max(convergence_difference_m)
        )

        print(
            "Maximum sampled step-halving difference: "
            f"{maximum_difference:.6f} m",
            flush=True,
        )

        if maximum_difference > CONVERGENCE_LIMIT_M:
            raise RuntimeError(
                f"{name}: step-halving difference exceeds "
                f"{CONVERGENCE_LIMIT_M} m. "
                "Reduce this satellite's step and rerun."
            )

        del reference_fine, times_fine

        satellite_results = {}

        for label, n_max, m_max, omitted in models:
            # GEO drag is disabled in the reference as well.
            if omitted == "drag" and not satellite["drag"]:
                satellite_results[label] = [0.0, 0.0]
                continue

            print(f"{name}: {label}...", flush=True)

            model_times, states = propagate(
                n_max=n_max,
                m_max=m_max,
                omitted=omitted,
            )

            np.testing.assert_allclose(
                times,
                model_times,
                rtol=0.0,
                atol=1e-8,
            )

            delta_r = states[:, :3] - reference[:, :3]
            difference_m = 1000.0 * np.linalg.norm(
                delta_r,
                axis=1,
            )

            maxima = []

            for end_time in intervals:
                # Maximum over all grid samples in this interval.
                mask = times <= end_time
                maximum = float(np.max(difference_m[mask]))

                # Include the exact interval endpoint approximately
                # when it falls between integration grid samples.
                if times[mask][-1] < end_time:
                    delta_at_end = np.array([
                        np.interp(end_time, times, delta_r[:, axis])
                        for axis in range(3)
                    ])
                    endpoint_difference = (
                        1000.0 * np.linalg.norm(delta_at_end)
                    )
                    maximum = max(maximum, endpoint_difference)

                maxima.append(maximum)

            satellite_results[label] = maxima

        results[name] = satellite_results

    print("\nExercise 3.4 — Orbit Perturbations")
    print("Maximum sampled position differences [m]")
    print("Reference: gravity 20/20 + Sun + Moon + SRP")
    print("Harris-Priester Drag included for LEO only.")
    print()
    print(
        f"{'Model':<12}"
        f"{'LEO 1 rev':>14}"
        f"{'LEO 1 day':>14}"
        f"{'GEO 1 day':>14}"
        f"{'GEO 2 days':>14}"
    )
    print("-" * 68)

    for label, _, _, _ in models:
        values = results["LEO"][label] + results["GEO"][label]
        print(
            f"{label:<12}"
            + "".join(f"{value:14.3f}" for value in values)
        )

    print("\nGEO No Drag = 0 by the model assumption.")
    print(
        "LEO one-revolution endpoint uses linear interpolation "
        "of the position difference."
    )


if __name__ == "__main__":
    main()
