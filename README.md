# Astrodynamics: Orbit Propagation, Determination & Trajectory Analysis

This repository documents an ongoing study and implementation project in astrodynamics based on two primary references:

- Oliver Montenbruck & Eberhard Gill — *Satellite Orbits: Models, Methods and Applications*
- Howard D. Curtis — *Orbital Mechanics for Engineering Students*

The project combines theoretical study with Python-based computational exercises and reusable implementations in orbital mechanics, flight dynamics, orbit propagation, and orbit determination.

The current computational work is primarily based on Montenbruck & Gill. Methods and applications from Curtis will be progressively integrated as the study expands toward broader trajectory and interplanetary analysis.

## Project Objective

The project progresses from fundamental orbital mechanics and physical modeling toward increasingly realistic numerical propagation, measurement modeling, orbit determination, and trajectory analysis.

Reusable methods developed through the exercises are collected in a common Python package, allowing later problems to extend and reuse previously implemented functionality.

The long-term objective is to develop the accumulated work into a compact astrodynamics analysis toolkit covering both Earth-orbit flight dynamics and interplanetary trajectory problems.

## Current Scope

The repository currently includes implementations related to:

- Keplerian orbital mechanics
- Orbital state and element transformations
- Orbital maneuvers and basic orbit design
- Coordinate and reference-frame transformations
- Time-system calculations
- Preliminary orbit determination
- Earth gravity and perturbation models
- Sun and Moon ephemerides
- Atmospheric drag and solar radiation pressure
- Numerical orbit propagation
- Satellite tracking and measurement modeling
- Light-time solutions and media corrections

## Repository Structure

```text
astrodynamics-orbit-analysis/
│
├── exercises/
│   ├── chapter_02/
│   ├── chapter_03/
│   ├── chapter_04/
│   ├── chapter_05/
│   └── chapter_06/
│
├── src/
│   └── satellite_orbits/
│       ├── astrodynamics/
│       │   ├── keplerian/
│       │   ├── orbit_design/
│       │   ├── orbit_determination/
│       │   ├── perturbations/
│       │   └── propagation/
│       ├── coordinates/
│       ├── ephemerides/
│       ├── observations/
│       ├── time_systems/
│       └── constants.py
│
├── README.md
├── LICENSE
├── .gitignore
└── pyproject.toml
```

### `exercises/`

Chapter-based computational exercises and numerical studies.

Each exercise contains its own implementation and documentation while reusing functionality from the common source package where appropriate. Results are verified against the values published in the textbook, and any discrepancy is investigated and documented in the exercise README rather than tuned away.

### `src/satellite_orbits/`

Reusable functions and models developed throughout the project.

- **`astrodynamics/keplerian/`** — Kepler's equation, conversions between state vectors and orbital elements, mean anomaly propagation.
- **`astrodynamics/orbit_design/`** — sun-synchronous conditions, draconic period relations, orbital maneuvers.
- **`astrodynamics/orbit_determination/`** — preliminary orbit determination from two position vectors.
- **`astrodynamics/perturbations/`** — Earth gravity field, third-body attraction, atmospheric drag, solar radiation pressure, secular J2 rates.
- **`astrodynamics/propagation/`** — analytical Keplerian propagation and numerical integrators, together with the equation of motion they advance.
- **`coordinates/`** — inertial, Earth-fixed, geodetic and topocentric frames.
- **`ephemerides/`** — Sun and Moon positions, Chebyshev evaluation, interpolation of tabulated ephemerides.
- **`observations/`** — ground-station states, light-time solutions, range and range rate models, ionospheric and tropospheric corrections.
- **`time_systems/`** — Modified Julian Date and sidereal time.

## Installation

Python 3.11 or newer is recommended.

From the repository root:

```bash
pip install -e .
```

This installs the local `satellite_orbits` package in editable mode so that the exercise scripts can use the reusable modules under `src/`.

## Conventions

- Lengths are in kilometres and times in seconds throughout the package, with conversions confined to the exercise scripts where a textbook table uses other units.
- Angles are in radians in all function signatures. Degrees appear only in printed output.
- Functions take explicit physical parameters rather than reading global constants, so that a model can be exercised with alternative values.
- Comments and documentation are in English.

## Status

Work in progress.

Chapters 2 to 6 of Montenbruck & Gill are implemented. The repository will continue to evolve alongside the study of Montenbruck & Gill and Curtis, with additional implementations and integrated flight-dynamics and interplanetary trajectory applications added progressively.

## License

The source code in this repository is released under the MIT License.

The referenced textbooks and their contents remain the intellectual property of their respective authors and publishers.
