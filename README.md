# Astrodynamics: Orbit Propagation, Determination & Trajectory Analysis

This repository documents an ongoing study and implementation project in astrodynamics based on two primary references:

- Oliver Montenbruck & Eberhard Gill — *Satellite Orbits: Models, Methods and Applications*
- Howard D. Curtis — *Orbital Mechanics for Engineering Students*

The project combines theoretical study with Python-based computational exercises and reusable implementations in orbital mechanics, flight dynamics, orbit propagation, and orbit determination.

The current computational work is primarily based on Montenbruck & Gill. Methods and applications from Curtis will be progressively integrated as the study expands toward broader trajectory and interplanetary analysis.

## Project Objective

The project progresses from fundamental orbital mechanics and physical modeling toward increasingly realistic numerical propagation, orbit determination, and trajectory analysis.

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

## Repository Structure

```text
astrodynamics-orbit-analysis/
│
├── exercises/
│   ├── chapter_02/
│   └── chapter_03/
│
├── src/
│   └── satellite_orbits/
│       ├── astrodynamics/
│       ├── coordinates/
│       ├── ephemerides/
│       └── time_systems/
│
├── README.md
├── LICENSE
├── .gitignore
└── pyproject.toml
```

### `exercises/`

Chapter-based computational exercises and numerical studies.

Each exercise contains its own implementation and documentation while reusing functionality from the common source package where appropriate.

### `src/satellite_orbits/`

Reusable functions and models developed throughout the project.

The package currently contains functionality for orbital mechanics, orbit design, orbit determination, perturbation modeling, numerical propagation, coordinate transformations, ephemerides, and time systems.

## Installation

Python 3.11 or newer is recommended.

From the repository root:

```bash
pip install -e .
```

This installs the local `satellite_orbits` package in editable mode so that the exercise scripts can use the reusable modules under `src/`.

## Status

Work in progress.

The repository will continue to evolve alongside the study of Montenbruck & Gill and Curtis, with additional implementations and integrated flight-dynamics and interplanetary trajectory applications added progressively.

## License

The source code in this repository is released under the MIT License.

The referenced textbooks and their contents remain the intellectual property of their respective authors and publishers.