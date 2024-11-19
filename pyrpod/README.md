# Core Modules

The `pyrpod/` directory contains core Python modules for simulations, mission planning, and modeling of vehicle dynamics. Each module performs a specific function in the simulation pipeline.

## Modules
- `IsentropicExpansion.py`: Calculates isentropic expansion properties.
- `JetFiringHistory.py`: Manages jet firing data and history.
- `MissionPlanner.py`: Core script for setting up and running missions.
- `TargetVehicle.py` and `VisitingVehicle.py`: Models for the vehicles in the simulation.

## Usage
Use test cases to run simulations:
```bash
cd ../tests 
python3 <module_test>/test_case.py
```

Example: Running all test cases.
```bash
cd ../tests
python3 test_all.py
```

## Dependencies
Install dependencies from `requirements.txt`. Additional dependencies may be needed for specific modules; see inline comments within each module.

## Adding New Modules
Place new Python scripts here with functions well-documented using docstrings. New functionality should include properly defined test cases. Ensure compatibility with existing modules by following established naming conventions and testing procedures.