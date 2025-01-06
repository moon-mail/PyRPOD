# Data Directory

This folder contains shared data files, flight plans, STL models, and thruster configuration data used across cases.

## Subdirectories
- `flight_plan/`: Contains CSV files with different flight plans (to be deleted?).
- `jfh/`: Jet Firing History files in `.A` format.
- `stl/`: Contains STL models for objects like cylinders, plates, and the ISS.
- `tcd/`: Thruster Configuration Data files.

## Usage
These files are referenced in cases and modules in `pyrpod/` and `cases/`. You can add new data here if needed for custom simulations.

## Adding New Data
To add new data:
1. Place the file in the relevant subfolder.
2. Update any relevant case or module configurations to reference the new file.