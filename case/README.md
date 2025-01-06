# Case Studies

This directory contains simulation cases, each representing a specific configuration. Each case has a configuration file (`config.ini`), STL files for geometry, and transformation scripts.

## Structure
Each subfolder represents a unique case study:
- `1d_approach/`: 1-dimensional approach configurations.
- `base_case/`: Baseline scenario for testing.
- `cant_sweep/`: Case with varying cant angles for thrusters.
- `multi_thrusters_square/`: Simulation with multiple thrusters in a square arrangement.

## Usage
To run a case:
1. Execute the corresponding test case (for now).
2. Example: running 1D Approach is hosed in file_name.py:
   ```bash
   python pyrpod/MissionPlanner.py --config case/<case_name>/config.ini
   ```

Results will be saved in each case's `results/` subdirectory.

## Adding a New Case
1. Create a new folder under `case/`.
2. Add your `config.ini`, `stl/` files, and transformations.
3. Run the case using the above commands.