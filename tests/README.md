# Tests

This directory contains unit, integration, and verification tests. The tests ensure that the modules and cases function correctly and produce accurate results.

## Structure
- `mdao/`: Tests for the MDAO (Multidisciplinary Design, Analysis, and Optimization) components.
- `mission/`: Tests for mission planning modules.
- `plume/`: Tests for plume impingement dynamics.
- `rpod/`: Tests related to the rendezvous and proximity operations docking (RPOD).
- `validation/`: Validation tests to cross-check physical models against known data.

## Running Tests
Use `pytest` to run all tests:
```bash
pytest
```

To run tests in a specific module:
```bash
pytest tests/<module_name>
```

## Adding New Tests
1. Place new tests in the relevant subdirectory.
2. Follow naming conventions (`test_<feature>.py`) and include descriptive docstrings.
3. Use `assert` statements to validate output.