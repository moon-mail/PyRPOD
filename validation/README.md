# Validation

The `validation/` directory contains scripts to validate simulation results against expected outcomes, ensuring the accuracy and integrity of the model.

## Files
Each file represents a validation case (e.g., `validation_case_01.py`), which compares simulation outputs against benchmark data.

## Running Validation
To run a validation case:
```bash
python validation/<validation_case>.py
```

Example:
```bash
python validation/validation_case_01.py
```

## Adding New Validation Scripts
1. Create a new script with a descriptive name (`validation_case_<number>.py`).
2. Include detailed comments on the validation process.
3. Use `assert` statements or comparison functions to check the output.