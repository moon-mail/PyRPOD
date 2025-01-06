# PyRPOD Test Dashboard

This dashboard provides an overview of all tests in the PyRPOD project, categorized by their respective modules. Each test is listed with its type, description, and current status.

---

## **MDAO Tests**
| Test Name                       | Type           | Description                                   | Status |
|---------------------------------|----------------|-----------------------------------------------|--------|
| `mdao_integration_test_01.py`   | Integration    | Tests the integration of MDAO subsystems.    | ✅     |
| `mdao_unit_test_01.py`          | Unit           | Verifies individual MDAO functionalities.    | ✅     |
| `mdao_unit_test_02.py`          | Unit           | Ensures MDAO parameter calculations work.    | ✅     |
| `mdao_verification_test_02.py`  | Verification   | Confirms outputs align with requirements.    | ❌     |
| `mdao_verification_test_03.py`  | Verification   | Validates design constraints in MDAO.        | ❌     |
| `mdao_verification_test_04.py`  | Verification   | Ensures optimization results are correct.    | ✅     |

---

## **Mission Tests**
| Test Name                       | Type           | Description                                   | Status |
|---------------------------------|----------------|-----------------------------------------------|--------|
| `mission_integration_test_01.py`| Integration    | Tests mission logic under nominal conditions.| ✅     |
| `mission_integration_test_02.py`| Integration    | Verifies inter-system mission data flow.     | ❌     |
| `mission_integration_test_03.py`| Integration    | Tests fault-tolerant mission execution.      | ✅     |
| `mission_integration_test_04.py`| Integration    | Simulates edge-case mission scenarios.       | ✅     |
| `mission_integration_test_05.py`| Integration    | Evaluates performance during long missions.  | ❌     |
| `mission_integration_test_06.py`| Integration    | Validates automated mission recovery.        | ✅     |
| `mission_integration_test_07.py`| Integration    | Tests mission re-planning capabilities.      | ❌     |
| `mission_integration_test_08.py`| Integration    | Assesses mission data collection.            | ✅     |
| `mission_integration_test_09.py`| Integration    | Tests payload integration in missions.       | ❌     |
| `mission_unit_test_01.py`       | Unit           | Verifies basic mission operations.           | ✅     |
| `mission_verification_test_01.py`| Verification  | Confirms compliance with mission objectives. | ✅     |

---

## **Plume Tests**
| Test Name                       | Type           | Description                                   | Status |
|---------------------------------|----------------|-----------------------------------------------|--------|
| `plume_integration_test_01.py`  | Integration    | Tests plume modeling in integrated systems.  | ✅     |
| `plume_unit_test_01.py`         | Unit           | Verifies individual plume calculation methods.| ✅     |
| `plume_verification_test_01.py` | Verification   | Validates plume outputs against benchmarks.  | ❌     |

---

## **RPOD Tests**
| Test Name                       | Type           | Description                                   | Status |
|---------------------------------|----------------|-----------------------------------------------|--------|
| `rpod_integration_test_01.py`   | Integration    | Simulates nominal rendezvous scenarios.      | ✅     |
| `rpod_integration_test_02.py`   | Integration    | Tests proximity operations under stress.     | ❌     |
| `rpod_integration_test_03.py`   | Integration    | Verifies docking sequence and safety.        | ✅     |
| `rpod_integration_test_04.py`   | Integration    | Evaluates advanced rendezvous maneuvers.     | ❌     |
| `rpod_unit_test_01.py`          | Unit           | Verifies core RPOD calculations.             | ✅     |
| `rpod_unit_test_02.py`          | Unit           | Tests position tracking during docking.      | ✅     |
| `rpod_unit_test_03.py`          | Unit           | Assesses collision-avoidance algorithms.     | ❌     |
| `rpod_verification_test_03.py`  | Verification   | Confirms compliance with docking requirements.| ✅     |
| `rpod_verification_test_04.py`  | Verification   | Validates safety protocols during docking.   | ❌     |

---

## **Legacy/Old Tests**
| Test Name                       | Type           | Description                                   | Status |
|---------------------------------|----------------|-----------------------------------------------|--------|
| `test_case_15.py`               | Miscellaneous  | Old test for deprecated functionality.       | ❌     |
| `test_case_17.py`               | Miscellaneous  | Tests legacy feature interactions.           | ✅     |
| `test_case_19.py`               | Miscellaneous  | Validates compatibility of old methods.      | ❌     |
| `test_case_sweep_cants.py`      | Sweep Test     | Evaluates various canting angles.            | ✅     |
| `test_case_sweep_coords.py`     | Sweep Test     | Tests coordinate system transformations.     | ❌     |

---

## **Status Legend**
- ✅ = Passed
- ❌ = Failed/Requires Attention

This dashboard serves as a quick reference for test organization and tracking within PyRPOD. Update the statuses regularly to ensure it reflects the latest testing outcomes.
