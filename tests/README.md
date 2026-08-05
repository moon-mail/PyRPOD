<!-- AUTO-GENERATED FILE. DO NOT EDIT DIRECTLY.
     Edit tests/test_manifest.yaml and rerun the inventory generator:
         python scripts/generate_test_dashboard.py -->

# PyRPOD Test Inventory

This inventory documents every test asset in `tests/`: what it is for,
which subsystem and category it belongs to, how it is executed, and how
far its development has progressed.

**It deliberately records no pass/fail results.** Execution outcomes belong
to a specific run and are published in the pytest HTML report
(`reports/pyrpod-pytest-report.html`); development status is long-lived
metadata that a green test run does not change.

Test files follow the `<subsystem>_<category>_test_NN.py` naming convention,
and `tests/conftest.py` uses it to tag every collected test with a subsystem
marker (`logging`, `mdao`, `mission`, `plume`, `rpod`, `tooling`) and a
category marker (`unit`, `integration`, `verification`).

Regenerate this file and the HTML report with:

```bash
python scripts/generate_test_dashboard.py
```

Source of truth for the metadata below: [`test_manifest.yaml`](test_manifest.yaml).

---

## Summary

| Metric | Count |
| --- | --- |
| Manifest entries | 102 |
| Collected by pytest (files) | 53 |
| Collected by pytest (test cases) | 440 |
| Manual verification scripts | 42 |
| Placeholder tests | 12 |
| Blocked tests | 5 |
| Archived / legacy tests | 5 |
| Needs review | 10 |

---

## Automated pytest tests

Files pytest collects and runs. The **Cases** column is the number of test
cases pytest currently collects from the file (parametrized variants count
separately).

### Logging

#### Unit

| Test file | Cases | Description | Development status | Reference |
| --- | --- | --- | --- | --- |
| [`logging_unit_test_01.py`](logging/logging_unit_test_01.py) | 14 | Unit tests for the centralized operational logging system (pyrpod.logging_utils): import side effects, handler ownership, configuration precedence, console/file toggles, runtime-log naming and location, configuration snapshots, input-asset logging and array summaries. | `implemented` | — |

#### Integration

| Test file | Cases | Description | Development status | Reference |
| --- | --- | --- | --- | --- |
| [`logging_integration_test_01.py`](logging/logging_integration_test_01.py) | 13 | Integration tests for operational logging of the plume-strike workflow: fail-fast validation of required inputs, progress cadence, serial/parallel event equivalence, DEBUG-vs-INFO artifact levels, parallel-to-serial fallback with successful_with_warning status and optional-visualization warn-and-continue behavior. | `implemented` | — |

### MDAO

#### Unit

| Test file | Cases | Description | Development status | Reference |
| --- | --- | --- | --- | --- |
| [`mdao_unit_test_01.py`](mdao/mdao_unit_test_01.py) | 1 | Placeholder MDAO unit test; the test body returns immediately and asserts nothing. | `placeholder` | — |
| [`mdao_unit_test_02.py`](mdao/mdao_unit_test_02.py) | 1 | Intended to build an array of cant-angle-swept thruster configurations (symmetric pitch/yaw canting) and visualize each sweep step. The entire body is currently commented out, so the test asserts nothing. | `placeholder` | — |
| [`mdao_unit_test_03.py`](mdao/mdao_unit_test_03.py) | 23 | Unit tests for prescribed firing generation and the exact meaning of n_firings: invalid counts are rejected, N requested firings produce exactly N JFH entries (read back through JetFiringHistory), an explicit firing list that disagrees with n_firings is an error, a whole-sweep sequence holds exactly poses x n_firings pose-tagged entries, the generated pose convention reproduces the committed flat-plate sweep JFH, and the dynamics-driven approach honors an exact count. | `implemented` | — |
| [`mdao_unit_test_04.py`](mdao/mdao_unit_test_04.py) | 20 | Unit tests for per-component surface-load integration on meshes with analytically known loading: pressure and shear force integration, moments about a user-defined reference point, center of pressure (recovered, moment-consistent, and unavailable for zero-load and near-cancellation cases), coefficient calculation and omission, and component face selection. | `implemented` | — |
| [`mdao_unit_test_05.py`](mdao/mdao_unit_test_05.py) | 34 | Unit tests for the YAML study configuration: the committed flat-plate examples parse, paths resolve against the configuration file, validation rejects an unsupported plume model, an unknown sweep mode, invalid firing counts, mismatched firing lists and bad geometry or normalization inputs, the per-pose n_firings semantics hold in both sweep modes, and the case's own config.ini keeps parsing unchanged. | `implemented` | — |
| [`mdao_unit_test_06.py`](mdao/mdao_unit_test_06.py) | 22 | Unit tests for the generic external-reference comparison (absolute and relative error, normalized RMSE, peak error, integrated-load error, center-of-pressure displacement; CSV / JSON / YAML datasets compared identically regardless of origin; unmatched cases and missing quantities reported rather than fabricated) and for the structured result schema's CSV and JSON serialization. | `implemented` | — |
| [`mdao_unit_test_07.py`](mdao/mdao_unit_test_07.py) | 63 | Unit tests for the ISS-panel study extensions: collisionless plume-model dispatch (both Cai variants selectable by name, unknown names rejected, both reduced to one common LocalFieldState feeding one shared Maxwellian gas-surface implementation); panel-local pose generation (parallel_to_normal translates the source without tilting the plume axis, and coincides with aim_at_reference at zero offset); offset-sweep enumeration and the refusal to combine incompatible pose definitions; derived Knudsen metadata for both reference-length modes; panel-local coordinate transforms and hand-computed normal-force and moment projections; and the distribution-CSV schema. | `implemented` | Cai & Wang 2012, J. Spacecraft and Rockets 49(2), 335-340, doi:10.2514/1.A32046 |

#### Integration

| Test file | Cases | Description | Development status | Reference |
| --- | --- | --- | --- | --- |
| [`mdao_integration_test_01.py`](mdao/mdao_integration_test_01.py) | 1 | Placeholder MDAO integration test; the test body returns immediately and asserts nothing. | `placeholder` | — |
| [`mdao_integration_test_02.py`](mdao/mdao_integration_test_02.py) | 8 | Baseline flat-plate case run end to end through the package-level TradeStudy.from_config API: one case with exactly one JFH entry, the standard per-face strike VTK at the advertised path, CSV and JSON summaries with full provenance, physically consistent head-on loads, available coefficients, and agreement of the integrated normal load with the independent Cai 2016 exact reference. | `implemented` | Cai 2016, Aerospace 3(4), 43, doi:10.3390/aerospace3040043, Eq. 15 |
| [`mdao_integration_test_03.py`](mdao/mdao_integration_test_03.py) | 13 | Multi-angle / multi-distance flat-plate sweep through the trade-study API (a reduced subset of the committed sweep configuration): one result per angle-distance case, mirror symmetry in +/- angle, load decay with distance, center-of-pressure travel, agreement with the independent Cai 2016 reference through the generic comparison interface, per-case VTK isolation, optional trend plots, and the same machinery running on the cylinder target with coefficients correctly unavailable. | `implemented` | Cai 2016, Aerospace 3(4), 43, doi:10.3390/aerospace3040043, Eq. 15 |
| [`mdao_integration_test_04.py`](mdao/mdao_integration_test_04.py) | 12 | The single-history sweep engine (sweep.mode: single_jfh): engine dispatch from the configuration, exactly one Jet Firing History holding poses x n_firings entries, one result record per firing keyed to its pose, per-firing equivalence with the per-case engine, a single results/strikes VTK series, the per-component sweep envelope (absent rather than fabricated without the full pipeline), and the same independent Cai 2016 reference comparison. | `implemented` | Cai 2016, Aerospace 3(4), 43, doi:10.3390/aerospace3040043, Eq. 15 |
| [`mdao_integration_test_05.py`](mdao/mdao_integration_test_05.py) | 41 | The committed ISS-representative solar-panel case (22 m x 12 m flat panel, source translated parallel to the panel) run end to end through TradeStudy.from_config().run(): panel-local resultants, derived Knudsen metadata and distribution export on the committed baseline; proof that selecting CollisionlessGasKinetics changes the ANSWER and not merely the metadata; the pressure force re-integrated from the exported per-face values matching the CaseResult; symmetry of a centered source and the documented moment signs for offset sources; offset-sweep case count, order and identifiers; headless plot generation; and backward compatibility of the existing flat-plate configurations, case identifiers and baseline results. | `implemented` | Cai & Wang 2012, J. Spacecraft and Rockets 49(2), 335-340, doi:10.2514/1.A32046 |

#### Verification

| Test file | Cases | Description | Development status | Reference |
| --- | --- | --- | --- | --- |
| [`mdao_verification_test_04.py`](mdao/mdao_verification_test_04.py) | 1 | Intended variable-sweep study of the maximum overshoot velocity a logistics module can absorb with a fixed thruster configuration and deceleration start distance. The entire body is commented out, so the test asserts nothing. | `placeholder` | — |
| [`mdao_verification_test_05.py`](mdao/mdao_verification_test_05.py) | 1 | Intended trade study sweeping the surface cant angle of the RCS pack. The entire body is commented out, so the test asserts nothing. | `placeholder` | — |
| [`mdao_verification_test_06.py`](mdao/mdao_verification_test_06.py) | 1 | Intended multi-variable trade study sweeping axial overshoot together with surface cant angle. The entire body is commented out, so the test asserts nothing. | `placeholder` | — |

### Mission

#### Unit

| Test file | Cases | Description | Development status | Reference |
| --- | --- | --- | --- | --- |
| [`mission_unit_test_01.py`](mission/mission_unit_test_01.py) | 1 | Placeholder mission unit test; the test body returns immediately and asserts nothing. | `placeholder` | — |

#### Integration

| Test file | Cases | Description | Development status | Reference |
| --- | --- | --- | --- | --- |
| [`mission_integration_test_01.py`](mission/mission_integration_test_01.py) | 1 | Calculates the 6DOF performance (normal vector, force, translational acceleration, torque) of each individual thruster in the logistics module. Expected-value helpers are present but no assertion is made against them, so the test only exercises the code path. | `needs_review` | — |
| [`mission_integration_test_02.py`](mission/mission_integration_test_02.py) | 1 | Intended to calculate the 6DOF performance of thruster working groups. The body only assigns mass-distribution constants; the header records that the analysis was disabled because it created unwanted data files, so the test asserts nothing. | `placeholder` | — |
| [`mission_integration_test_03.py`](mission/mission_integration_test_03.py) | 1 | Calculates RCS performance for a given flight plan, approximating the 1D delta-v requirements. Exercises the flight-evaluation pipeline end to end but makes no assertion. | `needs_review` | — |
| [`mission_integration_test_04.py`](mission/mission_integration_test_04.py) | 1 | Analyzes a notional 1D translation-plus-rotation approach (module header marks it "NEEDS TLC"). Exercises the flight-evaluation pipeline end to end but makes no assertion. | `needs_review` | — |
| [`mission_integration_test_05.py`](mission/mission_integration_test_05.py) | 1 | Graphs thrust versus time or distance for the given design requirements to establish an RCS thrust-requirement flight envelope. Exercises the plotting pipeline but makes no assertion. | `needs_review` | — |
| [`mission_integration_test_06.py`](mission/mission_integration_test_06.py) | 1 | Graphs propellant-mass requirements across the flight envelope for the given design requirements. Exercises the plotting pipeline but makes no assertion. | `needs_review` | — |
| [`mission_integration_test_07.py`](mission/mission_integration_test_07.py) | 1 | Contours the burn-time plot across a range of thrust and Isp values (module header marks it "NEEDS TLC"). Exercises the plotting pipeline but makes no assertion. | `needs_review` | — |
| [`mission_integration_test_08.py`](mission/mission_integration_test_08.py) | 1 | Contours propellant usage across the various delta-v legs of a given flight plan. Exercises the plotting pipeline but makes no assertion. | `needs_review` | — |

#### Verification

| Test file | Cases | Description | Development status | Reference |
| --- | --- | --- | --- | --- |
| [`mission_verification_test_01.py`](mission/mission_verification_test_01.py) | 1 | Placeholder mission verification test; the test body returns immediately and asserts nothing. | `placeholder` | — |
| [`mission_verification_test_02.py`](mission/mission_verification_test_02.py) | 1 | Builds and summarizes a chain of Hohmann transfers (LEO to MEO to GEO) through MissionPlanner.orbital_transfer. Exercises the orbital transfer path but makes no assertion. | `needs_review` | — |

### Plume

#### Unit

| Test file | Cases | Description | Development status | Reference |
| --- | --- | --- | --- | --- |
| [`plume_unit_test_01.py`](plume/plume_unit_test_01.py) | 1 | Placeholder plume unit test; the test body returns immediately and asserts nothing. | `placeholder` | — |
| [`plume_unit_test_02.py`](plume/plume_unit_test_02.py) | 9 | Verifies the closed form of the special factor Q against a direct 50-term truncation of the printed Legendre series, its centerline reduction, and the 0 < Q <= 1 bound that keeps the exp(-S0^2 (1-Q)) combination overflow-safe. | `implemented` | Cai & Wang 2012, JSR 49(1), DOI 10.2514/1.A32046, Eq. 9 |
| [`plume_unit_test_03.py`](plume/plume_unit_test_03.py) | 13 | Tests the Simons cosine-law plume model after the gamma generalization: isentropic throat ratios from gamma, zero density beyond the limiting angle, parameterizable beaming exponent kappa and exit-referenced density scaling. | `implemented` | Cai & Wang 2012, JSR 49(1), DOI 10.2514/1.A32046, Sec. II.C, Eqs. 25-26 |
| [`plume_unit_test_04.py`](plume/plume_unit_test_04.py) | 3 | Asserts that the NumPy-vectorized plume-strike detection in compute_plume_strikes() reproduces the scalar reference implementation exactly for every firing of case/rpod/1d_approach (geometry only) and case/rpod/multi_thrusters_square (Simplified kinetics, multiple thrusters). | `implemented` | — |

#### Integration

| Test file | Cases | Description | Development status | Reference |
| --- | --- | --- | --- | --- |
| [`plume_integration_test_01.py`](plume/plume_integration_test_01.py) | 1 | Placeholder plume integration test; the test body returns immediately and asserts nothing. | `placeholder` | — |
| [`plume_integration_test_02.py`](plume/plume_integration_test_02.py) | 3 | Asserts that jfh_plume_strikes() keeps its default serial, return-compatible behavior, produces identical per-firing and cumulative output when the optional process-parallel path is enabled, and rejects invalid worker counts with a clear error. | `implemented` | — |

#### Verification

| Test file | Cases | Description | Development status | Reference |
| --- | --- | --- | --- | --- |
| [`plume_verification_test_01.py`](plume/plume_verification_test_01.py) | 1 | Intended to plot simple isentropic radial-expansion profiles. Both plotting calls are commented out, leaving only construction of an IsentropicExpansion object, so the test asserts nothing. | `placeholder` | — |
| [`plume_verification_test_02.py`](plume/plume_verification_test_02.py) | 23 | Pinning tests for SimplifiedGasKinetics against reference values computed independently from the paper's equations with mpmath at 40 significant digits, plus far-field asymptote convergence and the near-field divergence of the corrected Eq. 21 centerline temperature quadrature. | `implemented` | Cai & Wang 2012, JSR 49(1), DOI 10.2514/1.A32046, Eqs. 13-19, 22-24 |
| [`plume_verification_test_03.py`](plume/plume_verification_test_03.py) | 27 | Verification of the full collisionless model CollisionlessGasKinetics against the paper's internal anchors: centerline reduction to the Eq. 18/19 closed forms, vanishing W, the Eq. 30 far-field error bound, monotonic centerline density decay and asymptote convergence. | `implemented` | Cai & Wang 2012, JSR 49(1), DOI 10.2514/1.A32046, Eqs. 5-12, 20, 30 |

### RPOD

#### Unit

| Test file | Cases | Description | Development status | Reference |
| --- | --- | --- | --- | --- |
| [`rpod_unit_test_01.py`](rpod/rpod_unit_test_01.py) | 1 | Converts cylinder STL data to VTK and asserts the resulting file has the proper VTK data format, cell counts and point counts. | `implemented` | — |
| [`rpod_unit_test_02.py`](rpod/rpod_unit_test_02.py) | 1 | Reads a jet firing history and asserts that all required keys exist and that every value has the expected data type. | `implemented` | — |
| [`rpod_unit_test_03.py`](rpod/rpod_unit_test_03.py) | 3 | Captures the exact file output of the three JFH printing helpers in pyrpod.util.io.file_print to a temporary directory and compares it against the committed tests/rpod/jfh_outputs snapshot. | `implemented` | — |
| [`rpod_unit_test_04.py`](rpod/rpod_unit_test_04.py) | 31 | Unit tests for the PR #116 plume-mesh consolidation: VisitingVehicle.transform_plume_mesh, the cached thruster-id map, multi-digit cluster ids, stl.compose_meshes and the reusable stl.transform_mesh mesh-object API, pinned against hand-computed golden values. | `implemented` | — |

#### Integration

| Test file | Cases | Description | Development status | Reference |
| --- | --- | --- | --- | --- |
| [`rpod_integration_test_01.py`](rpod/rpod_integration_test_01.py) | 1 | Asserts the expected number of cell strikes on a flat-plate STL for a notional "sweep" trajectory above it across 20 distinct firings. This is the established base case for RPOD plume-impingement analysis. | `implemented` | — |
| [`rpod_integration_test_02.py`](rpod/rpod_integration_test_02.py) | 1 | Asserts the expected number of cell strikes on a flat-plate STL as a notional visiting vehicle approaches it along a 1D-physics trajectory, across 15 distinct firings. | `implemented` | — |
| [`rpod_integration_test_03.py`](rpod/rpod_integration_test_03.py) | 1 | Analyzes keep-out-zone impingement and asserts the expected strike counts against the committed expected-strikes fixture (module header marks the case "WIP"). | `implemented` | — |
| [`rpod_integration_test_04.py`](rpod/rpod_integration_test_04.py) | 1 | Produces hollow-cube target data and asserts the expected strike counts against the committed expected-strikes fixture. | `implemented` | — |
| [`rpod_integration_test_05.py`](rpod/rpod_integration_test_05.py) | 1 | Runs the plume gas-kinetic models through a JFH with multiple thrusters per firing (case/rpod/multi_thrusters_square). Exercises graph_jfh and jfh_plume_strikes end to end but makes no assertion. | `needs_review` | — |
| [`rpod_integration_test_06.py`](rpod/rpod_integration_test_06.py) | 3 | Performance benchmark for plume-strike computation (scalar vs NumPy-vectorized geometry, serial vs process-parallel). Elapsed times are printed but no speedup threshold is asserted; each benchmark asserts only that the compared paths produce identical strike arrays. | `implemented` | — |
| [`rpod_integration_test_07.py`](rpod/rpod_integration_test_07.py) | 1 | Runs the Cai 2016 inclined-plate case through PlumeStrikeEstimationStudy, converts per-face dimensional loads to the paper's coefficient normalization and compares them face by face against the exact reference functions evaluated at face centroids. | `implemented` | Cai 2016, Aerospace 3(4), 43, doi:10.3390/aerospace3040043, Eqs. 9-13 |
| [`rpod_integration_test_08.py`](rpod/rpod_integration_test_08.py) | 4 | End-to-end geometry-equivalence regressions for the PR #116 plume-mesh consolidation: pins the visualization pipeline's produced geometry, artifact counts and file names against an inline reproduction of the legacy transform/compose sequence, and proves the thruster transform is applied exactly once. | `implemented` | — |

#### Verification

| Test file | Cases | Description | Development status | Reference |
| --- | --- | --- | --- | --- |
| [`rpod_verification_test_04.py`](rpod/rpod_verification_test_04.py) | 1 | Visualizes STLs after decoupling the thruster configuration data and is intended to verify that plume strikes line up with the thrusters. Exercises graph_jfh and jfh_plume_strikes but makes no assertion, so the alignment must still be checked by eye in ParaView. | `needs_review` | — |
| [`rpod_verification_test_06.py`](rpod/rpod_verification_test_06.py) | 1 | Cai 2016 flat-plate SWEEP verification: runs the 95-firing sweep JFH (19 approach angles x 5 stand-off distances) through the strike pipeline, reduces each pose to the Eq.-15 plate-averaged coefficients, checks mirror symmetry in +/- alpha and compares against the exact reference envelope. | `implemented` | Cai 2016, Aerospace 3(4), 43, doi:10.3390/aerospace3040043, Eq. 15 |
| [`rpod_verification_test_07.py`](rpod/rpod_verification_test_07.py) | 1 | Cylinder-target sweep smoke test: exercises the full PlumeStrikeEstimationStudy pipeline on a curved, closed target (14036-face cylinder) over 19 approach angles x 5 orbit radii. | `implemented` | — |

### Tooling

#### Unit

| Test file | Cases | Description | Development status | Reference |
| --- | --- | --- | --- | --- |
| [`tooling_unit_test_01.py`](tooling/tooling_unit_test_01.py) | 27 | Unit tests for the test-inventory tooling itself: manifest schema validation, detection of missing and stale entries, and deterministic rendering of the generated tests/README.md. | `implemented` | — |

---

## Manual verification scripts

Run directly from the repository root; they are **not** pytest tests and
never appear in the HTML execution report. Most reproduce a published
figure and write a PNG under `tests/plume/output/`.

### MDAO

| Script | Description | Development status | Command | Reference |
| --- | --- | --- | --- | --- |
| [`mdao_verification_test_01.py`](mdao/mdao_verification_test_01.py) | OpenMDAO axial-positioning optimizer: minimizes the maximum heat-flux load for a 1D approach by varying the axial position of the thruster packs. | `blocked` | `python tests/mdao/mdao_verification_test_01.py` | — |
| [`mdao_verification_test_02.py`](mdao/mdao_verification_test_02.py) | OpenMDAO cant-angle optimizer: minimizes the maximum heat-flux load for a 1D approach by varying the cant angle of the deceleration thrusters. | `blocked` | `python tests/mdao/mdao_verification_test_02.py` | — |
| [`mdao_verification_test_03.py`](mdao/mdao_verification_test_03.py) | OpenMDAO axial-plus-cant evaluator: reports maximum cumulative heat-flux load and JFH propellant expenditure for a 1D approach given the axial position and cant angle of the deceleration thrusters. | `blocked` | `python tests/mdao/mdao_verification_test_03.py` | — |

### Plume

| Script | Description | Development status | Command | Reference |
| --- | --- | --- | --- | --- |
| [`plume_impingement_error_summary.py`](plume/plume_impingement_error_summary.py) | Maximum and mean relative differences between the exact Cai 2016 reference surface coefficients (Eqs. 9-14) and the current PyRPOD approximation chain for the inclined-plate impingement study. | `implemented` | `python tests/plume/plume_impingement_error_summary.py` | Cai 2016, Aerospace 3(4), 43, doi:10.3390/aerospace3040043, Figs. 17-21 conditions |
| [`plume_verification_error_summary.py`](plume/plume_verification_error_summary.py) | Model-vs-model analog of Table 1: maximum relative differences in density and mass flux between the analytical, simplified and Simons models along the centerline and along r/D = 10 at S0 = 2.0. The DSMC reference columns are emitted as "pending digitized data". | `implemented` | `python tests/plume/plume_verification_error_summary.py` | Cai & Wang 2012, JSR 49(1), DOI 10.2514/1.A32046, Table 1 |
| [`plume_verification_test_04.py`](plume/plume_verification_test_04.py) | Reproduces Fig. 2: plume boundaries (the n/n0 = 0.001 contour of the full analytical model) for exit speed ratios S0 = 1, 2, 3. | `implemented` | `python tests/plume/plume_verification_test_04.py` | Cai & Wang 2012, JSR 49(1), DOI 10.2514/1.A32046, Fig. 2 |
| [`plume_verification_test_05.py`](plume/plume_verification_test_05.py) | Reproduces Fig. 3: normalized analytical centerline number density (Eq. 18) versus X/D for S0 = 1, 2, 3. | `implemented` | `python tests/plume/plume_verification_test_05.py` | Cai & Wang 2012, JSR 49(1), DOI 10.2514/1.A32046, Fig. 3 |
| [`plume_verification_test_06.py`](plume/plume_verification_test_06.py) | Reproduces Fig. 4: normalized analytical centerline U-velocity (Eq. 19) versus X/D for S0 = 1, 2, 3. | `implemented` | `python tests/plume/plume_verification_test_06.py` | Cai & Wang 2012, JSR 49(1), DOI 10.2514/1.A32046, Fig. 4 |
| [`plume_verification_test_07.py`](plume/plume_verification_test_07.py) | Reproduces Fig. 5: normalized analytical centerline temperature (Eq. 21, exact quadrature) versus X/D for S0 = 1, 2, 3. | `implemented` | `python tests/plume/plume_verification_test_07.py` | Cai & Wang 2012, JSR 49(1), DOI 10.2514/1.A32046, Fig. 5 |
| [`plume_verification_test_08.py`](plume/plume_verification_test_08.py) | Reproduces Fig. 6: normalized number-density contours at Kn = 100, S0 = 2.0; analytical plus simplified in the upper half-plane, Simons plus the DSMC overlay slot in the lower half-plane. | `implemented` | `python tests/plume/plume_verification_test_08.py` | Cai & Wang 2012, JSR 49(1), DOI 10.2514/1.A32046, Fig. 6 |
| [`plume_verification_test_09.py`](plume/plume_verification_test_09.py) | Reproduces Fig. 7: normalized number-density contours at Kn = 0.1, S0 = 2.0. Same analytic content as Fig. 6; the Kn distinction lives in the DSMC overlay slot. | `implemented` | `python tests/plume/plume_verification_test_09.py` | Cai & Wang 2012, JSR 49(1), DOI 10.2514/1.A32046, Fig. 7 |
| [`plume_verification_test_10.py`](plume/plume_verification_test_10.py) | Reproduces Fig. 8: normalized number-density contours at Kn = 0.01, S0 = 2.0. Same analytic content as Fig. 6; the Kn distinction lives in the DSMC overlay slot. | `implemented` | `python tests/plume/plume_verification_test_10.py` | Cai & Wang 2012, JSR 49(1), DOI 10.2514/1.A32046, Fig. 8 |
| [`plume_verification_test_11.py`](plume/plume_verification_test_11.py) | Reproduces Fig. 9: relative density error between the analytical and Simons cosine-law solutions, \|n_A/n_Simons - 1\| (Eq. 28), S0 = 2.0. | `implemented` | `python tests/plume/plume_verification_test_11.py` | Cai & Wang 2012, JSR 49(1), DOI 10.2514/1.A32046, Fig. 9 |
| [`plume_verification_test_12.py`](plume/plume_verification_test_12.py) | Reproduces Fig. 10: relative density error between the analytical and simplified analytical solutions, \|n_A/n_As - 1\| (Eq. 29), S0 = 2.0. | `implemented` | `python tests/plume/plume_verification_test_12.py` | Cai & Wang 2012, JSR 49(1), DOI 10.2514/1.A32046, Fig. 10 |
| [`plume_verification_test_13.py`](plume/plume_verification_test_13.py) | Reproduces Fig. 11: normalized pressure contours p1/p0 = (n1/n0)(T1/T0) at Kn = 100, S0 = 2.0. | `implemented` | `python tests/plume/plume_verification_test_13.py` | Cai & Wang 2012, JSR 49(1), DOI 10.2514/1.A32046, Fig. 11 |
| [`plume_verification_test_14.py`](plume/plume_verification_test_14.py) | Reproduces Fig. 12: normalized pressure contours p1/p0 at Kn = 0.1, S0 = 2.0. Same analytic content as Fig. 11. | `implemented` | `python tests/plume/plume_verification_test_14.py` | Cai & Wang 2012, JSR 49(1), DOI 10.2514/1.A32046, Fig. 12 |
| [`plume_verification_test_15.py`](plume/plume_verification_test_15.py) | Reproduces Fig. 13: normalized temperature contours T1/T0 at Kn = 100, S0 = 2.0 - the figure the paper uses to argue that p = n k T0 is invalid because T1 < T0 everywhere downstream. | `implemented` | `python tests/plume/plume_verification_test_15.py` | Cai & Wang 2012, JSR 49(1), DOI 10.2514/1.A32046, Fig. 13 |
| [`plume_verification_test_16.py`](plume/plume_verification_test_16.py) | Reproduces Fig. 14: normalized U-velocity contours U1*sqrt(beta0) at Kn = 100, S0 = 2.0. | `implemented` | `python tests/plume/plume_verification_test_16.py` | Cai & Wang 2012, JSR 49(1), DOI 10.2514/1.A32046, Fig. 14 |
| [`plume_verification_test_17.py`](plume/plume_verification_test_17.py) | Reproduces Fig. 15: normalized transverse-velocity contours at Kn = 100, S0 = 2.0. In the plotted XOZ plane the y-component is zero by axisymmetry, so the model's W (Eq. 7) is plotted. | `implemented` | `python tests/plume/plume_verification_test_17.py` | Cai & Wang 2012, JSR 49(1), DOI 10.2514/1.A32046, Fig. 15 |
| [`plume_verification_test_18.py`](plume/plume_verification_test_18.py) | Reproduces Fig. 16: normalized radial-velocity Vr contours at Kn = 100, S0 = 2.0. | `implemented` | `python tests/plume/plume_verification_test_18.py` | Cai & Wang 2012, JSR 49(1), DOI 10.2514/1.A32046, Fig. 16 |
| [`plume_verification_test_19.py`](plume/plume_verification_test_19.py) | Reproduces Fig. 17: normalized radial-velocity Vr contours at Kn = 0.1, S0 = 2.0. Same analytic content as Fig. 16. | `implemented` | `python tests/plume/plume_verification_test_19.py` | Cai & Wang 2012, JSR 49(1), DOI 10.2514/1.A32046, Fig. 17 |
| [`plume_verification_test_20.py`](plume/plume_verification_test_20.py) | Reproduces Fig. 18: normalized radial-velocity Vr contours at Kn = 0.01, S0 = 2.0. Same analytic content as Fig. 16. | `implemented` | `python tests/plume/plume_verification_test_20.py` | Cai & Wang 2012, JSR 49(1), DOI 10.2514/1.A32046, Fig. 18 |
| [`plume_verification_test_21.py`](plume/plume_verification_test_21.py) | Reproduces Fig. 19: centerline density profiles at Kn = 100, S0 = 2.0 - analytical (Eq. 18), simplified (Eq. 14), Simons (exit-referenced cosine law) and the DSMC overlay slot. | `implemented` | `python tests/plume/plume_verification_test_21.py` | Cai & Wang 2012, JSR 49(1), DOI 10.2514/1.A32046, Fig. 19 |
| [`plume_verification_test_22.py`](plume/plume_verification_test_22.py) | Reproduces Fig. 20: centerline density profiles at Kn = 0.1, S0 = 2.0. Same analytic content as Fig. 19. | `implemented` | `python tests/plume/plume_verification_test_22.py` | Cai & Wang 2012, JSR 49(1), DOI 10.2514/1.A32046, Fig. 20 |
| [`plume_verification_test_23.py`](plume/plume_verification_test_23.py) | Reproduces Fig. 21: centerline density profiles at Kn = 0.01, S0 = 2.0. Same analytic content as Fig. 19. | `implemented` | `python tests/plume/plume_verification_test_23.py` | Cai & Wang 2012, JSR 49(1), DOI 10.2514/1.A32046, Fig. 21 |
| [`plume_verification_test_24.py`](plume/plume_verification_test_24.py) | Reproduces Fig. 22: density profiles along r/D = 10 at Kn = 100, S0 = 2.0 - analytical, simplified, Simons with kappa = 1.5/2/3 under a single Boyton normalization, and the DSMC overlay slot. | `implemented` | `python tests/plume/plume_verification_test_24.py` | Cai & Wang 2012, JSR 49(1), DOI 10.2514/1.A32046, Fig. 22 |
| [`plume_verification_test_25.py`](plume/plume_verification_test_25.py) | Reproduces Fig. 23: density profiles along r/D = 10 at Kn = 0.1, S0 = 2.0. Same analytic content as Fig. 22. | `implemented` | `python tests/plume/plume_verification_test_25.py` | Cai & Wang 2012, JSR 49(1), DOI 10.2514/1.A32046, Fig. 23 |
| [`plume_verification_test_26.py`](plume/plume_verification_test_26.py) | Reproduces Fig. 24: density profiles along r/D = 10 at Kn = 0.01, S0 = 2.0. Same analytic content as Fig. 22. | `implemented` | `python tests/plume/plume_verification_test_26.py` | Cai & Wang 2012, JSR 49(1), DOI 10.2514/1.A32046, Fig. 24 |
| [`plume_verification_test_27.py`](plume/plume_verification_test_27.py) | Reproduces Fig. 25: normalized mass flux along r/D = 10 versus theta at S0 = 2.0, together with three DSMC overlay slots. The module header records the flux-normalization convention adopted to match the printed magnitudes. | `implemented` | `python tests/plume/plume_verification_test_27.py` | Cai & Wang 2012, JSR 49(1), DOI 10.2514/1.A32046, Fig. 25 |
| [`plume_verification_test_28.py`](plume/plume_verification_test_28.py) | Reproduces Fig. 17: diffuse-plate surface pressure contours Cp,d(s, tau) for the round argon jet on the inclined plate (Eq. 9), overlaid with the current PyRPOD approximation. | `implemented` | `python tests/plume/plume_verification_test_28.py` | Cai 2016, Aerospace 3(4), 43, doi:10.3390/aerospace3040043, Fig. 17 |
| [`plume_verification_test_29.py`](plume/plume_verification_test_29.py) | Reproduces Fig. 18: specular-plate surface pressure contours Cp,s(s, tau) (Eq. 14), overlaid with the current PyRPOD approximation at sigma = 0. | `implemented` | `python tests/plume/plume_verification_test_29.py` | Cai 2016, Aerospace 3(4), 43, doi:10.3390/aerospace3040043, Fig. 18 |
| [`plume_verification_test_30.py`](plume/plume_verification_test_30.py) | Reproduces Fig. 19: diffuse-plate friction coefficient Cf1,d(s, tau) along the inclined direction (Eq. 11), whose zero contour marks the stagnation line, overlaid with the current PyRPOD approximation. | `implemented` | `python tests/plume/plume_verification_test_30.py` | Cai 2016, Aerospace 3(4), 43, doi:10.3390/aerospace3040043, Fig. 19 |
| [`plume_verification_test_31.py`](plume/plume_verification_test_31.py) | Reproduces Fig. 20: diffuse-plate friction coefficient Cf2,d(s, tau) along the horizontal direction (Eq. 12), antisymmetric in s, overlaid with the current PyRPOD approximation. | `implemented` | `python tests/plume/plume_verification_test_31.py` | Cai 2016, Aerospace 3(4), 43, doi:10.3390/aerospace3040043, Fig. 20 |
| [`plume_verification_test_32.py`](plume/plume_verification_test_32.py) | Reproduces Fig. 21: diffuse-plate heat-flux coefficient Cq,d(s, tau) (Eq. 13), peaking just beneath the plate center, overlaid with the current PyRPOD approximation. | `implemented` | `python tests/plume/plume_verification_test_32.py` | Cai 2016, Aerospace 3(4), 43, doi:10.3390/aerospace3040043, Fig. 21 |
| [`plume_verification_test_33.py`](plume/plume_verification_test_33.py) | Reproduces Fig. 5: temperature contours T/T0 for the 2D slot jet impinging on an inclined DIFFUSE planar plate, the combined free-jet plus wall-emission field of Section 3. | `implemented` | `python tests/plume/plume_verification_test_33.py` | Cai 2016, Aerospace 3(4), 43, doi:10.3390/aerospace3040043, Fig. 5 |
| [`plume_verification_test_34.py`](plume/plume_verification_test_34.py) | Reproduces Fig. 6: temperature contours T/T0 for the 2D slot jet impinging on an inclined SPECULAR planar plate, using the paper's virtual-nozzle construction (Eq. 7). | `implemented` | `python tests/plume/plume_verification_test_34.py` | Cai 2016, Aerospace 3(4), 43, doi:10.3390/aerospace3040043, Fig. 6 |
| [`plume_verification_test_35.py`](plume/plume_verification_test_35.py) | Reproduces Fig. 7: 2D diffuse-plate surface pressure profiles Cp,d(s) (Eq. 2) for four (S0, alpha0) combinations at Tw/T0 = 1.5. | `implemented` | `python tests/plume/plume_verification_test_35.py` | Cai 2016, Aerospace 3(4), 43, doi:10.3390/aerospace3040043, Fig. 7 |
| [`plume_verification_test_36.py`](plume/plume_verification_test_36.py) | Reproduces Fig. 8: 2D specular-plate surface pressure profiles Cp,s(s) (Eq. 8) for four (S0, alpha0) combinations. | `implemented` | `python tests/plume/plume_verification_test_36.py` | Cai 2016, Aerospace 3(4), 43, doi:10.3390/aerospace3040043, Fig. 8 |
| [`plume_verification_test_37.py`](plume/plume_verification_test_37.py) | Reproduces Fig. 9: 2D diffuse-plate surface friction profiles Cf,d(s) (Eq. 3), including the near-s/(2H) = -2 zero crossing the paper identifies as a possible flow-separation spot. | `implemented` | `python tests/plume/plume_verification_test_37.py` | Cai 2016, Aerospace 3(4), 43, doi:10.3390/aerospace3040043, Fig. 9 |
| [`plume_verification_test_38.py`](plume/plume_verification_test_38.py) | Reproduces Fig. 10: 2D diffuse-plate surface heat-flux profiles Cq,d(s) (Eq. 4) for four (S0, alpha0) combinations at Tw/T0 = 1.5. | `implemented` | `python tests/plume/plume_verification_test_38.py` | Cai 2016, Aerospace 3(4), 43, doi:10.3390/aerospace3040043, Fig. 10 |
| [`plume_verification_test_39.py`](plume/plume_verification_test_39.py) | Reproduces Fig. 15: static-pressure contours p/p0 in the Y = 0 plane for the round jet impinging on an inclined DIFFUSE rectangular plate. | `implemented` | `python tests/plume/plume_verification_test_39.py` | Cai 2016, Aerospace 3(4), 43, doi:10.3390/aerospace3040043, Fig. 15 |
| [`plume_verification_test_40.py`](plume/plume_verification_test_40.py) | Reproduces Fig. 16: static-pressure contours p/p0 in the Y = 0 plane for the round jet impinging on an inclined SPECULAR rectangular plate, using the paper's 3D virtual-nozzle construction. | `implemented` | `python tests/plume/plume_verification_test_40.py` | Cai 2016, Aerospace 3(4), 43, doi:10.3390/aerospace3040043, Fig. 16 |

---

## Placeholder tests

Collected by pytest but containing no assertion or verification behavior -
typically an empty body or a fully commented-out one. Each is explicitly
skipped in its source file so that the HTML report shows it as **skipped**
rather than passed; a placeholder must never be mistaken for coverage.

| Test file | Subsystem | Category | Description |
| --- | --- | --- | --- |
| [`mdao_unit_test_01.py`](mdao/mdao_unit_test_01.py) | mdao | unit | Placeholder MDAO unit test; the test body returns immediately and asserts nothing. |
| [`mdao_unit_test_02.py`](mdao/mdao_unit_test_02.py) | mdao | unit | Intended to build an array of cant-angle-swept thruster configurations (symmetric pitch/yaw canting) and visualize each sweep step. The entire body is currently commented out, so the test asserts nothing. |
| [`mdao_integration_test_01.py`](mdao/mdao_integration_test_01.py) | mdao | integration | Placeholder MDAO integration test; the test body returns immediately and asserts nothing. |
| [`mdao_verification_test_04.py`](mdao/mdao_verification_test_04.py) | mdao | verification | Intended variable-sweep study of the maximum overshoot velocity a logistics module can absorb with a fixed thruster configuration and deceleration start distance. The entire body is commented out, so the test asserts nothing. |
| [`mdao_verification_test_05.py`](mdao/mdao_verification_test_05.py) | mdao | verification | Intended trade study sweeping the surface cant angle of the RCS pack. The entire body is commented out, so the test asserts nothing. |
| [`mdao_verification_test_06.py`](mdao/mdao_verification_test_06.py) | mdao | verification | Intended multi-variable trade study sweeping axial overshoot together with surface cant angle. The entire body is commented out, so the test asserts nothing. |
| [`mission_unit_test_01.py`](mission/mission_unit_test_01.py) | mission | unit | Placeholder mission unit test; the test body returns immediately and asserts nothing. |
| [`mission_integration_test_02.py`](mission/mission_integration_test_02.py) | mission | integration | Intended to calculate the 6DOF performance of thruster working groups. The body only assigns mass-distribution constants; the header records that the analysis was disabled because it created unwanted data files, so the test asserts nothing. |
| [`mission_verification_test_01.py`](mission/mission_verification_test_01.py) | mission | verification | Placeholder mission verification test; the test body returns immediately and asserts nothing. |
| [`plume_unit_test_01.py`](plume/plume_unit_test_01.py) | plume | unit | Placeholder plume unit test; the test body returns immediately and asserts nothing. |
| [`plume_integration_test_01.py`](plume/plume_integration_test_01.py) | plume | integration | Placeholder plume integration test; the test body returns immediately and asserts nothing. |
| [`plume_verification_test_01.py`](plume/plume_verification_test_01.py) | plume | verification | Intended to plot simple isentropic radial-expansion profiles. Both plotting calls are commented out, leaving only construction of an IsentropicExpansion object, so the test asserts nothing. |

---

## Ignored or blocked tests

Excluded from pytest collection (see `collect_ignore` in
[`conftest.py`](conftest.py)) or otherwise unable to run.

| Test file | Subsystem | Development status | Collection status | Reason | Command |
| --- | --- | --- | --- | --- | --- |
| [`test_case_25.py`](test_case_25.py) | rpod | `blocked` | `ignored` | Listed in tests/conftest.py collect_ignore: imports the flat `pyrpod.<Module>` layout that no longer exists after the rpod -> plume/vehicle refactor. Excluded from collection until it is fixed or removed outright. | — |
| [`rpod_verification_test_05.py`](rpod_verification_test_05.py) | rpod | `blocked` | `ignored` | Listed in tests/conftest.py collect_ignore: imports the flat `pyrpod.<Module>` layout that no longer exists after the rpod -> plume/vehicle refactor. Excluded from collection until it is fixed or removed outright. | — |
| [`mdao_verification_test_01.py`](mdao/mdao_verification_test_01.py) | mdao | `blocked` | `manual` | Defines an openmdao ExplicitComponent rather than a pytest test, so pytest collects nothing from it. The module header records "This has been having issues, disregard for now (4/2/24)". | `python tests/mdao/mdao_verification_test_01.py` |
| [`mdao_verification_test_02.py`](mdao/mdao_verification_test_02.py) | mdao | `blocked` | `manual` | Defines an openmdao ExplicitComponent rather than a pytest test, so pytest collects nothing from it. The module header records "This has been having issues, disregard for now (4/3/24)", and the __main__ block contains an incomplete call (prob.set_val with a missing value). | `python tests/mdao/mdao_verification_test_02.py` |
| [`mdao_verification_test_03.py`](mdao/mdao_verification_test_03.py) | mdao | `blocked` | `manual` | Defines an openmdao ExplicitComponent rather than a pytest test, so pytest collects nothing from it. The module header records "This has been having issues, disregard for now (4/3/24)", and the __main__ block contains an incomplete call (prob.set_val with a missing value). | `python tests/mdao/mdao_verification_test_03.py` |

---

## Archived or legacy tests

Kept for historical reference under `tests/old/`. They are not repaired,
renamed, or deleted as part of routine work.

| Test file | Subsystem | Category | Description | Reason |
| --- | --- | --- | --- | --- |
| [`old/test_case_sweep_cants.py`](old/test_case_sweep_cants.py) | mdao | unit | Legacy sweep building an array of cant-angle-swept thruster configurations with symmetric pitch/yaw canting. Superseded by tests/mdao/mdao_unit_test_02.py. | tests/old is listed in tests/conftest.py collect_ignore. Imports the flat `pyrpod.<Module>` layout removed by the rpod -> plume/vehicle refactor, plus a `test_header` module that no longer exists. |
| [`old/test_case_sweep_coords.py`](old/test_case_sweep_coords.py) | mdao | unit | Legacy sweep building an array of axially swept thruster configurations on a common ring x-coordinate. | tests/old is listed in tests/conftest.py collect_ignore. Imports the flat `pyrpod.<Module>` layout removed by the rpod -> plume/vehicle refactor, plus a `test_header` module that no longer exists. |
| [`old/test_case_17.py`](old/test_case_17.py) | rpod | unit | Legacy test of STL-to-VTK conversion, checking the VTK data format. Superseded by tests/rpod/rpod_unit_test_01.py. | tests/old is listed in tests/conftest.py collect_ignore. Imports the flat `pyrpod.<Module>` layout removed by the rpod -> plume/vehicle refactor, plus a `test_header` module that no longer exists. |
| [`old/test_case_15.py`](old/test_case_15.py) | rpod | integration | Legacy test of the plume gas-kinetic models in JFH firings. | tests/old is listed in tests/conftest.py collect_ignore. Imports the flat `pyrpod.<Module>` layout removed by the rpod -> plume/vehicle refactor, plus a `test_header` module that no longer exists. |
| [`old/test_case_19.py`](old/test_case_19.py) | rpod | integration | Legacy test producing hollow-cube data. Superseded by tests/rpod/rpod_integration_test_04.py. | tests/old is listed in tests/conftest.py collect_ignore. Imports the flat `pyrpod.<Module>` layout removed by the rpod -> plume/vehicle refactor, plus a `test_header` module that no longer exists. |

---

## Legend

Three independent axes. A test can be `implemented` and still fail today;
a `placeholder` can be green in CI only because it is skipped.

**Execution outcome** — owned by pytest, one value per run, published only
in `reports/pyrpod-pytest-report.html`:

| Outcome | Meaning |
| --- | --- |
| passed | The test ran and every assertion held. |
| failed | The test ran and an assertion or the code under test failed. |
| error | The test could not run to completion (setup or teardown raised). |
| skipped | The test was not executed (placeholder, or a skip condition). |

**Development status** — maintained in `test_manifest.yaml`, long-lived:

| Status | Meaning |
| --- | --- |
| `implemented` | Complete: exercises the code and asserts a result. |
| `placeholder` | No assertion or verification behavior yet; skipped. |
| `needs_review` | Runs real code but asserts nothing, or its purpose is insufficiently documented. |
| `blocked` | Cannot run against the current architecture. |
| `archived` | Superseded; kept for historical reference only. |
| `deprecated` | Slated for removal. |

**Collection status** — whether pytest picks the file up:

| Status | Meaning |
| --- | --- |
| `collected` | Collected and run by pytest; appears in the HTML report. |
| `manual` | Run by hand; defines no pytest tests. |
| `ignored` | Listed in `collect_ignore` in `conftest.py`. |
| `archived` | Under `tests/old/`, excluded from collection wholesale. |
