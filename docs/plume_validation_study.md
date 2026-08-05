# Prescribed plume-validation trade studies

This document covers the repaired `TradeStudy` architecture and the
package-level API for **prescribed** plume/target validation sweeps: studies
in which the firing poses are placed by the engineer (or generated from a
swept approach angle and source distance about a stationary target) rather
than flown from vehicle dynamics.

- [Quick start](#quick-start)
- [Architecture](#architecture)
- [Plume models](#plume-models)
- [Sweep modes: one JFH per case, or one for the sweep](#sweep-modes-one-jfh-per-case-or-one-for-the-sweep)
- [Source axis modes: aiming vs translating](#source-axis-modes-aiming-vs-translating)
- [Panel-local u/v offsets](#panel-local-uv-offsets)
- [Derived Knudsen metadata](#derived-knudsen-metadata)
- [YAML configuration](#yaml-configuration)
- [`n_firings` and prescribed firings](#n_firings-and-prescribed-firings)
- [Integrated loads](#integrated-loads)
  - [Force](#force)
  - [Moment](#moment)
  - [Center of pressure](#center-of-pressure)
  - [Thermal quantities](#thermal-quantities)
  - [Coefficients](#coefficients)
- [Result schema](#result-schema)
- [VTK outputs](#vtk-outputs)
- [Panel-local surface distributions](#panel-local-surface-distributions)
- [Optional plots](#optional-plots)
- [External reference data](#external-reference-data)
- [Adding another validation geometry](#adding-another-validation-geometry)
- [ISS-panel example](#iss-panel-example)
- [Known limitations](#known-limitations)

---

## Quick start

A complete flat-plate example ships with the repository. It runs the Cai 2016
Section-4 conditions (argon round jet, `D = 1 m`, `S0 = 2.0`, `T0 = 200 K`,
`Tw = 300 K`, fully diffuse) against the 8 m x 8 m plate of
`case/plume/plume_flat_plate_sweep`, with the plume source head-on at
`L = 4D`:

```python
from pyrpod.mdao.TradeStudy import TradeStudy

study = TradeStudy.from_config(
    'case/plume/plume_flat_plate_sweep/study/flat_plate_baseline.yaml')
results = study.run()

case = results.cases[0]
print(case.force)                 # [~0, ~0, -2.763]  N  (plume pushes -Z)
print(case.moment)                # ~0 about the plate center
print(case.center_of_pressure)    # ~[0, 0, 0]
print(case.max_pressure)          # 0.340 Pa
print(case.max_heat_flux)         # 50.1 W/m^2
print(case.coefficients['CF'])    # 0.0391
print(case.vtk_path)              # .../cases/case000_.../results/strikes/firing-0.vtu
print(results.summary_csv_path)   # .../case_results.csv
```

The multi-angle / multi-distance sweep (19 angles x 5 distances) is the same
call against `flat_plate_sweep.yaml`:

```python
study = TradeStudy.from_config(
    'case/plume/plume_flat_plate_sweep/study/flat_plate_sweep.yaml')
results = study.run()
study.plot()                       # optional trend figures
```

The **ISS-representative panel** example (22 m x 12 m panel, source
translated parallel to it, either plume model, panel-local distributions and
derived Knudsen metadata) is the same call — see
[ISS-panel example](#iss-panel-example):

```python
study = TradeStudy.from_config(
    'case/plume/iss_panel_thesis/study/iss_panel_offset_distance_sweep.yaml')
results = study.run()

case = results.cases[0]
print(case.normal_force)          # + into the panel
print(case.local_moment_v)        # panel moment about the transverse axis
print(case.center_of_pressure_u)  # panel-local, from the panel center
print(case.knudsen_number)        # derived metadata; not a model input
print(case.surface_distribution_path)
```

Outputs land under the configured `output_dir` (by default
`<case>/results/studies/<study name>/`, which is gitignored):

```
results/studies/flat_plate_sweep/
  jfh/case000_alpha0p0_d4.A            one prescribed JFH per case
  cases/case000_alpha0p0_d4/results/strikes/firing-0.vtu
  sweep_results.csv                    one row per case x component x firing
  sweep_metadata.json                  provenance + nested case records
  plots/                               optional trend figures
  reference_comparison.csv             written when reference data is configured
```

To run a committed configuration into a scratch location (this is what the
automated tests do), pass `output_dir`:

```python
TradeStudy.from_config(path, output_dir='/tmp/my_study').run()
```

---

## Architecture

`TradeStudy` is a thin façade. The work lives in small modules with one
responsibility each, all under `pyrpod/mdao/`:

| Module | Responsibility |
| --- | --- |
| `study_config.py` | Parse and validate the YAML study configuration |
| `firing_plan.py` | Build prescribed firings; write Jet Firing Histories |
| `plume_validation.py` | **Engine**: one JFH per case (`mode: per_case`) |
| `parameter_sweep.py` | **Engine**: one JFH for the sweep (`mode: single_jfh`) |
| `study_runtime.py` | Plumbing both engines share (assets, geometry, strikes, records) |
| `surface_loads.py` | Integrate per-face fields into component loads |
| `study_results.py` | Result schema; CSV + JSON output |
| `surface_distribution.py` | Panel-local per-face distribution CSV + sidecar |
| `reference_data.py` | Generic external-reference comparison and metrics |
| `study_plots.py` | Optional angle-sweep trend figures |
| `panel_plots.py` | Optional offset-sweep and panel-pressure figures |

Plume-model dispatch lives under `pyrpod/plume/` in
`gas_kinetics_models.py`, next to the models themselves.

Everything domain-specific is delegated to the existing PyRPOD objects —
`JetFiringHistory`, `TargetVehicle`, `VisitingVehicle`, `MissionEnvironment`
and `PlumeStrikeEstimationStudy` — so a study inherits the pipeline's
validation, logging, plume physics and VTK conventions instead of
reimplementing them.

`TradeStudy.from_config` picks the engine from `sweep.mode` (see the next
section); both do the same seven things, differing only in how many Jet
Firing Histories the sweep is decomposed into:

1. build the prescribed firings;
2. write the Jet Firing History and read it back with the normal
   `JetFiringHistory` parser;
3. run the plume-strike calculation;
4. integrate per-component loads;
5. export the per-face VTK fields;
6. record a structured result;
7. optionally compare against external reference data.

### Repairs to the legacy sweeps

The dynamics-driven sweeps (`run_axial_overshoot_sweep`,
`run_surface_cant_sweep`, `run_multi_var_sweep`) were broken before any
physics ran. Fixed on this branch:

- `init_trade_study` constructed `PlumeStrikeEstimationStudy.RPOD(case_dir)`
  — a class that does not exist — and passed a case directory where the study
  takes a `MissionEnvironment`. It now builds the environment and the real
  `PlumeStrikeEstimationStudy`.
- The sweeps called `jfh_plume_strikes(trade_study=True)`, which that
  method's `(parallel, workers)` signature rejects with `TypeError`.
- `print_mission_report` read impingement maxima from study attributes that
  are never set; it now summarizes them from the per-firing arrays
  `jfh_plume_strikes` returns (the old attribute path still works for
  existing callers).
- `interpret_mission_report` read `self.rpod.config`, which does not exist;
  the configuration lives on the study's `MissionEnvironment`.

### Execution paths

With VTK output enabled (the default) each case runs the full
`PlumeStrikeEstimationStudy.jfh_plume_strikes()` pipeline, with the target's
output root temporarily redirected to that case's folder so sweep cases
cannot overwrite one another's artifacts. With VTK output disabled the study
calls the pipeline's own per-firing core, `compute_plume_strikes()` — the
same function `jfh_plume_strikes` calls per firing, so the numbers are
identical — which keeps large sweeps and automated tests fast and
artifact-free. This applies to both engines.

---

## Plume models

Two **collisionless** analytical plume models are selectable, both from
Cai & Wang 2012 and both already verified in this repository:

| `plume_model.name` | `[pm] kinetics` | Model |
| --- | --- | --- |
| `SimplifiedGasKinetics` | `Simplified` | Far-field simplification (the `Q'` of Eq. 13) — a closed form, fast |
| `CollisionlessGasKinetics` | `Collisionless` | Full model: the exact factor `Q` (Eq. 9) integrated over the finite exit disk, valid in the near field |

`SimplifiedGasKinetics` is the default, so **a configuration that names no
model behaves exactly as before**. Any other name raises `StudyConfigError`.

**The name selects the model that computes the plume field — it is not
metadata.** `study_runtime.load_case_assets` applies it by setting the
environment's *in-memory* `[pm] kinetics` key, which is the one input both
strike paths read, so the selection reaches the calculation itself. The
case's `config.ini` on disk is never modified, and a study naming the model
its case already configures changes nothing. A case with
`[pm] kinetics = None` is rejected: a validation study needs surface loads.

Dispatch is a small registry in `pyrpod/plume/gas_kinetics_models.py`, not a
plugin framework — `CollisionlessGasKinetics` subclasses
`SimplifiedGasKinetics` with an identical constructor, so selecting between
them is a class lookup. Both are reduced to one common `LocalFieldState`:

```
number_density, mass_density, axial_velocity, radial_velocity,
velocity_magnitude, temperature, speed_ratio
```

and a single `maxwellian_surface_loads()` applies the Shen gas-surface
interaction to it, so **pressure, shear and heat-transfer logic exists once**
rather than once per model. On the plume centerline both models use the same
exact closed forms and therefore agree bit-for-bit; off-axis they differ,
which is what makes the selection meaningful.

`pyrpod/plume/CaiImpingement2016.py` stays **independent of the strike
pipeline** and is not a model backend. It is used only as an external
analytical reference (see [External reference data](#external-reference-data)
and `case/plume/iss_panel_thesis/cai2016_reference.py`).

---

## Sweep modes: one JFH per case, or one for the sweep

The same configuration can be decomposed two ways, chosen with `sweep.mode`.
Per-firing numbers are **identical** either way — an automated test asserts
it pose by pose — so the choice is about what you want to look at, not about
accuracy.

| | `per_case` (default) | `single_jfh` |
| --- | --- | --- |
| Engine | `PlumeValidationStudy` | `ParameterSweepStudy` |
| Jet Firing Histories | one per angle-distance case | **one**, spanning every pose |
| `case_id` / `firing_id` | one case per pose, firings 1..`n_firings` | one case, firings 1..`poses x n_firings` |
| Strike runs | one per case | one, over the whole history |
| VTK layout | `cases/<case_id>/results/strikes/firing-<i>.vtu` | `results/strikes/firing-<i>.vtu` — one series, scrubbable in ParaView |
| Cumulative fields | restart each case | accumulate across the sweep |
| Sweep envelope | n/a | recorded per component |

Use `per_case` when each pose is an independent experiment (a validation
matrix compared pose-by-pose against reference data): nothing one pose does
can affect another's artifacts. Use `single_jfh` when the poses form one
sequence and you want the pipeline's cumulative fields to mean something:

- `max_pressures` / `max_shears` — the worst load each face saw **anywhere**
  in the sweep;
- `cum_strikes` — a coverage map over the whole sweep;
- `cum_heat_flux_load` — accumulated heat-flux dose.

Those per-face fields live in the last firing's VTK file; their per-component
summary is exposed as `study.validation_study.envelope` and written into the
metadata document under `provenance.sweep_envelope`:

```python
study = TradeStudy.from_config('.../flat_plate_sweep_single_jfh.yaml')
results = study.run()
study.validation_study.envelope['plate']
# {'max_pressure': 0.340, 'max_shear_stress': 0.120,
#  'max_heat_flux_load': 126.95, 'total_strike_events': 31104.0,
#  'unique_struck_faces': 10368, 'swept_affected_area': 64.0,
#  'component_area': 64.0}
```

The envelope needs the pipeline's cumulative arrays, which only the full
strike path produces: with `output.vtk.enabled: false` it is reported as
empty rather than reconstructed from the per-firing records.

`n_firings` keeps one meaning in both modes — **entries contributed by each
pose** — so a configuration switches modes by changing one line. In
`single_jfh` mode the single history is required to hold exactly
`len(poses) x n_firings` entries, checked after writing and again after
reading the file back.

---

## Source axis modes: aiming vs translating

`sweep.source_axis_mode` chooses how the plume axis is oriented at each
generated pose. These are **different experiments**, not two spellings of one.

For a target reference point `C` with normal `n` (toward the source),
longitudinal axis `u` (= `target.tangent`) and transverse axis `v = n x u`:

### `aim_at_reference` (default — existing behavior)

```
d(alpha)  = cos(alpha)*n + sin(alpha)*u
position  = C + L*d(alpha)
axis      = -d(alpha)                     # re-aimed at C at every angle
```

The source rides an arc of radius `L` about `C` and **always points at
`C`**. `plate_angles_deg` sweeps the approach angle; `alpha = 0` is head-on.
This is what every existing study uses, and it is unchanged.

### `parallel_to_normal` (ISS-panel studies)

```
position  = C + L*n + u_offset*u + v_offset*v
axis      = -n                            # fixed; never re-aimed
```

The source is **translated parallel to the surface** with its axis held
fixed, so the plume centerline meets the panel at `(u_offset, v_offset)`
instead of always at `C`. `plate_angles_deg` has no meaning here and is
rejected if set.

**At zero offset the two modes coincide exactly** (`parallel_to_normal`
equals `aim_at_reference` at `alpha = 0`, position and DCM), so the new mode
extends the old convention rather than redefining it. In both modes the JFH
DCM carries the thruster axis in its **first column**, the existing
repository convention.

Incompatible combinations are refused, never silently merged: offsets in
`aim_at_reference` mode, or a non-zero `plate_angles_deg` in
`parallel_to_normal` mode, each raise `StudyConfigError`.

---

## Panel-local u/v offsets

The surface-local basis comes from the **existing** `target.normal` and
`target.tangent` keys — no new geometry keys:

| Axis | Definition | ISS panel |
| --- | --- | --- |
| `u` | `target.tangent`, longitudinal | +X, the 22 m dimension |
| `v` | `n x u`, transverse | +Y, the 12 m dimension |
| `n` | `target.normal`, toward the source | +Z |

The triad is right-handed (`u x v = n`), and `v` is exactly the binormal
`firing_plan.pose_for` already used for the second DCM column, so the pose
convention and the reporting convention share one basis.
`TargetSpec.local_basis()` returns it.

`sweep.source_offsets_u` and `sweep.source_offsets_v` sweep the source along
`u` and `v`. Both default to `[0.0]`, so an existing configuration's poses
are untouched. Cases are enumerated **distance-major, then u, then v, then
angle**:

```
for distance: for u_offset: for v_offset: for angle
```

which collapses to the historical "all angles at the first distance, then
all angles at the next" at the default offsets. For an offset sweep the case
count is exactly `n_distances x n_u_offsets x n_v_offsets`. Non-finite
offsets are rejected.

Case identifiers name what actually varies:

| Mode | Identifier |
| --- | --- |
| `aim_at_reference` | `case000_alpha0p0_d4` (unchanged) |
| `parallel_to_normal` | `case000_modelCollisionless_L4_u0_v0` |

Every offset value is recorded in each result (`source_offset_u`,
`source_offset_v`, `source_axis_mode`).

---

## Derived Knudsen metadata

**PyRPOD's plume models are collisionless, and the optional `knudsen` block
does not change that.** No solution, field value or surface load anywhere in
the pipeline reads Kn back; nothing is corrected for rarefaction. The block
exists so an analytical case can be *labelled* with the regime it is meant to
represent, which is what a later, entirely separate workflow needs to line
PyRPOD cases up with externally generated data.

`Kn = mean_free_path_m / L_ref`, with `L_ref` chosen by **exactly one** of two
mutually exclusive modes:

```yaml
knudsen:
  mean_free_path_m: 1.0
  reference_length: source_distance          # the case's own swept distance
  definition: lambda_over_source_distance
```

```yaml
knudsen:
  mean_free_path_m: 1.0
  reference_length_m: 0.5                    # a fixed length
  definition: lambda_over_nozzle_diameter
```

Rules, all enforced with a specific `StudyConfigError`:

- `mean_free_path_m` is **required**, positive and finite. It is **never
  inferred** from the gas properties in the thruster definition file — a
  free-molecular model carries no collision rate to infer it from.
- Exactly one reference-length mode. Supplying both, or neither, is an error.
- `reference_length` accepts only the symbolic value `source_distance`; any
  other fixed length goes in `reference_length_m`.
- `definition` is a free-text label, defaulted from the mode.

Omit the block entirely and every Knudsen field is simply absent — empty CSV
columns and `null` in JSON, never a fabricated value. The study metadata
records the block with an explicit `role` field stating it is derived
metadata only.

To sweep Kn itself, write **one study per label** (the mean free path is a
single scalar per study, not a swept axis); see the ISS-panel README.

---

## YAML configuration

A study configuration is a **layer on top of an existing PyRPOD case**. The
case's `config.ini` keeps owning the STL, TCF, TDF, plume model and gating
geometry, so every existing case and public API is unaffected; the YAML adds
only what a trade study needs and an INI cannot express.

```yaml
study:
  name: cai2016_flat_plate_baseline        # required; identifies every result row
  description: ...
  case_dir: ..                             # required; relative to THIS file
  output_dir: ../results/studies/flat_plate_baseline

thruster:
  id: T1                                   # optional; validated against the case TCF

plume_model:
  name: SimplifiedGasKinetics              # | CollisionlessGasKinetics
  parameters:                              # recorded for provenance only
    gas: argon
    speed_ratio_S0: 2.0

target:
  geometry_id: flat_plate_transformed.stl  # defaults to the case's [tv] stl
  reference_point: [0.0, 0.0, 0.0]         # sweep is built about this point
  normal: [0.0, 0.0, 1.0]                  # n: toward the plume-source side
  tangent: [1.0, 0.0, 0.0]                 # u: longitudinal (v = n x u)
  components:
    - name: plate
      selector: all                        # or face_indices: [...] / bounds: {...}

sweep:
  mode: per_case                           # per_case (default) | single_jfh
  source_axis_mode: aim_at_reference       # (default) | parallel_to_normal
  plate_angles_deg: [0.0]                  # 0 = head-on; aim_at_reference only
  source_distances: [4.0]                  # from `reference_point`
  source_offsets_u: [0.0]                  # panel-local; parallel_to_normal only
  source_offsets_v: [0.0]                  # panel-local; parallel_to_normal only
  n_firings: 1                             # EXACT JFH entries per pose
  firing_duration_s: 1.0
  thrusters: [1]
  # firings: [...]                         # optional explicit poses; see below

loads:
  moment_reference_point: [0.0, 0.0, 0.0]
  normalization:                           # optional; omit for no coefficients
    reference_area: 64.0
    reference_length: 4.0
    dynamic_pressure: 1.1044652197738332
    reference_heat_flux: 637.3520362956127

knudsen:                                   # optional; DERIVED METADATA ONLY
  mean_free_path_m: 1.0                    # required; never inferred
  reference_length: source_distance        # XOR reference_length_m: <number>
  definition: lambda_over_source_distance  # free-text label

output:
  vtk:      {enabled: true}
  surface_distribution:
            {enabled: false, subdir: distributions}
  summary:  {csv: case_results.csv, metadata: study_metadata.json}
  plots:    {enabled: false, subdir: plots, per_case_distribution: false}

reference:
  path: null                               # optional external reference data
  label: null

metadata:
  coordinate_system: case global frame
  units: {}                                # overrides the SI defaults
```

**Schema additions on this branch**, all optional and defaulted so every
existing YAML file parses and behaves unchanged:

| Key | Default | Meaning |
| --- | --- | --- |
| `plume_model.name` | `SimplifiedGasKinetics` | Now also accepts `CollisionlessGasKinetics`; selects the model that computes the field |
| `sweep.source_axis_mode` | `aim_at_reference` | `aim_at_reference` \| `parallel_to_normal` |
| `sweep.source_offsets_u` | `[0.0]` | Longitudinal source offsets (m) |
| `sweep.source_offsets_v` | `[0.0]` | Transverse source offsets (m) |
| `knudsen` | absent | Whole block; see [Derived Knudsen metadata](#derived-knudsen-metadata) |
| `knudsen.mean_free_path_m` | — | Required within the block; positive, finite |
| `knudsen.reference_length` | — | Only `source_distance`; XOR with the next |
| `knudsen.reference_length_m` | — | A fixed positive length; XOR with the previous |
| `knudsen.definition` | from the mode | Free-text label |
| `output.surface_distribution.enabled` | `false` | Export panel-local per-face CSVs |
| `output.surface_distribution.subdir` | `distributions` | Where they land |
| `output.plots.per_case_distribution` | `false` | Per-case panel pressure maps |

Validation is strict and specific: a missing case directory, a case without a
`config.ini`, an unknown plume model, an empty or non-positive sweep axis,
non-finite offsets, an unknown axis mode, incompatible pose definitions,
non-orthogonal target axes, duplicate component names, a non-positive
normalization value, an invalid or ambiguous `knudsen` block, or a firing
list whose length disagrees with `n_firings` each raise `StudyConfigError`
naming the offending key.

Angles, distances and offsets are geometry, not physics: `normal` and
`tangent` define both the plane the source is swept in and the surface-local
reporting basis, so a curved target simply supplies the axes its sweep
should use (see [Adding another validation
geometry](#adding-another-validation-geometry)).

---

## `n_firings` and prescribed firings

**`n_firings` is the exact number of Jet Firing History entries contributed
by each swept pose.** It is validated before anything is generated (positive
integer only — zero, negative, fractional and non-numeric values are
rejected), the generated sequence is asserted to match it, and the JFH read
back from disk is asserted to report the same count. In `per_case` mode that
is the exact length of every case's history; in `single_jfh` mode the one
shared history must hold exactly `len(poses) * n_firings` entries.

Two ways to supply the poses:

- **Generated** (the usual case): the pose comes from the swept
  `plate_angles_deg` / `source_distances` entry and is repeated for
  `n_firings` successive firing intervals — one JFH entry per firing, each of
  `firing_duration_s`.
- **Explicitly prescribed**: list the poses under `sweep.firings`. Their
  count must agree with `n_firings` (in `single_jfh` mode, with
  `len(poses) * n_firings`, the sequence being assigned to the poses in sweep
  order); a mismatch is an error, never a silent truncation or extension.

```yaml
sweep:
  n_firings: 2
  firings:
    - position: [0.0, 0.0, 4.0]
      dcm: [[0, 0, 1], [0, 1, 0], [-1, 0, 0]]   # first COLUMN = thruster axis
      thrusters: [1]
      duration_s: 0.5
    - position: [0.5, 0.0, 4.0]
      dcm: [[0, 0, 1], [0, 1, 0], [-1, 0, 0]]
```

The generated pose convention depends on `sweep.source_axis_mode` — see
[Source axis modes](#source-axis-modes-aiming-vs-translating). In the default
`aim_at_reference` mode, for a target reference point `C` with outward normal
`n_hat` and in-plane tangent `t_hat`:

```
d_hat(alpha) = cos(alpha) * n_hat + sin(alpha) * t_hat
source position = C + L * d_hat(alpha)
thruster axis   = -d_hat(alpha)                 (aimed at C)
```

`alpha = 0` is head-on. The JFH DCM carries the thruster axis as its first
**column**, which is what the strike pipeline reads as the plume normal. This
reproduces the committed sweep-JFH generators
(`case/plume/plume_flat_plate_sweep/jfh/generate_sweep_jfh.py`) to file
precision — a regression test pins all 95 poses of that file.

The dynamics-driven path was corrected to the same meaning:
`PlumeStrikeEstimationStudy.print_jfh_1d_approach_n_fire(..., n_firings=N)`
previously accepted `N` and ignored it. It now validates the count, stops the
approach simulation at exactly `N` entries, and raises if the approach
completes in fewer firings than requested.

---

## Integrated loads

`pyrpod/mdao/surface_loads.py` turns the strike pipeline's per-face arrays
into component-level resultants. Sign and direction conventions follow the
pipeline's own face-selection rule: a face is struck only when its unit
normal `n_hat` opposes the plume direction, i.e. points back toward the
source.

### Force

```
F_p,i = -p_i * A_i * n_hat_i                         (pressure, into the surface)
t_hat_i = normalize(u_hat_i - (u_hat_i . n_in_i) * n_in_i),  n_in_i = -n_hat_i
F_s,i = tau_i * A_i * t_hat_i                        (shear, along the tangential flow)
```

`u_hat_i` is the local flow direction — the unit vector from the plume source
to the face centroid, since the collisionless flow is radial. Faces at exactly
normal incidence have no tangential direction and contribute no shear force.
Reported: `pressure_force`, `shear_force`, `force`, `force_magnitude`.

### Moment

```
M_ref = sum (r_i - r_ref) x F_i
```

about the user-defined `loads.moment_reference_point`. The pressure and shear
contributions are reported separately (`pressure_moment`, `shear_moment`)
alongside the combined `moment` and `moment_magnitude`.

### Center of pressure

A unique three-dimensional center of pressure does not exist in general: the
resultant of a distributed load is a force **plus a couple**, and only the
moment component perpendicular to the force can be represented by shifting
the force's line of action. The definition used here is stated explicitly:

```
r_cop = r_ref + (F x M_ref) / |F|^2
```

— the point **on the line of action of the resultant force that lies closest
to the moment reference point**. The force-parallel moment component is
reported separately as `residual_couple`; no choice of `r_cop` can remove it.

Degenerate cases are handled explicitly and never return a misleading value:

| Situation | `center_of_pressure` | `center_of_pressure_status` |
| --- | --- | --- |
| Non-zero resultant | the point above | `ok` |
| Nothing loaded | `None` | `zero_load` |
| Resultant is cancellation noise (below `1e-6` of the summed face-force magnitudes, e.g. a closed target loaded symmetrically) | `None` | `ill_conditioned` |

`pressure_weighted_centroid` — `sum(p_i A_i r_i) / sum(p_i A_i)` — is
reported as an auxiliary, always-defined location; it coincides with the
classical center of pressure for a planar component under unidirectional
pressure.

### Panel-local resultants and their signs

The same force and moment vectors are also reported on the target's
surface-local basis `(u, v, n)` (see [Panel-local u/v
offsets](#panel-local-uv-offsets)). This is a pure projection — it adds no
physics — and `surface_loads.project_to_panel_frame()` computes it.

**The sign conventions are stated and tested explicitly:**

| Field | Definition | Sign |
| --- | --- | --- |
| `normal_force` | `-F . n` | **Positive when the load presses INTO the panel** (away from the plume source). The target normal points *toward* the source and the plume pushes against it, so a compressive impingement load is positive. A negative value would mean the resultant pulls the panel toward the source. |
| `local_force_u`, `local_force_v` | `F . u`, `F . v` | In-surface components of the same resultant. |
| `local_moment_u/v/n` | `M_ref . u`, `. v`, `. n` | The **same** moment vector the global-frame fields report, about the **same** reference point, projected on each axis. Positive follows the right-hand rule about that axis. |
| `center_of_pressure_u/v` | `(r_cop - r_ref) . u`, `. v` | Panel-local coordinates of the center of pressure, **measured from the moment reference point** (the panel center in the ISS-panel studies). |

Worked moment sign: a pressure patch centred at `+u` pushes along `-n`, so
`M = (u_off * u) x (-F_n * n) = +u_off * F_n * v` (using `u x n = -v`). **A
source displaced toward `+u` therefore gives a positive `local_moment_v`**,
growing with the offset until the patch starts leaving the panel. Likewise a
source displaced toward `+v` gives a **negative** `local_moment_u`.

When no center of pressure exists (`zero_load` / `ill_conditioned`),
`center_of_pressure_u/v` are `None` rather than zero.

### Thermal quantities

Instantaneous and peak quantities only (time-integrated heat dose is
deliberately out of scope on this branch): `max_heat_flux` (peak per-face
heat flux, W/m^2) and `total_heat_load` (`sum q_i A_i` over the component, W),
plus `max_pressure`, `max_shear_stress`, `affected_area` and `struck_faces`.
The per-face heat-flux load the pipeline already accumulates stays in the VTK
output.

### Coefficients

Coefficients are computed **only** when the configuration supplies every
input a given coefficient needs. Nothing is defaulted or invented.

| Coefficients | Required inputs |
| --- | --- |
| `Cp_max`, `Cf_max` | `dynamic_pressure` |
| `Cq_max` | `reference_heat_flux` |
| `CF`, `CFx`, `CFy`, `CFz` | `dynamic_pressure`, `reference_area` |
| `CM`, `CMx`, `CMy`, `CMz` | + `reference_length` |

With no `normalization` block at all, `coefficients` is empty and
`coefficients_available` is `False` — which is what the cylinder example
does, since no reference values exist for it yet.

---

## Result schema

`StudyResults` holds one `CaseResult` per **case x component x firing**, each
carrying enough metadata to reproduce and later compare the calculation. In
`per_case` mode there is one `case_id` per pose; in `single_jfh` mode every
record shares the study's single `case_id` and is distinguished by
`firing_id`, but both carry the same `plate_angle_deg` / `source_distance`,
so downstream code (CSV, plots, reference comparison) is mode-agnostic:

- **identity** — `study_name`, `case_id`, `component`, `firing_id`
- **geometry** — `geometry_id`, `mesh_faces`, `component_faces`,
  `component_area`, `coordinate_system`, `units`
- **pose and sweep** — `plume_source_position`, `plume_source_orientation`
  (9 DCM values), `target_normal`, `target_tangent`,
  `target_reference_point`, `plate_angle_deg`, `source_distance`,
  `source_offset_u`, `source_offset_v`, `source_axis_mode`,
  `firing_duration_s`, `thrusters`
- **model** — `plume_model`, `plume_model_parameters`, and the derived
  `model_variant` (`Simplified` / `Collisionless`, for plot legends and
  grouping)
- **loads** — `pressure_force`, `shear_force`, `force`, `force_magnitude`,
  `moment_reference_point`, `pressure_moment`, `shear_moment`, `moment`,
  `moment_magnitude`, `center_of_pressure`, `center_of_pressure_status`,
  `residual_couple`, `pressure_weighted_centroid`
- **panel-local loads** — `normal_force`, `local_force_u`, `local_force_v`,
  `local_moment_u`, `local_moment_v`, `local_moment_n`,
  `center_of_pressure_u`, `center_of_pressure_v` (see [Panel-local
  resultants and their signs](#panel-local-resultants-and-their-signs))
- **surface fields** — `max_pressure`, `max_shear_stress`, `max_heat_flux`,
  `total_heat_load`, `affected_area`, `struck_faces`
- **derived Knudsen metadata** — `knudsen_number`, `mean_free_path`,
  `knudsen_reference_length`, `knudsen_definition`
- **coefficients** — `coefficients`, `coefficients_available`
- **artifacts and provenance** — `vtk_path`, `jfh_path`,
  `surface_distribution_path`, `config_path`, `case_dir`, `code_version`
  (git commit when available), `generated_at`

**Every field added on this branch is optional or safely defaulted**, so an
older result file still loads and the reference-comparison API is unchanged.
Offsets default to `0.0` and the axis mode to `aim_at_reference`; the
panel-local, Knudsen and distribution fields default to `None` and serialize
as **empty CSV columns**, never as a fabricated zero.

`CaseResult.quantity()` exposes the new comparable scalars —
`normal_force`, `local_force_u/v`, `local_moment_u/v/n`,
`center_of_pressure_u/v` and `knudsen_number` — alongside the existing ones,
returning `None` when a quantity is unavailable for that record.

Two machine-readable artifacts are written, in formats the repository already
uses (no Parquet, no new dependency):

- **CSV** (`StudyResults.write_csv`) — one flat row per record; vectors are
  expanded to `<name>_x/_y/_z` columns and each coefficient gets its own
  `coeff_<name>` column, so the file is directly plottable. **No per-face
  array is ever embedded in a CSV row** — those go to the VTK files and the
  distribution CSVs;
- **JSON** (`StudyResults.write_metadata`) — `schema`, study-level
  `provenance` (including the plume model, the source axis mode, the panel
  basis, the Knudsen block with its explicit "derived metadata only" role,
  mesh size, component list, code version and known limitations) and the
  nested per-case records. This is the document an externally generated
  dataset is later transformed into for comparison.

---

## VTK outputs

Per-face data are never reduced away. With `output.vtk.enabled` (the default)
the standard pipeline writer produces one `.vtu` per firing, laid out
according to the sweep mode:

- `per_case`: `<output_dir>/cases/<case_id>/results/strikes/firing-<i>.vtu`
- `single_jfh`: `<output_dir>/results/strikes/firing-<i>.vtu` — one numbered
  series over the whole sweep, which ParaView opens as a time sequence

Both carry the pipeline's own cell fields:

```
strikes, cum_strikes, pressures, max_pressures, shear_stress, max_shears,
heat_flux_rate, heat_flux_load, cum_heat_flux_load
```

Every case result advertises its own `vtk_path`, and that path is also in the
CSV and JSON summaries, so artifacts are discoverable from the returned
results. Set `output.vtk.enabled: false` for a fast, artifact-free run (the
numbers are identical; only the files are skipped).

---

## Panel-local surface distributions

**The VTK files remain the primary full-resolution visualization output.**
Enabling `output.surface_distribution.enabled` adds a *second view of the
same numbers*: one CSV per case and component, holding every face of that
component in the target's own `(u, v)` coordinates.

The reason is practical. A `.vtu` is the right artifact for ParaView and the
wrong one for a plotting script, a spreadsheet, or a later comparison
workflow that would otherwise need a VTK reader.

Columns (stable order, `pyrpod/mdao/surface_distribution.py`):

```
face_index, centroid_x, centroid_y, centroid_z, local_u, local_v,
area, pressure, shear_stress, heat_flux, strike_count
```

- `local_u` / `local_v` are `(centroid - target.reference_point)` projected
  on the target's `u` and `v` axes.
- `face_index` is the index into the **full** target mesh, so a row traces
  back to both the mesh and the VTK file even for a component subset.
- **No interpolation, resampling, smoothing or structured-grid projection.**
  Every row is one native mesh face carrying the value the strike pipeline
  computed for it. An unstructured triangle mesh is exported as unstructured
  triangles.
- **No common-grid projection onto any external mesh.** Producing a shared
  grid is a separate workflow's job and is not implemented here.

Units and case metadata go in a sidecar `<name>.meta.json` beside the CSV —
column units, the panel basis, the pose, the plume model and any derived
Knudsen fields — so a distribution file is self-describing without the study
metadata document. The CSV path is recorded in
`CaseResult.surface_distribution_path`, and therefore in the summary CSV and
JSON too.

Files land in the per-case directory
(`<output_dir>/cases/<case_id>/distributions/`) for `per_case` mode, and in
`<output_dir>/distributions/` for `single_jfh`.

---

## Optional plots

Plot generation is entirely optional — no automated test requires graphical
output — and matplotlib's non-interactive `Agg` backend is pinned when a
plotting module is imported (lazily, only when plots are asked for), so a
headless or CI run never opens a window. No plotting dependency beyond
matplotlib is introduced. Enable with `output.plots.enabled: true` or call
`study.plot()`.

`study_runtime.study_plots_for()` picks the families a given study needs, so
neither engine nor the `TradeStudy` façade holds plotting logic.

**Angle sweeps** (`pyrpod/mdao/study_plots.py`, always produced):

```
force_vs_angle.png        moment_vs_angle.png      heat_flux_vs_angle.png
force_vs_distance.png     moment_vs_distance.png   center_of_pressure.png
reference_comparison.png  (when a comparison report exists)
```

**Offset sweeps** (`pyrpod/mdao/panel_plots.py`, added when the study sweeps
offsets or uses `parallel_to_normal`). Every series is grouped by stand-off
distance **and** plume model, so results merged from a Simplified run and a
Collisionless run of the same geometry plot as separate labelled series
rather than being averaged:

```
normal_force_vs_offset_u.png     moment_v_vs_offset_u.png
peak_pressure_vs_offset_u.png    cop_u_vs_offset_u.png
normal_force_vs_distance.png
```

When `v` offsets are swept, the analogous transverse figures appear too
(`normal_force_vs_offset_v.png`, `moment_u_vs_offset_v.png`,
`cop_v_vs_offset_v.png`).

**Per-case panel pressure maps** (`output.plots.per_case_distribution: true`,
which also needs `output.surface_distribution.enabled`):
`panel_pressure_<case_id>.png` — panel-local `u` horizontal, `v` vertical,
pressure in Pa, with the panel edges drawn and the plume-centerline
intersection marked. A flat plate is an **unstructured triangle mesh**, so no
structured grid is fabricated: the field is drawn as a Delaunay
triangulation of the face centroids, falling back to a face-coloured scatter
when the face count is too small for contours to be honest.

---

## External reference data

`pyrpod/mdao/reference_data.py` compares study results against independently
generated data **without knowing where it came from**. A reference record is
just named quantities attached to matching keys; DSMC, an analytical
solution, an experiment and another code are all handled identically, and no
producer-specific importer exists.

Supported formats:

- **CSV** — one row per record. `case_id`, `component`, `plate_angle_deg` and
  `source_distance` are recognized as keys; every other numeric column is a
  quantity, and `<name>_x/_y/_z` triplets are folded into vectors. This is
  exactly the layout `StudyResults.write_csv` emits, so an external producer
  can be transformed into it column-for-column.
- **JSON / YAML**:

```yaml
label: DSMC run 12
source: /runs/dsmc-12
units: {force: N}
records:
  - key: {plate_angle_deg: -40.0, source_distance: 2.0, component: plate}
    quantities:
      force: [1.117, 0.0, -2.376]
      max_pressure: 0.34
      center_of_pressure: [0.031, 0.0, 0.015]
```

```python
from pyrpod.mdao.reference_data import load_reference_dataset

report = study.compare(load_reference_dataset('dsmc_run_12.json'))
print(report.max_relative_error('force'))
report.write_csv('comparison.csv')
```

Metrics are applied only where they are mathematically meaningful: absolute
error, relative error (`None` when the reference is zero, never infinity),
normalized RMSE (by reference range or mean), peak-value error,
integrated-load error (compared as vectors, so a load of the right magnitude
pointing the wrong way is an error) and center-of-pressure displacement. A
quantity the reference supplies but the result does not is reported as
`missing_candidate`; a case with no matching record is listed in
`unmatched_cases`. Nothing is fabricated or defaulted.

For the flat-plate tests the independent reference is the exact Cai 2016
solution (`pyrpod/plume/CaiImpingement2016.py`, Eq. 15 quadrature),
converted to dimensional loads with the case's own normalization — PyRPOD's
own output is never used as the validation baseline. Measured agreement over
the sweep: integrated normal load within ~1.5% mean / 2.2% max, component
heat load within ~15%, matching the documented accuracy of the Maxwellian
engineering chain against the exact collisionless solution.

### Adding DSMC results later

> **No DSMC handling exists in PyRPOD, and none is added by this branch.**
> Nothing here launches OpenFOAM, writes OpenFOAM dictionaries, imports or
> parses DSMC fields, interpolates a DSMC mesh, manages DSMC jobs, digitizes
> published DSMC data, or produces a DSMC comparison report. Both plume
> models are **collisionless**: there is no collisional correction, no wake,
> shadowing or secondary-collision physics, and the Knudsen number is derived
> metadata that no solution reads. What PyRPOD produces is **analytical
> datasets**; comparing them with externally generated DSMC results is a
> separate workflow.

Nothing DSMC-specific is needed. Transform the external results into either
supported format, keyed by `plate_angle_deg` / `source_distance` /
`component` (or `case_id`, which the study metadata records), then call
`study.compare(...)` or point `reference.path` at the file. Because the
result schema records the mesh, coordinate system, units, source pose, model
and configuration provenance, the transformation can be done from the stored
metadata without rerunning PyRPOD.

---

## Adding another validation geometry

1. Create (or reuse) a PyRPOD case directory with its `config.ini`, target
   STL, TCF and TDF — the study layer adds nothing here.
2. Add a study YAML next to it under `study/`, setting `target.reference_point`
   to the geometry's reference (plate center, cylinder centroid, ...) and
   `target.normal` / `target.tangent` to the plane the source should be swept
   in.
3. Break the target into components if the loads should be reported
   separately, using `face_indices` or a `bounds` box on face centroids.
4. Supply `loads.normalization` only if valid reference values exist for that
   geometry; otherwise leave it out and the coefficients are reported as
   unavailable.
5. Point `reference.path` at external data whenever it exists.

`case/plume/plume_cylinder_sweep/study/cylinder_baseline.yaml` is a worked
example on a curved, closed target. Its results are a pipeline smoke case,
not a validated physical answer — no cylinder reference data exists yet — and
it deliberately supplies no normalization inputs.

---

## ISS-panel example

`case/plume/iss_panel_thesis/` is the worked example of everything above: an
idealized **22 m x 12 m** flat panel standing in for one ISS solar-array
wing, with the source translated parallel to it. See that directory's
`README.md` for the full conventions, expansion guidance and Knudsen labels.

```bash
# Centered source at L = 4 m, SimplifiedGasKinetics (seconds)
python case/plume/iss_panel_thesis/run.py baseline-simplified

# The same case with CollisionlessGasKinetics (~20 s, 2112 faces)
python case/plume/iss_panel_thesis/run.py baseline-full-cai

# 3 distances x 5 longitudinal offsets = 15 cases (~20 s)
python case/plume/iss_panel_thesis/run.py sweep

# All three, into a scratch tree, with progress logging
python case/plume/iss_panel_thesis/run.py all \
    --output-dir /tmp/iss_panel --verbose

# Skip figure generation
python case/plume/iss_panel_thesis/run.py sweep --no-plots
```

Or from Python, through the ordinary package API:

```python
from pyrpod.mdao.TradeStudy import TradeStudy

study = TradeStudy.from_config(
    'case/plume/iss_panel_thesis/study/iss_panel_offset_distance_sweep.yaml')
results = study.run()
paths = study.plot()
```

The two baseline configurations differ in **exactly one line** —
`plume_model.name` — so the pair isolates the far-field simplification: 2.824
N vs 2.802 N centered normal force, 0.3317 Pa vs 0.3136 Pa peak pressure.
That the numbers differ at all is what demonstrates that model selection
reaches the calculation rather than only the metadata.

An **independent** analytical cross-check ships alongside it:
`case/plume/iss_panel_thesis/cai2016_reference.py` evaluates the exact Cai
2016 surface solution at the same face centroids and exports the same
distribution schema. It is a reference generator, **not** a plume-model
backend — PyRPOD never imports it and `PlumeStrikeCalculator` cannot reach
it. At `L = 4 m` the peaks agree to 2.8% in pressure, 3.4% in shear and 7.8%
in heat flux, which is the Maxwellian wall chain against a direct
integration of the wall fluxes.

---

## Known limitations

- **No shadowing or occlusion.** Face selection is the existing pipeline
  behavior: a face is struck when it lies inside the plume wedge/radius and
  its normal faces the source. Plume shadowing, occlusion by intervening
  geometry, self-shadowing of concave or closed targets and back-facing
  surfaces are **not** modeled, so a face hidden behind other geometry still
  receives load. This branch changes nothing about that; results on closed
  targets must be read with it in mind.
- **Edge-on poses are degenerate.** At exactly `+/-90 deg` the facing test's
  dot product is zero in exact arithmetic, so float32 STL normals decide
  strike membership arbitrarily. The plate-averaged load there is ~zero, but
  individual near-nozzle faces can carry large grazing values. The committed
  sweep keeps those poses; the symmetry checks exclude them.
- **Cumulative fields mean different things per mode.** In `single_jfh` mode
  they are a genuine sweep envelope; in `per_case` mode they only restate the
  case they belong to, since each case's history starts fresh. Read
  `max_pressures` / `cum_strikes` in a per-case VTK accordingly.
- **Coefficients require explicit normalization.** By design: nothing is
  inferred from the geometry.
- **One thruster per study.** The workflow prescribes a single firing source.
  Multi-thruster and multi-group behavior is untouched elsewhere in PyRPOD
  but is not exercised here.
- **Both plume models are collisionless.** `SimplifiedGasKinetics` and
  `CollisionlessGasKinetics` are free-molecular: no intermolecular
  collisions, no continuum or transitional correction, no wake, and no
  secondary-collision physics. A configured Knudsen number is **derived
  metadata** and never enters the solution, so labelling a case `Kn = 0.01`
  does **not** make it a continuum result — it records the regime the case is
  intended to represent. A study is only physically meaningful where the
  free-molecular assumption holds.
- **No DSMC.** No OpenFOAM execution, dictionary generation, field import,
  mesh interpolation, job management or DSMC comparison report exists in
  PyRPOD. The outputs are analytical datasets for a later, separate
  comparison workflow.
- **Distribution CSVs are native faces only.** No interpolation and no
  common-grid projection onto an external mesh; producing a shared grid is
  out of scope here.
- **Cylinder validation is not quantitative yet.** The architecture and the
  comparison interface are ready for it; the reference data is not.
- **Heat loading is instantaneous/peak only.** Time-integrated heat dose and
  impulse are deliberately out of scope on this branch.
