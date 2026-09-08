# Build-time selected 0D reactors

These inputs exercise `pelef0d_selected`, which exists only when CMake is
configured with `PELEF_MECHANISM_BUNDLE`. Composition entries are mapped by
exact species name, so their input order does not need to match the bundle.

The fixture input belongs to `tests/fixtures/mechanism_bundle.json`. The full
input belongs to `mechanisms/h2o2_full.json` and its pinned Cantera YAML.
`short_final_step.nml` verifies that the scheduler can land exactly on a final
time whose last requested fragment is shorter than the adaptive minimum; the
fragment must be accepted whole rather than adaptively reduced again.

`integrator` defaults to `implicit`. The `_cvode.nml` inputs require a build
configured with `PELEF_ENABLE_SUNDIALS=ON` and the pinned SUNDIALS 7.2.0
official Fortran modules. A non-SUNDIALS build rejects those inputs before
creating the CSV. `maximum_steps` applies to the native adaptive loop;
`cvode_max_internal_steps` is CVODE's cumulative internal-step limit.

The CVODE cases use the same output schedule and CSV schema. The full H2/O2
case requests tighter mass-fraction tolerances than the native case because
CVODE's scalar tolerances control only the reduced differential state;
temperature is recovered algebraically from fixed internal energy.

CSV `wdot_<species>` columns are molar production rates in `kmol/(m^3 s)`.
