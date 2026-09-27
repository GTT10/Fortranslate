# Reference KHRT coupling and application input — 2026-09-27

Implementation `2bc627e1`, based on `2f0a509c`, local main, version 0.246.0.

## Change

`spray_options%khrt_parcel_factor` selects the child policy: -1 retains the
four-direction quadrature; [0,1] selects reference count and sampled azimuth.
It requires KHRT children. The release APIs take the option by default and
retain their explicit factor override. A reference-policy call without a gas
momentum output is rejected; default calls may request a zero momentum array.

Serial single/mixed regular coupling and serial EB, MPI regular, AMR/EB and
recursive SDC source paths now deposit returned momentum together with energy
at parent source locations. Recursive virtual/ghost populations and non-moving
stages disable the child policy along with child creation. The configured
regular constructor checks variable-count capacity and ID overflow before
publication. Existing candidate-state transactions contain failed deposition.

Regular serial/MPI input accepts flat `khrt_parcel_factor`. Serial EB and MPI
AMR use `options%khrt_parcel_factor`; MPI EB tree uses
`particle_options%khrt_parcel_factor` in `eb_spray`. Enabled policy joins the
checkpoint context. Disabled derived-type serialization omits the new field,
preserving the preceding canonical context. Particle IDs and the global ID
watermark supply the deterministic sample key; no rank-local RNG state is added.
These count/basis formulas follow the pinned reference, while the random stream
and conservative gas compensation remain PeleF choices.

## Independently completed checks

Existing Ninja Debug fixed-mechanism build `build/spray-amr-reactive-debug`,
GNU Fortran 13.3.0, Open MPI 4.1.6 (`/usr/bin/mpifort`, `/usr/bin/mpiexec`).

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_eb_spray_3d pelef_spray_3d pelef_mpi_spray_3d test_spray_mpi --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_(reference_khrt_app_(common|recursive)|mpi_migration_4|mpi_eb_sdc_recursive_restart_regrid)$' --output-on-failure
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_reference_khrt_app_(common|recursive)$' --output-on-failure
python3 -m py_compile tests/application/check_reference_khrt_app.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final application run: **2/2 passed**, common 15.51 s, recursive 31.57 s,
total 47.14 s. The preceding invocation also passed MPI migration/release
(11.31 s) and legacy recursive regrid/restart (73.88 s). Its two new application
checks failed because no children formed in the initial short, low-density
fixture. The final fixture uses the existing KHRT-case density 1 kg/m³ and
pressure 435483.87096774194 Pa; the test also uses exported IDs to identify
new children, since this CSV has no `parent_id` column. No tolerance changed.

Both final cases start with two water parcels, use factor 0.8, and finish with
26 records. Chemistry and droplet drag/heat/evaporation are disabled to isolate
breakup; gas molecular transport remains enabled. Combined gas/particle
inventory residuals are 8.3697297422764337e-16 (common) and
6.4711408465033872e-16 (recursive), below 3e-12. A step-1 checkpoint on 2 ranks
continues on 4 ranks with identical gas/budget CSV bytes and sorted particle
rows to continuous execution. Changing the factor or returning to -1 rejects
the checkpoint context. The symmetric-policy population differs as expected.
The legacy regrid/restart checkpoint contains no `KHRT_PARCEL_FACTOR` field.

Build, provenance, Python syntax, contract and whitespace checks passed. The
initial compile hit a long namelist line, fixed by wrapping it. Existing
real-comparison, static-storage and executable-stack warnings remain.
Mechanism source SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
Artifacts: `build/khrt-reference-app-checks/` includes initial/final logs and
`LastTest.log`. Final inputs/CSVs/checkpoints:

- `build/spray-amr-reactive-debug/reference-khrt-common/khrt_28huxpny/`
- `build/spray-amr-reactive-debug/reference-khrt-recursive/khrt_aj0fmnft/`

## Limits

The new mode was executed through fixed-mechanism MPI EB SDC common/recursive
applications. Other coupled routes were compiled, not independently exercised
with this option. Selected mechanisms, new-mode regridding/films/LES/active
chemistry, random-distribution quality, installed binaries, full configuration
matrix, external full-solver parity and physical validation were not checked.
Legacy canonical-field omission was inspected; an old binary's checkpoint was
not independently restored in this increment.
