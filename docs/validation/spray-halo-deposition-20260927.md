# Matched particle halo sampling and deposition — 2026-09-27

Implementation `c2754c6b`, based on `8c5124a5`, local `main`, version 0.246.0.
PeleC reference `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
H2O2 provenance preflight SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

## Operation

`deposit_mpi_amr_particle_halo` shares the sampler's CIC weights, fluid-fraction
normalization and optional EB centroid visibility. Extensive particle inputs
become fresh source densities per fluid volume in the supplied halo layout.
The gas halo remains read-only. Source fields publish only after collective
success; invalid support or arithmetic overflow retains the previous output.
The sampler also accepts paired `contributions` and `source_halo` arguments
for combined sampling/deposition, with joint publication of both outputs.

Source halos can feed `sum_mpi_amr_particle_source_halo`; this change does not
wire particle evolution, population accounting or recursive time scheduling
into the application. It establishes matched spatial operators, not full
PeleC algorithm parity or physical validation.

## Completed verification

Fixed-mechanism Ninja Debug build `build/spray-amr-reactive-debug`, GNU Fortran
13.3.0, OpenMPI/OpenRTE 4.1.6 through `/usr/bin/mpifort` and `/usr/bin/mpiexec`.

`spray_mpi_amr_sampling_2` and `_4` passed (2/2, 0.49 s and 0.53 s;
1.06 s total). The added assertions use two halo queries, each contributing
[2,-3], and fluid fraction 0.5. The global fluid-volume integral equals [4,-6]
within 1e-13. The componentwise inner product of source density and gas field,
weighted by fluid volume, equals the particle contribution/sample product
within 1e-13. Removing all fluid support rejects and retains the saved source
fields exactly. A manufactured slab blocks donors across the solid while
preserving the total deposited inventory within 1e-13. Existing native/packed
halo sampling, periodic-image queries and rejection tests also execute.

Provenance, target build, Python syntax, project contract and whitespace checks
passed. No build/test failures occurred in this change. No full matrix,
installed-binary audit, recursive application run or physical validation was
performed.

## Reproduction

Artifacts: `build/spray-halo-deposition-checks/`, including `test-detail.log`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_amr_sampling --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_amr_sampling_[24]$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
