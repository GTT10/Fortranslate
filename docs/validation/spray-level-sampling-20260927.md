# Level-bounded particle gas sampling — 2026-09-27

Implementation `23cdce8a`, based on `1d2e1d62`, local `main`, version 0.246.0.
PeleC reference `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
H2O2 provenance preflight SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

## Change

Regular and EB cloud stencils accept optional `maximum_level`. The containing
cell/cloud width is selected at or below this ceiling. Deeper nodes are skipped
and do not mask their parents. Existing overlap, periodic wrap and EB
fluid-volume/visibility weighting then operate on this capped composite view.
The serial regular sampler and MPI regular/EB sampler forward the limit.
MPI agrees on the effective limit, including empty-query ranks, and fetches
only cells in the resulting stencils. Failed requests retain sample output.
The default remains the original finest-grid composite interpolation.

This is a gas view for coarse virtual particles. It does not reconstruct
fine-level ghost cells or interpolate parent states in time: outside the capped
fine patches the view still uses coarse cells. Separate level source deposition,
transfer and recursive spray application integration remain open. Physical
feedback APIs retain their existing composite path. No upstream interpolation
or trajectory parity is claimed.

## Completed checks

Fixed-mechanism Ninja Debug build `build/spray-amr-reactive-debug`,
GNU Fortran 13.3.0, OpenMPI/OpenRTE 4.1.6 through `/usr/bin/mpifort` and
`/usr/bin/mpiexec`.

Six CTests passed, 3.03 s final run:
`spray_mpi_amr_sampling_2`, `_4`, `spray_mpi_eb_sample_2`, `_4`,
`spray_mpi_eb_periodic_sample_2`, `_4`.
Regular sampling uses distinct constant fields at three levels and checks
analytic coarse/intermediate values, exact serial/MPI equality, ignored NaNs
in unrequested deeper fields, invalid-level rollback and rank-disagreeing
ceilings. Serial empty queries check valid and invalid ceilings. EB tests
verify coarse sampling returns the covered parent's value at cut/periodic
interfaces. Existing composite interpolation and collective rejection cases
also execute. The final rerun followed an empty-query validation fix found
during review; no build/test failures occurred.

Provenance, target builds, Python syntax, project contract and whitespace checks
passed. Existing compiler warnings remain. No full matrix, selected-mechanism
build, installed audit, recursive application run or physical validation was
performed.

## Reproduction

Artifacts: `build/spray-level-sampling-checks/`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_amr_sampling test_spray_mpi_eb_sample --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_(amr_sampling|eb_(periodic_)?sample)_[24]$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
