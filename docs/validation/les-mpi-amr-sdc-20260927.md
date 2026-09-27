# MPI AMR LES SDC — 2026-09-27

Implementation `be7b4931`, based on `687d4cc9`, local main, version 0.246.0.

## Change

`compute_mpi_amr_tree_les_rhs_3d` returns instantaneous owner-local LES rates
and face flux densities using existing MPI ghost exchange and static/dynamic
LES kernels. It supports level selection, parent old/new packed endpoints and
stage fractions. Inactive owned rates are zero; fluxes allocate only on owned
active nodes. The maximum diffusivity is reduced across all ranks. Controls,
boundaries and temporal-donor presence are checked collectively. Outputs remain
private until success. The RHS does not advance, average down or reflux state.

`advance_mpi_amr_tree_les_sdc_3d` connects this evaluator to common-step and
recursive MPI gas SDC. Native molecular transport can be enabled alongside LES.
The existing SDC driver integrates and combines fluxes and performs coarse/fine
synchronization. The wrapper checks rank agreement on common/recursive selection
and validates optional count storage collectively. Diffusivity and counts are
published only with an accepted step.

## Independently completed checks

Artifacts: `build/les-mpi-amr-sdc-checks/` (`provenance.log`, `build.log`,
`build-test.log`, `ctest.log`, `LastTest.log`). Existing fixed-mechanism Ninja
Debug configuration, GNU Fortran 13.3.0 and `/usr/bin/mpiexec` Open MPI 4.1.6.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_spray_tree_runtime --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_les_mpi_amr_sdc --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^les_mpi_amr_sdc_(common|recursive)_4$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**2/2 passed**, no failures or skips: common 14.92 s, recursive 29.96 s,
CTest total 44.93 s. Builds, mechanism provenance, Python syntax, project
contract and whitespace checks passed. Existing real-comparison/conversion
warnings and executable-stack linker warning remain.
Pinned mechanism source SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

Each test uses a periodic two-level tree with 4 cubed root cells, refinement
ratio 2, four MPI ranks and two idle ranks. Smagorinsky, WALE, Vreman and Dynamic
all produce positive diffusivity. Molecular transport is enabled except in the
Vreman case. MPI state, temperature, boundary/source inventories, correction,
limiter, diffusivity and level counts match the serial wrapper exactly. The
serial composite inventory check uses a 3e-12 normalized tolerance; both runs
printed residual zero. An invalid dynamic filter on the last, idle rank rejects
the collective step and retains state, temperature, inventories, diffusivity
and counts.

## Limits

This is focused library integration evidence. Active chemistry, physical
boundaries, deeper trees, selected mechanisms, application routing, installed
binaries, a full configuration matrix, external PeleC parity and physical LES
validation were not exercised. Existing SSPRK2 LES operators were not changed.
