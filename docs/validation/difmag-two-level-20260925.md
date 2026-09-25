# Two-level AMR artificial-viscosity propagation — 2026-09-25

## Implementation

Code `ec46e26`, parent `8b05c74`. Serial and sparse MPI two-level hydro,
reaction/hydro Strang splitting, and reaction/transport/hydro composition accept
a trailing optional `difmag`. The value reaches periodic/physical coarse CTU
and each fine substep's halo CTU before integrated fluxes enter reflux. Omitted
values retain the existing solver default. Serial hydro rejects nonfinite or
negative values and positive coefficients outside PPM. MPI hierarchy, physical
coarse and periodic coarse entry points resolve the default for collective
control comparison, enforcing coefficient agreement before publishing state.

Application input and checkpoint binding for this older two-level family are
not enabled yet. Tree gas/spray application support has separate earlier
records. No PCM/PLM artificial-viscosity kernel was added here.

## Completed checks

GNU Fortran 13.3.0, OpenMPI 4.1.6, Ninja Debug, Python 3.12.3. Existing fixed
mechanism, tests-enabled build `build/spray-amr-reactive-debug`.

**2/2 passed**, sparse hierarchy test 20.15 s and coarse/halo test 2.26 s,
22.45 s total elapsed. Build, provenance, Python syntax, project contract and
whitespace checks passed. No failed checks or skips in this increment;
compiler warnings remain.

The existing eight-rank PPM fixtures now pass coefficient 0.037 into serial and
MPI references, including periodic coarse/fine hierarchy, physical boundary
patch placements and reacting molecular-transport composition. Coarse/fine
conserved state, temperature and maximum reflux correction match serial
exactly. Coarse face fluxes also match exactly. Existing boundary ledger,
composite conservation and elemental-balance tolerances remain unchanged.
One-rank coefficient mismatch rejects with zero reported correction and
bit-identical hierarchy state. Existing CFL and invalid-boundary rollback
checks remain active. Empty fine ownership remains covered by the sparse
fixture layout.

No selected build, application/restart run, full matrix, installed audit,
external PeleC oracle or physical validation was performed.

## Commands and artifacts

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_amr_reactive_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^mpi_amr_ppm_(coarse|hierarchy)_3d_np8$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Logs: `build/difmag-two-level-checks/{build.log,test.log,LastTest.log}`.
