# MPI EB spray viscosity — 2026-09-25

Code `3a14a7e`, parent `747362a`. The MPI EB spray coupling API accepts optional
finite nonnegative `difmag`, collectively validates it and requires exact rank
agreement before private particle/film work. It forwards the resolved value
(default zero) to the gas full step. Both application impact branches now pass
the tree input coefficient. Gas patch checkpoints already bind that value;
changing it rejects gas snapshot restoration.

## Completed checks

GNU Fortran 13.3.0, OpenMPI 4.1.6, Python 3.12.3, Ninja Debug; reused
`build/spray-amr-reactive-debug`. Fixed MPI executable built. Focused CTest
**1/1 passed**, 55.42 s (55.46 s elapsed), no failures or skips. Compiler warnings
remain. Provenance, Python syntax, project contract and whitespace checks passed.

The coefficient 0.037 fixture uses a level-set cut wall, Mach inlet, automatic
timestep selection, recursive gas stepping, implicit chemistry, transport and
particle wall reflection. Existing 2/4-rank gas/particle comparisons pass.
Two-rank checkpoint continuation on four ranks matches uninterrupted gas and
sorted particle CSV records exactly; inventory ledgers satisfy the existing
3e-12 relative bound. Changed coefficient, negative and NaN inputs reject
before result CSV publication. Existing immutable checkpoint, missing restart,
changed context/inlet and truncated-particle checks also pass.

Rank-inconsistent coefficients were not separately injected into the coupling
API here. Selected, film-impact and common-step positive configurations were
not independently run. This is bounded application/restart evidence, not
external PeleC parity, nozzle qualification or physical validation. EB default
policy and coarse/fine geometry fidelity remain open.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_eb_spray_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_difmag$' --output-on-failure
python3 -m py_compile tests/application/check_mpi_eb_spray.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/mpi-eb-spray-difmag-checks/build.log`, `LastTest.log`;
fixture directory `build/spray-amr-reactive-debug/tests/spray_mpi_eb_difmag/eb_spray_6t1eb_jm/`.
