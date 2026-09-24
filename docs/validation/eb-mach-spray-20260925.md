# EB spray Mach inlet coupling — 2026-09-25

Implementation `51ef939`, base `e1fa35f`.
`advance_mpi_eb_spray_reactive` now forwards optional six-face
Mach, temperature and pressure controls to its gas reaction/transport/hydro
stage. Both application dispatch branches (with and without wall-impact film
state) forward them. The application rejection of active Mach inlets is removed.
Automatic timestep and checkpoint envelope binding already use the same gas
controls. Particle boundary rules and film evolution are unchanged; a gas inlet
does not implicitly select a different particle boundary model.

## Focused verification

GNU Fortran 13.3.0, Open MPI 4.1.6 with `/usr/bin/mpiexec`, Ninja Debug,
Python 3.12.3. Fixed MPI EB spray and selected H2/O2 MPI EB spray built in their
existing separate directories. Commands:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_eb_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_reactive_eb_tree_3d_selected_spray --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_mach$' --output-on-failure --output-log build/eb-mach-spray-checks/ctest.log
python3 -m py_compile tests/application/check_mpi_eb_spray.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

The level-set EB fixture enables high-y Mach 0.1 at 400 K and 120000 Pa,
automatic timesteps and recursive gas stepping. It uses the existing
gas/particle 2/4-rank comparison, wall reflection, 2-to-4-rank restart, coupled
budget comparison, checkpoint overwrite protection and truncated-particle
payload rejection. It adds changed-inlet restart rejection, requiring
`Restart inlet mismatch` and no result CSV. Rank comparisons retain 3e-12
tolerance; continued parsed gas and particle CSVs must equal the uninterrupted
four-rank result exactly. No tolerance was changed.

**1/1 passed**, test **50.95 s**, total **50.99 s**. Both builds, provenance,
Python syntax, project contract and whitespace passed. No failures or skips;
compiler warnings remain. Logs:
`build/eb-mach-spray-checks/{build-fixed.log,build-selected.log,ctest.log,LastTest.log}`.
Persistent fixture: `build/spray-amr-reactive-debug/tests/spray_mpi_eb_mach/eb_spray_ua4de3ol/`.

The new execution uses the reflection branch, not empirical deposition/splash
or wall-film inlet qualification. Selected execution, full matrix, production
scaling, external solver comparison and physical nozzle/diesel validation were
not performed for this increment.
