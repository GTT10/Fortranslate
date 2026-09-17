# Coupled MPI AMR gas/particle remapping — 2026-09-18

Implementation: `4c0d4899e4ae23cb4d3ce67372aed969e7a47404`, based on
`457d2f3`. `remap_mpi_amr_spray` validates original gas thermodynamics,
migrates private parcel/component arrays, then conservatively remaps gas and
recovers temperatures. Tree, owner map, gas, temperature and particle payloads
publish together only after success. The physical root domain and root
resolution remain fixed; topology and owners may change. Periodic flags govern
particle lookup, and nonperiodic exterior particles reject. The caller selects
the target tree/owners and receive capacity and retains its unchanged ID watermark.

## Completed focused verification

Reused the preceding clean GNU Fortran 13.3.0 Debug/Ninja builds, with explicit
`/usr/bin/gfortran`, Open MPI `/usr/bin/mpifort` and `/usr/bin/mpiexec`. Compiler,
MPI and backend configuration did not change. Fixed tests were enabled and
selected H2/O2 tests disabled; no installed prefix was used. Configuration
commands are in [the preceding evidence](spray-amr-reactive-20260918.md).

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_amr_remap --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_amr_remap_(2|4)$' --output-on-failure --output-log build/spray-amr-remap-checks/ctest.log
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final **2/2 passed**, no skips: two ranks **0.52 s**, four ranks **0.54 s**,
**1.10 s** total. A nonuniform three-level source is coarsened to the root,
refined with moved patches and returned to the original geometry, using changed
owners. Every owned conserved field and recovered temperature matches serial
remapping exactly. Composite gas inventories remain within 3e-13 of each initial
inventory scale. Particle count, phase-space, component mass, multiplicity,
material/lineage and deformation/KH-RT payloads are preserved; all particles
reach their target finest owners without nonowner gas storage.

Receive-capacity failure and invalid original gas preserve all published state.
A valid target with incompatible root resolution permits private particle
migration but rejects gas remapping, verifying late failure rollback of topology,
owners, gas, temperature and particle outputs.

The first run failed both topology-success checks because the fixture assigned
parent IDs greater than child IDs. Correcting the fixture to valid lineage and
adding an explicit initial-parcel validity check produced the final results;
no solver tolerance or validation rule changed. Initial Open MPI output also
contained display-authorization messages, but processes ran and failed at the
explicit fixture assertions.

The selected library compiled successfully; this is compile evidence, not
selected executable qualification. Provenance, Python syntax, project contract
and whitespace checks passed. No Python source changed. No executable regrid
scheduler, restart, external solver parity, convergence, scaling or physical
spray validation was exercised.

Logs: `build/spray-amr-remap-checks/{build,build-selected,ctest-initial,ctest}.log`.
Executable: `build/spray-amr-reactive-debug/tests/test_spray_mpi_amr_remap`.
