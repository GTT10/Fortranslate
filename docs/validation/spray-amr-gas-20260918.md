# Transactional AMR gas feedback — 2026-09-18

Implementation: `28b059fc274f9ab2f947c9a53d9e72050ec99864`, based on
`4ffc9e8`. `deposit_amr_tree_spray_gas` wraps the existing composite field
depositor with initial-gas validation and post-update EOS temperature recovery
on every patch. It retains private copies until every cell, including averaged
parents, accepts density, composition and energy. U/T publish together;
failure preserves both original fields. Empty feedback validates inputs and
preserves state and temperature.

The wrapper does not change conserved values during temperature recovery. The
existing cloud-overlap stencil, volume scaling, covered-cell exclusion and
bottom-up averaging remain unchanged. This is serial gas feedback; owner-local
MPI transfer and coupled parcel stepping remain separate work.

## Completed focused verification

GNU Fortran 13.3.0 Debug/Ninja, `/usr/bin/gfortran`. Reused compatible fixed and
selected H2/O2 serial build directories; compiler/backend settings unchanged.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-partition-debug --target test_spray_amr_tree_gas --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-partition-debug -R '^spray_amr_tree_gas$' --output-on-failure --output-log build/spray-amr-gas-checks/ctest.log
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-partition-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**1/1 passed**, no failures/skips, test time **0.08 s**, total **0.11 s**.
The three-level fixture deposits moving 1500 K water-vapor contributions into
1000 K nitrogen across coarse/fine, nested and periodic interfaces. Conserved
fields exactly match the generic field depositor; refreshed parent values
exactly match another bottom-up averaging pass. Every cell temperature matches
a separate conserved-to-primitive call. Composite mass, momentum, energy and
species changes equal supplied contributions within 3e-13 of the component's
initial-inventory/contribution scale.

Negative-density, out-of-range-energy and negative-species feedback each reject
and preserve every original U/T cell. The test also checks empty feedback and
rejection of an invalid initial gas state. This is an admissibility and
conservation test with prescribed feedback, not an evaporation-rate or reacting
spray benchmark.

The selected library compiled successfully. Provenance, Python syntax, project
contract and whitespace checks passed. No failed build/test attempt. No Python
source changed. No selected executable, MPI feedback, full matrix, installation,
AMR convergence or physical validation is claimed.

Logs: `build/spray-amr-gas-checks/{build,ctest,build-selected}.log`.
Executable: `build/spray-amr-partition-debug/tests/test_spray_amr_tree_gas`.
