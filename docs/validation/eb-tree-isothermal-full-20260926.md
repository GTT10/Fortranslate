# Reacting EB tree isothermal walls — 2026-09-26

Implementation `09eddf5`, parent `a703e97`, local main, version 0.246.0.
Serial/MPI full R/2 T/2 H T/2 R/2 APIs now accept optional wall temperatures
and forward them through both transport half-steps, common or recursive.
Isothermal faces require temperature, enabled transport and conduction.
Supplied temperatures must be finite and positive and require face labels.
MPI validation is collective; transport agrees temperature values before use.
Candidate state and ledgers remain private until the complete split succeeds.

## Completed focused verification

GNU Fortran 13.3.0, Open MPI 4.1.6, Python 3.12.3, Ninja Debug;
existing fixed-mechanism tests-enabled `build/spray-amr-reactive-debug`.
No compiler/backend configuration changed.

Two CTest cases passed, with no failures or skips: recursive 108.62 s,
common 85.54 s, total 194.21 s. The five-patch, three-level planar fixture
uses three ranks, 2000 K x/y walls, z outflow, PLM/CTU and explicit/implicit
chemistry. Explicit chemistry uses distributed owners; implicit chemistry puts
all patches on the last rank, leaving two ranks empty.

Serial/MPI conserved state, temperature, boundary ledger and limiter match
exactly. Wall and chemistry ledgers meet the unchanged 3e-14 check.
Composite normalized conservation residual is at most 4.094700e-15.
A rank-zero wall temperature of 2001 K rejects after private chemistry and
preserves state, temperature and all public ledger/limiter sentinels.
Invalid hydro also rejects after private reaction/transport without publishing
partial state. The existing timestep comparison runs with wall temperatures.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_eb_tree_subcycle_full_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_mpi_eb_tree_ctu-isothermal(-common)?_full_3d_3$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Provenance, syntax, contract and whitespace checks passed. No tolerances or
golden values changed. Compiler warnings and runtime X authorization messages
remain in logs; neither test failed.
Artifacts: `build/mpi-eb-isothermal-full-checks/{build.log,ctest.log,LastTest.log}`.
Mechanism SHA: `0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
PeleC: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
PelePhysics: `6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`.

## Remaining scope

Tree application input, restart temperature binding and spray forwarding remain
open. No full matrix, installed audit, external solver comparison or physical
validation ran. These focused checks do not establish complete PeleC parity.
