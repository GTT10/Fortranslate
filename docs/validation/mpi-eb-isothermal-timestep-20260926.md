# MPI tree thermal timestep — 2026-09-26

Implementation f944ec0, parent 4a91ca7, local main, version 0.246.0.
MPI timestep selection now accepts wall temperatures, collectively agrees
presence, values and isothermal activity, and computes the existing heat-transfer
bound on owned domain-wall patches. Recursive transport scaling matches serial;
a final global minimum chooses the root interval. Invalid/mismatched controls
reject before publishing dt. Empty ranks participate in the same collectives.

## Completed focused check

GNU Fortran 13.3.0, Open MPI 4.1.6, Python 3.12.3, Ninja Debug. Existing fixed
tests-enabled build/spray-amr-reactive-debug, no backend changes.
The three-rank isothermal transport fixture passed **1/1**, 4.61 s (CTest total
4.66 s), no failures or skips. It now tests common and recursive thermal bounds
with distributed and all-on-rank-zero ownership, including empty ranks.
Serial/MPI dt values match exactly. With 2000 K walls and conduction-only
comparison, recursive baseline is 5.6633521555547262e-11 s and the thermal limit
is 1.7812646913305850e-11 s. Common-step dt remains controlled by its existing
3.5395950972217039e-12 s bound. Rank-zero temperature 2001 K rejects and preserves
the sentinel dt. Existing transport state/ledger parity and rollback checks
remain active in the same fixture.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_eb_tree_transport_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_mpi_eb_tree_isothermal_transport_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Provenance, syntax, project contract and whitespace passed. Compiler warnings
remain. No tolerance or golden-data changes.
Artifacts: build/mpi-eb-isothermal-timestep-checks/build.log and LastTest.log.
Mechanism SHA: 0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4.
Pinned PeleC: bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8.
Pinned PelePhysics: 6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b.

## Limits

Full reacting wall-temperature forwarding and application input/restart remain
to be connected. No full matrix, installed audit, external solver parity or
physical validation ran. Complete PeleC parity remains unproven.
