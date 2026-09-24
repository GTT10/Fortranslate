# PeleMP cap input in EB tree frontends — 2026-09-25

Implementation `19a71c6c130da03091913ea7204bb04e8b100a75`, base `aebdb94`.

The preceding integration added the normal boiling temperature to the regular,
MPI AMR and single-level EB drivers but missed the separate EB tree frontend.
Its shared fixed/selected serial/MPI source now accepts
`mixture_normal_boiling_temperatures`, rejects unused entries, assigns each
liquid record and includes the field in its serialized namelist context.
No property equation changed in this correction.

GNU Fortran 13.3.0, Open MPI 4.1.6, Debug/Ninja, Python 3.12.3; unchanged fixed
tests-enabled and selected H2/O2 tests-disabled configurations were reused.
Fixed MPI EB and selected MPI EB applications built. Selected was not executed.

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_eb_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_reactive_eb_tree_3d_selected_spray --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_pelemp_cap$' --output-on-failure --output-log build/pelemp-eb-input-checks/ctest.log
python3 -m py_compile tests/application/check_mpi_eb_spray.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

CTest **1/1 passed**, **47.07 s** test / **47.11 s** total. The single-liquid
fixture uses Tref=300 K, normal Tb=330 K, Tc=647 K and initial droplets at 350 K,
so the cap is active. It exercises two/four-rank gas/particle agreement at the
existing relative 3e-12 bound, wall reflection, two-rank checkpoint/four-rank
continuation with exactly equal parsed gas/particle CSVs, inventory/ledger
agreement, missing/truncated checkpoint and overwrite rejection. The context
mismatch case changes normal Tb from 330 to 331 K and is rejected.

No failures or skips in this correction. Mechanism provenance passed earlier
in the same work session; Python syntax, project contract and whitespace passed.
Compiler warnings remain. Artifacts: `build/pelemp-eb-input-checks/` contains
fixed/selected build logs, `ctest.log` and copied `LastTest.log`; fixture paths
are recorded in those logs under
`build/spray-amr-reactive-debug/tests/spray_mpi_eb_pelemp_cap/`.
No full matrix, new selected execution, regrid/scaling test, external solver
parity, physical boiling validation or shared-vapor MPI execution was performed.
