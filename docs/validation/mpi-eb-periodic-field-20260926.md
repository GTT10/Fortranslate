# MPI periodic EB particle field exchange — 2026-09-26

Source `57bb1b9`, parent `185bb51`, local main. MPI EB field sampling and
deposition forward periodic axes to the composite EB stencil. Collective
geometry validation checks root face metrics and nodal endpoint agreement
even when all ranks have zero particle queries. Existing owner-routed cell
requests, conservative feedback, parent restriction and EOS recovery are reused.

## Completed checks

GNU Fortran 13.3.0, OpenMPI 4.1.6, Python 3.12.3, Ninja Debug; existing
tests-enabled build.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_sample test_spray_mpi_eb_deposit --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_periodic_(sample|deposit)_4$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**2/2 passed**, sample 0.50 s, deposit 0.52 s, total 1.06 s; no failures or
skips. Four ranks own two patches on ranks zero and three. A query at the
x/z periodic corner reaches the opposite refined patch; sampling agrees with
serial stencil weights to the existing 1e-13 bound. Deposition preserves the
analytical composite inventory and fluid-weighted parent restriction to the
existing 1e-13 bound. Covered gas and temperature slots remain unchanged.
Nonphysical feedback preserves gas and temperature outputs. Sampling retains
capacity, rank-geometry mismatch and solid-query rollback checks; globally
empty queries still reject mismatched periodic nodal endpoints.

Two-rank test entries were added alongside four-rank entries; only the
four-rank pair above was run. Provenance, Python syntax, project contract and
whitespace checks passed.
Artifacts: `build/mpi-eb-periodic-field-checks/build.log` and `LastTest.log`.
Mechanism SHA:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

MPI particle stepping and coupled periodic spray application input remain
pending. No full matrix, installed audit, external solver parity or physical
validation was run.
