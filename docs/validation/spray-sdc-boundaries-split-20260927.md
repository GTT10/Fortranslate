# MPI EB spray SDC boundaries and splitting — 2026-09-27

Implementation `61bdba24`, based on `3f050e4c` on local `main`.
PeleC reference `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`;
PelePhysics reference `6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`.
Mechanism SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

## Implementation

`advance_mpi_eb_spray_sdc` now forwards optional `inflow_state`,
`inflow_mach`, `inflow_temperature`, `inflow_pressure` and `wall_temperature`
to the existing gas SDC driver. Particle inflow crossings retain outflow
behavior, and isothermal walls retain specular particle reflection; wall
thermal data applies to the gas operator.

After all SDC iterations succeed, computational multiplicity splitting runs
once on the private final parcel population. It uses the existing global
parent-ID ordering and component-mass handling. Optional `split_created`
returns the number of new IDs, or zero on failure. Capacity failure leaves
gas, temperature, reaction history, particles, ID and budgets unpublished.
The maximum multiplicity must agree across ranks, including enabled versus
disabled splitting; invalid spray options reject collectively before work.

This extends the existing common-step library API. Application dispatch,
spray SDC restart, film integration, particle acceleration and recursive
particle scheduling remain open. No upstream scheduler parity is claimed.

## Focused evidence

GNU Fortran 13.3.0, Ninja Debug, fixed H2O2 tests-enabled build
`build/spray-amr-reactive-debug`, OpenMPI 4.1.6 through `/usr/bin/mpifort`
and `/usr/bin/mpiexec`. Existing comparison and executable-stack warnings remain.

Final CTest **3/3 passed**, 18.07 s:

| Test | Seconds | Combined relative inventory residual |
| --- | ---: | ---: |
| `spray_mpi_eb_sdc_2` | 5.92 | 1.1968665763966046e-15 |
| `spray_mpi_eb_sdc_inflow_2` | 6.12 | 1.1968665763966046e-15 |
| `spray_mpi_eb_sdc_isothermal_2` | 5.99 | 1.7767815781073639e-15 |

All use two ranks, two AMR levels, two SDC iterations, chemistry and all
molecular transport components, mixed-liquid evaporation/heating/drag,
EB reflection and final computational splitting. A multiplicity-100 parcel
becomes two parcels with one new ID. The combined gas/parcel/boundary/wall/
chemical budget retains the existing 3e-12 relative tolerance. The isothermal
case additionally exercises an outer reflecting wall and 2000 K gas wall data.
The inflow case rejects an invalid rank-local inlet density without publication.
Mach data forwarding is implemented but not independently exercised here.

Each mode rejects capacity one after the coupled candidate has been computed,
retaining original gas, temperature, component masses, population count, ID,
source and boundary budgets. Enabling splitting on only one rank is also
rejected without publication. Existing invalid-hydro rollback checks remain.

The initial three tests failed because the new test fixture enabled splitting
while leaving `breakup_model='none'`, which existing option validation forbids.
The fixture now explicitly uses ETAB with high surface tension, following the
existing computational-splitting test setup; no production guard or tolerance
was loosened. After this correction all three passed, and passed again after
adding collective multiplicity agreement. No build failures occurred.
Provenance, Python syntax, project contract and `git diff --check` passed.
No full matrix, selected-mechanism build, installed audit, physical validation
or upstream trajectory comparison was run.

## Reproduction

Artifacts: `build/spray-sdc-boundary-checks/`. Initial failure is retained as
`ctest-initial-failure.log`; final evidence is `ctest-final.log`,
`detail-final.log`, `build-agreement.log`, `provenance.log` and `contract.log`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_reactive --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_sdc(_inflow|_isothermal)?_2$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
