# Regular MPI spray PLM viscosity — 2026-09-25

## Source and scope

Code `e7fafae`, parent `796e2e4`. Regular MPI spray PLM hydro applies explicit
positive artificial-viscosity coefficients at both SSPRK2 stages, before cell
updates and boundary accounting. Periodic calls with a positive coefficient
use the existing full boundary ghost fill so nodal divergence has transverse
ghost cells. Physical calls use the same path. No-slip masks are restricted to
physical domain ends along the decomposed axis. Omitted/zero PLM policy remains
unchanged. PPM dispatch remains unchanged.

The serial spray wrapper now accepts positive PCM/PLM coefficients, and the
MPI mixture wrapper accepts positive PLM coefficients, forwarding them through
their existing hydro calls. Those coupled paths were compiled but not separately
run here. Regular MPI spray currently supports PLM/PPM, not PCM. Non-PPM
application input/default/restart policy remains outstanding.

## Completed checks

GNU Fortran 13.3.0, OpenMPI 4.1.6, Ninja Debug, Python 3.12.3; reused
`build/spray-amr-reactive-debug`. Focused CTest **1/1 passed**, 11.49 s
(11.53 s elapsed). No failures or skips; compiler warnings remain.

The existing four-rank, two-cell-slab fixture gained a PLM selector. Its 24
combinations cover HLLC coefficient zero and PeleC coefficients 0.1/0.037,
periodic boundaries, slip-wall/outflow and six inflow orientations. Conserved
states and temperatures match serial exactly; surface budgets meet the existing
2e-13 scaled bound. Large-step, selector mismatch, coefficient mismatch and
negative-coefficient checks retain state. This fixture does not separately test
no-slip viscosity or moving parcels. It is library parity evidence, not external
PeleC parity, a full configuration matrix or physical validation.

Mechanism provenance, Python syntax, project contract and whitespace checks
passed. No tolerance changes.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_ppm --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_plm_viscosity_4$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/plm-mpi-spray-viscosity-checks/build.log` and copied
`LastTest.log`.
