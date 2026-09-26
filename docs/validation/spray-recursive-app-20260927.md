# MPI EB recursive spray application — 2026-09-27

Implementation `34979b0fc4ab4eb2741c252a33bf95ec05ed6917`, based on `c714bcd7`,
local main, version 0.246.0.

## Behavior

`pelef_mpi_eb_spray_3d` now dispatches to the conservative recursive library
when `/eb_tree/` selects `time_integrator='sdc'` and `subcycling=.true.`.
Existing CTU/PeleC and explicit-face requirements remain. Chemistry, molecular
transport, gas force density, inflow/isothermal metadata, source history and
accepted step ledgers are forwarded. The computational child-creation count is
zero for this core. Normal outer-step scheduling, retry and output remain in
the application.

Recursive breakup, wall films/film impacts, nonzero particle acceleration and
periodic migration reject at input. These are unfinished features, not a claim
that the full port is complete. The checkpoint context adds
`RECURSIVE_SDC_PHYSICAL_SYNC_V1` for recursive SDC so common-step and recursive
particle histories cannot share an identical context. Changed-rank recursive
continuation and regrid were not exercised in this check.

The [physical-source synchronization](spray-recursive-source-sync-20260927.md)
is a conservative PeleF extension to the pinned transfer algorithm. Exposing
it in an application does not establish external solver parity or physical
spray validation.

## Completed checks

Artifacts: `build/spray-recursive-app-checks/`; final generated cases:
`build/spray-amr-reactive-debug/spray-sdc-recursive-app/eb_spray_4vxn8hrj/`.
Existing Ninja Debug fixed-mechanism build, GNU Fortran 13.3.0,
`/usr/bin/mpifort` and `/usr/bin/mpiexec` OpenMPI 4.1.6.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_eb_spray_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_sdc_recursive_app$' --output-on-failure
python3 -m py_compile tests/application/check_mpi_eb_spray.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final CTest **1/1 passed**, 105.80 s (105.85 s total), no skips. Builds,
provenance, Python syntax, contract and whitespace checks passed. Existing
compiler/linker warnings remain.

The checker runs the planar EB water-parcel fixture for three outer steps to
30 ns on 2 and 4 ranks, with chemistry, all molecular transport and gas force
density. It compares composite gas CSVs and ID-sorted physical particle CSVs
with the existing 3e-12 rank-comparison criterion. It checks wall reflection,
heating, positive masses and component sums. Both runs process one injection
event. Nonzero particle acceleration and film-impact input are rejected before
result output; a missing restart manifest is rejected separately. This is five
MPI launches within one CTest, not five independently qualified solver modes.

The new budget check uses 3e-12 times a component reference scale. Mass, energy
and species use the maximum initial/final/injected magnitude. Momentum also
uses sqrt(mass_scale * energy_scale), with momentum units, because this fixture
starts with zero net momentum. The maximum scaled residual is
2.026379654940504e-15 on both ranks layouts. The application reports its separate
existing budget metric as 1.3665261161999555e-17.

The initial checker failed after both physical runs and rank comparison passed
(104.46 s). Its newly added budget assertion divided a momentum residual
-2.0964514981663939e-24 by the tiny final net momentum 1.7318055652586352e-13.
The checker was corrected to the stated characteristic momentum normalization;
this changes that new assertion's effective bound. No solver tolerance or
existing library conservation criterion was changed. The initial failure is
retained in `ctest.log` and `LastTest-first.log`; the final run is in
`ctest-final.log` and `LastTest-final.log`.

Only the fixed-mechanism executable was built and run here. Selected-mechanism
execution, recursive checkpoint continuation, scheduled regrid, three-level
coupled evolution, temporal order and physical validation remain separate work.
