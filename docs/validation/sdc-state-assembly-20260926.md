# SDC state assembly — 2026-09-26

Implementation `70240091`, parent `1e8daa63`, local main, version 0.246.0.

Public `construct_reactive_sdc_state_3d` implements the arithmetic and component
mapping of pinned `Advance.cpp::construct_Snew`. Starting from the old state,
it adds half-dt new then old rates for each source in list order, optional full-dt
hydro rates, and optional full-dt reaction rates for species and total energy.
Reaction layout is species rates followed by the energy rate. Empty source
lists are supported. Shape, interval and finiteness failures preserve the output
state; the candidate is published only after assembly succeeds.

This routine assembles a candidate, not an accepted physical state. As in the
reference separation, EOS recovery and iteration control belong to its caller.
It is not yet connected to an SDC time-step driver.

## Verification

`unit_reactive_ppm_step_3d` passed **1/1**, 5.06 s test / 5.11 s total. Added
exact binary-arithmetic fixtures check two source families plus hydro/reaction,
an empty-list identity, reaction exclusion from density/momentum, and unchanged
outputs for negative dt or inconsistent source lists. Existing CTU, source and
split checks also ran. No failures/skips. Provenance, Python syntax, project
contract and whitespace passed. GNU Fortran 13.3.0, Python 3.12.3, existing
tests-enabled Ninja Debug build; existing compiler warnings remain.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_reactive_ppm_step_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_reactive_ppm_step_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/sdc-state-assembly-checks/{provenance.log,build.log,ctest.log,unit.log}`.
Reference inspected with `git show` at
`bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8` in
`/mnt/sn850/Projects/gtt2/external/PeleC`.
This is a focused source-grounded check. SDC iteration, coupled chemistry,
application dispatch, external executable parity and physical validation remain
open.
