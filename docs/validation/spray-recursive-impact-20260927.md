# Recursive impact-to-film thermal handoff — 2026-09-27

Integration `5e14911a01e1b0ec6e28e2390546e46d0dc61acc`, following reflection fix
`00ac85da3b861807f686a25ee809b4d63c290254`, based on `9cf1dca1`, local main,
version 0.246.0. Numerical reflection and integration changes are separate commits.

## Change

The recursive particle coordinator forwards optional film inventories and impact
settings to its SDC stages. `advance_mpi_eb_spray_sdc_subcycle` accepts
`film_impacts` and `model_rejected`; enabled impacts require the complete film
argument group on every MPI rank.

At the old callback, particles see the existing physical film inventory for
impact thickness. After both particle and film old stages succeed, only the
appended impact-film suffix is attached to the physical thermal midpoint.
Existing films retain their old half-kick result and newborns receive only the
final thermal half kick. Particle-side snapshots are never assigned over the
evolving film hierarchy. Intermediate callbacks reuse cached sources; final
callbacks combine particle and film rates and physical wall inventories.

All state remains private through recursive gas synchronization and final
collection. Failed capacity/model/stage checks retain caller gas, temperature,
particles, films, component masses, IDs and inventories. The API supports this
coupled path; application input routing remains disabled in this increment.

## Reflection failure found during integration

The initial impact test failed in an old root virtual-particle reflection,
before physical splash evolution. Re-interpolating the contact at
`y=0.00125 m` produced a level-set value of `-4.7704895589343460e-21`, causing
the remaining reflected segment to reject as starting in solid.

A new dimensional-plane regression independently reproduced that rejection.
The fix moves each contact coordinate with nonzero fluid-normal component to
the adjacent representable value in that component's direction before continuing
the reflected segment. It changes neither velocity nor inventory and does not
relax the intersection acceptance criterion. Existing oblique, cavity, outflow,
periodic and rollback checks in the drift regression also pass.

## Independently completed checks

Artifacts: `build/spray-recursive-impact-checks/`. Existing fixed-mechanism Ninja
Debug build `build/spray-amr-reactive-debug`; GNU Fortran 13.3.0 and OpenMPI 4.1.6
(`/usr/bin/mpifort`, `/usr/bin/mpiexec`).

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_reactive --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_sdc_recursive_(impact|film)_2$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_eb_drift test_spray_mpi_eb_reactive --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^(spray_eb_drift|spray_mpi_eb_sdc_recursive_impact_2)$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final source: **2/2 passed**, drift 0.00 s and two-rank coupled impact 8.36 s,
8.41 s total (`ctest-drift-fix.log`, preserved `LastTest.log`). Provenance,
build, Python syntax, project contract and whitespace checks passed. No build
failures or skips; existing warnings remain. The earlier pre-fix run passed
stationary-film coupling in 7.30 s and failed impact coupling in 0.81 s. Three
diagnostic reruns isolated the same failure; the new standalone drift regression
also failed before the fix. Logs are retained as `ctest.log`, `ctest-debug*.log`
and `ctest-drift-repro.log`; temporary diagnostics were removed.

The coupled fixture starts with two mixed films and one impacting mixed parcel,
with gas chemistry, molecular transport, particle transfer and film evolution
enabled. Insufficient parcel capacity rejects transactionally, followed by a
successful retry. The accepted interval has **four live splash particles and
three films**. Combined gas/particle/film/wall inventory has maximum normalized
residual **1.5972576137871884e-15**, below the unchanged `3e-12` criterion.
One root and two child intervals produce six particle and six film stages.

No application run, four-rank coupled impact, restart, cross-level impact,
selected mechanism, physical validation or full external-solver parity was run.
