# MPI EB impacts after domain reflection — 2026-09-26

Implementation `58814385`, parent `ee12947d`, local main, version 0.246.0.
The MPI EB impact interval accepts collectively agreed nonperiodic outer faces.
It clips each flight to the earliest outer face, resolves any earlier EB impact,
then removes outflow or reflects at slip/no-slip faces and continues the remaining
flight. Domain impulse is retained in the exterior ledger; EB impact/film exchange
stays in its wall ledger. The event-depth cap covers both event types and all
candidate state/IDs/films publish transactionally. EB film options are validated
against the embedded surface separately from the exterior labels.

The particle step and application now permit EB films with explicit nonperiodic
domain walls. Outer walls remain specular; this adds no outer-domain film model.
Periodic films still reject. Isothermal gas labels use the existing particle
no-slip mapping while EB film temperatures remain separate inputs.

## Focused evidence

GNU Fortran 13.3.0, Open MPI 4.1.6, Python 3.12.3, Ninja Debug, fixed H2/O2,
existing tests-enabled `build/spray-amr-reactive-debug`.
The direct application checker passed **6 successful runs and 2 expected context
rejections**: deposition and splash, each with uninterrupted MPI4, MPI2 checkpoint
and MPI4 continuation. The CTU fixture reflects on the upper x domain wall before
reaching the planar EB, includes three regrids and measured-rebalance controls.
Each case creates two EB films; deposition consumes the parents, while splash
retains eight children. Gas, particle, film and budget continuation meet the
unchanged 3e-12 comparison. Film temperatures, components and restart IDs retain
the existing checks. Additional saved-output assertions verify both films lie
on the reflected side of the injection x position in full and resumed outputs.
Maximum absolute budget residual is 1.452831e-17.

The initial build failed on a new overlength Fortran signature; it was wrapped
and the final build passed. No runtime failures or skips occurred. A CTest entry
was registered; the checker was run directly, not counted twice. The reflection
position assertions were added while the checker ran and then evaluated against
all four saved full/resumed outputs in `reflection-check.log`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_eb_spray_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none python3 tests/application/check_mpi_eb_impact.py --exe build/spray-amr-reactive-debug/pelef_mpi_eb_spray_3d --mpiexec /usr/bin/mpiexec --source . --work build/mpi-eb-film-domain-checks/application --reflect-domain
python3 -m py_compile tests/application/check_mpi_eb_impact.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Provenance, Python syntax, contract and whitespace passed; tolerances and golden
data unchanged. Build warnings remain. Artifacts under
`build/mpi-eb-film-domain-checks`: `build-initial.log`, `build.log`,
`application.log`, `reflection-check.log`, `application/impact__ap_lwmj`.
Mechanism SHA: `0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
PeleC: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
PelePhysics: `6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`.

No selected mechanism, full matrix, external solver comparison, scaling or
physical impact/film validation ran. Complete PeleC parity remains unproven.
