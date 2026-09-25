# Periodic MPI EB spray application — 2026-09-26

Source `3200545`, parent `9a1978d`, local main. The MPI EB spray application
accepts periodic tree faces through the connected gas/particle library path.
Explicit boundaries still require `particle_outflow=.true.`. Periodic cases
with wall impacts/films or KH child generation reject at input. Existing
boundary checkpoint identity rejects changed periodic labels on restart.

## Completed checks

GNU Fortran 13.3.0, OpenMPI 4.1.6, Python 3.12.3, Ninja Debug; existing fixed
tests-enabled build.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_eb_spray_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none python3 tests/application/check_mpi_eb_spray.py --exe build/spray-amr-reactive-debug/pelef_mpi_eb_spray_3d --mpiexec /usr/bin/mpiexec --source . --work build/mpi-eb-periodic-application-checks/application --periodic --restart
python3 -m py_compile tests/application/check_mpi_eb_spray.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Checker **PASS**: two-/four-rank continuous runs, a two-rank checkpoint stop,
and four-rank continuation. The modified nozzle starts just below the high-z
periodic face with positive z velocity; final z lies near the low face.
Parcels also reflect from the y-normal EB and change temperature. Gas and
particle data agree across rank counts within the existing 3e-12 bound;
continued gas and particle rows equal the uninterrupted four-rank rows exactly.
Restart budget fields satisfy the existing 3e-12 bound. Particle component
masses sum to parcel mass within 1e-14 relative tolerance.

Expected failures cover missing restart, attempted checkpoint overwrite,
changed spray context, changed hydro method, changed periodic labels and
truncated particle payload. Checkpoint overwrite leaves all files unchanged.
No unexpected failures or skips. Provenance, Python syntax, project contract
and whitespace checks passed.

Artifacts: `build/mpi-eb-periodic-application-checks/build.log` and
`application/eb_spray_6dltsnf9/` below that directory.
Mechanism SHA:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

Selected-mechanism periodic spray, periodic regridding/rebalancing, films and
KH child generation were not qualified here. No full matrix, installed audit,
external solver parity or physical validation was run.
