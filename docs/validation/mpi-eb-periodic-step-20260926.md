# Periodic MPI EB particle and reacting intervals — 2026-09-26

Source `cdddd3d`, parent `90e01ea`, local main. The mixed-particle interval
accepts explicit periodic/outflow EB faces when an exterior inventory output
is supplied. Faces reach level-set drift, which wraps parcels before existing
EB owner migration. Particle field sampling/deposition use the same periodic
axes. The reacting wrapper maps gas inflow to particle outflow and forwards
periodic axes through both particle halves. Legacy calls retain their existing
boundary path. Explicit EB faces with films or KH child generation reject
until their event handling is connected.

## Completed checks

GNU Fortran 13.3.0, OpenMPI 4.1.6, Python 3.12.3, Ninja Debug; existing
tests-enabled build.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_step --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_(periodic_step|step)_4$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_reactive --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_reactive_ctu-periodic(-subcycle)?$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Particle interval: **2/2 passed**, periodic and legacy four-rank cases each
0.57 s, 1.18 s total. The periodic parcel crosses high x into the opposite
refined patch, migrates from rank zero to rank three, reflects from the EB,
heats and evaporates. Exterior inventory is exactly zero. Gas/parcel/wall
conservation meets the unchanged 3e-12 relative bound. Limited substep count
still preserves the original gas, temperature and particle population.

Reacting interval: **2/2 passed**, common CTU 2.33 s, recursive CTU 4.73 s,
7.11 s total. Both use four ranks and two patches, periodic crossing, nonzero
chemical species exchange, heating/evaporation and EB reflection. The full
gas/parcel/wall/boundary/chemical inventory meets the existing 3e-12 bound.
Invalid gas solver rejects after private particle work without publishing
gas, parcel or diagnostic changes. No failures or skips in either test run.

Provenance, Python syntax, project contract and whitespace checks passed.
Artifacts under `build/mpi-eb-periodic-step-checks/`: `build.log`,
`build-reactive.log`, `step-LastTest.log`, `reactive-LastTest.log`.
Mechanism SHA:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

Application periodic input/restart, periodic films and periodic KH child
generation remain pending. No full matrix, installed audit, external solver
parity or physical validation was run.
