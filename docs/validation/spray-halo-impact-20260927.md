# EB impacts in frozen particle halo stages — 2026-09-27

Implementation `0e6e66a15cc667eac98d8242da133a055d80229e`, based on `7ec7a648`,
local main, version 0.246.0.

## Change

`kick_mpi_amr_mixture_halo` accepts optional `wall_films`, `film_mass` and
`film_settings` together, with `last_id`, nodal geometry and drift controls.
Presence is checked collectively. When supplied, it uses the existing
`advance_mpi_level_set_impacts` flight operator after frozen-gas exchange and
breakup, enabling deposition, splash children and remaining-time flight.

The impact pressure comes from the stage's original gas sample; KHRT shedding
children inherit their parent's sampled pressure using the same ID lookup as
the common-step path. Gas feedback retains its original sampling/deposition
sites. The impact routine's globally summed wall-received inventory is negated
once into the wrapper's change-in-system convention; component outflow entries
are mapped to vapor species without a second MPI sum.

Film payloads and the ID watermark stay private until source deposition and
optional level/physical source export succeed. A failed transfer retains the
input parcels, films, component inventories, source and boundary/wall ledgers.
The optional `model_rejected` diagnostic distinguishes unsupported impact models.

This is the frozen halo kernel connection. Population old/new stages still need
to route physical film state, distinguish physical/copy impacts, and incorporate
new deposits into film clocks. The recursive application remains disabled for
wall films until that lifecycle is connected.

## Independently completed checks

Artifacts: `build/spray-halo-impact-checks/`. Existing fixed-mechanism Ninja
Debug build `build/spray-amr-reactive-debug`; GNU Fortran 13.3.0 and OpenMPI 4.1.6
(`/usr/bin/mpifort`, `/usr/bin/mpiexec`).

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_sample --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_sample_[24]$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final source: **2/2 passed**, 0.89 s on each rank layout, 1.83 s total.
Provenance, build, Python syntax, project contract and whitespace checks passed.
The first build failed because an automatic array declaration preceded its
`components` dummy declaration; moving it after that declaration fixed the error.
The failed build is preserved in `build.log`; final build and test output are
`build-final.log`, `ctest-final.log`, and `LastTest.log`. No test failures or
skips; existing compiler warnings remain.

The planar EB fixture uses 1 m/s and 100 m/s normal impacts. Both produce films;
the faster impact produces multiple live children and advances the ID watermark.
With drag/heating/evaporation disabled, particle plus film mass, normal momentum
and energy balance the wall ledger within `1e-12` of the initial inventory.
There is no outflow. An invalid receive limit forces source routing to reject
after private impact flight, preserving the pre-impact position, component
inventory, empty film output, source, ledgers and ID. Retrying succeeds on both
two and four MPI ranks. Existing particle and film stage checks also pass.

This does not qualify simultaneous evaporation/impact, KHRT shedding plus splash,
outer-boundary impacts, recursive population coupling, restart, selected
mechanisms, physical validation or exact external-solver parity.
