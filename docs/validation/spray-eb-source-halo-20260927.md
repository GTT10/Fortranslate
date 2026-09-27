# Inclined-wall SDC source donor repair — 2026-09-27

Implementation: `d1670517`, based on `0a8583f1`. The pinned Spray-EB
8³ connection failed before hydro flux evaluation. GDB traced the first
rejection to `fill_sdc_mpi_eb_tree_source_halo_3d`: target index `[8,-2,-2]`
was mapped to clamped donor `[8,1,6]`, whose fluid volume was zero. The state
halo and extended geometry already reflected physical walls, so that target
required a fluid donor at the mirrored interior position.

Serial and MPI SDC source halos now mirror donor positions at slip, no-slip
and isothermal walls, including repeated reflections for wide halos. Periodic
mapping and other physical-boundary clamping are retained. Source components
keep the existing scalar extension convention; this change introduces no
momentum sign flip, EOS recovery, tolerance change or covered-donor fallback.
The MPI same-level-only path still leaves physical exterior values zero.

## Completed verification

GNU Fortran 13.3.0 / OpenMPI 4.1.6, existing Debug MPI/spray build, selected
pinned reaction-free heptane Simple bundle; tests and SUNDIALS disabled.

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/heptane-simple-mpi-debug --target pelef_mpi_reactive_eb_tree_3d_selected_spray pelef_mpi_reactive_eb_tree_3d_selected pelef_mpi_support --parallel 2 > build/spray-eb-halo-fix-build.log 2>&1
env -u LD_LIBRARY_PATH -u LIBRARY_PATH /usr/bin/mpifort -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -ffree-line-length-none -I build/heptane-simple-mpi-debug/mpi_modules -I build/heptane-simple-mpi-debug/modules tests/unit/test_mpi_eb_tree_hydro_3d.F90 build/heptane-simple-mpi-debug/libpelef_mpi_support.a build/heptane-simple-mpi-debug/libpelef_core.a -o build/spray-eb-halo-unit/check
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none /usr/bin/mpiexec -n 2 build/spray-eb-halo-unit/check ctu-wall > build/spray-eb-halo-unit/run.log 2>&1
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

All listed commands passed. The focused two-rank wall test includes serial/MPI
temporal source comparison, failure rollback and CTU hydro. MPI inventory
residual is `4.6333897726434849e-15` in both owner layouts.

The unchanged failed connection inputs were copied into
`build/spray-eb-halo-fixed` and run there with:

```sh
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none /usr/bin/mpiexec -n 1 /home/gtt13/projects/pelec/build/heptane-simple-mpi-debug/pelef_mpi_reactive_eb_tree_3d_selected_spray tree.nml spray.nml > fortran.log 2>&1
```

The selected spray target had linked before this launch; the remaining build
targets subsequently completed. The executable SHA-256 is
`8baf0d635514dad5fd2f446fe0c2b31b4092a28418401321c16758ab95081c82`.
Two steps reach `2e-7 s` with all five parcels. Injected mass agrees with the
source-derived physical droplets: `9.325346553283266e-11 kg`. Total mass
residual is `1.8048394664492586e-19 kg`; energy residual is
`-6.879929568756604e-14 J` against an initial `2053.265037516425 J`.
Artifacts include `fortran.log`, gas/particle/budget CSVs and `result.json`.
GDB traces and original failed runs remain in `build/spray-eb-debug`,
`build/spray-eb-connection` and `build/spray-eb-gas-diagnostic`.

One diagnostic launch initially used the repository working directory and
failed to find the relative geometry; its log is preserved as
`build/spray-eb-gas-diagnostic/wrong-working-directory.log`. It is not counted
as a solver check. Existing array-temporary, X11 and executable-stack warnings
did not prevent the successful runs. No Python implementation changed, full
matrix, full 64³ duration or external-solver parity was run in this increment.

## Separate full-jet status

The earlier decane jet process ended with exit 1 after `4368.8009266946465 s`,
at `2026-09-27T09:18:53.434278+00:00`. Its error is `Tree step retry limit
exhausted`. The completion-gated inventory checker rejects this incomplete
run. `build/decane-jet-full-run/execution.json`, logs and checkpoints remain
intact; no restart was launched. That executable predates this change and
uses all-fluid geometry, so this wall-donor diagnosis does not establish the
cause of its failure.
