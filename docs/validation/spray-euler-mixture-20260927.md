# Euler mixture parcel update — 2026-09-27

Implementation: `fddddcda`.

`spray_options%mixture_time_scheme='euler'` selects explicit old-state drag
and temperature updates in `advance_mixture_frozen_gas`, matching the update
form in pinned PelePhysics Drag.H. Temperature uses old mass, old mixture heat
capacity and component latent energies. Mass loss retains the existing explicit
update. The gas receives the actual accepted particle momentum/energy change,
preserving the existing discrete conservation contract. This feedback policy
is not asserted identical to upstream's explicit source integration.

The default `'conservative'` scheme is unchanged. This option applies to the
mixed-component free-particle path, not the separate frozen-film parcel API or
wall-film integrator. An Euler drag step exceeding its relaxation time is
rejected, allowing the caller's existing retry/substep mechanism to act.
Invalid temperature updates retain the original particle/component inventory
and publish no gas exchange. Checkpoint identity includes a nondefault scheme;
the default field is omitted from the EB spray identity to retain old snapshots.
The evaporation comparison fixture now explicitly requests Euler stepping.

## Completed checks

GNU Fortran 13.3.0 Debug, OpenMPI 4.1.6, MPI/spray ON, tests/SUNDIALS OFF:
rebuilt `pelef_mpi_reactive_eb_tree_3d_selected_spray` in
`build/heptane-simple-mpi-debug`. Logs: `build/spray-euler-{,final-}build.log`.
The final rebuild includes all changes. Mechanism preflight passed for the
source-pinned fitted heptane bundle from the
[transport integration](simple-transport-integration-20260927.md).

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/heptane-simple-mpi-debug --target pelef_mpi_reactive_eb_tree_3d_selected_spray --parallel 2
python3 tools/prepare_evaporation_reference.py --bundle build/heptane-3sp-simple-reference/heptane_3sp.json --output build/spray-evaporation-euler
# From the prepared directory:
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none /usr/bin/mpiexec -n 1 /home/gtt13/projects/pelec/build/heptane-simple-mpi-debug/pelef_mpi_reactive_eb_tree_3d_selected_spray tree.nml spray.nml
```

One focused unit program passed: `tests/spray/test_spray_mixture_euler.F90`.
It checks an independently calculated 353 K update, remaining/gas mass,
total energy balance, rollback on an invalid temperature and scheme validation.
It was compiled directly (not run through CTest) with:

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH /usr/bin/gfortran -g -fcheck=all -ffpe-trap=invalid,zero,overflow -Ibuild/heptane-simple-mpi-debug/selected_runtime_modules -Ibuild/heptane-simple-mpi-debug/pelef_selected_spray_runtime_modules tests/spray/test_spray_mixture_euler.F90 build/heptane-simple-mpi-debug/libpelef_selected_spray_runtime.a build/heptane-simple-mpi-debug/libpelef_selected_reactive_1d_runtime.a build/heptane-simple-mpi-debug/libpelef_selected_runtime_core.a -o build/spray-evaporation-euler/test_spray_mixture_euler
build/spray-evaporation-euler/test_spray_mixture_euler
```

Four application runs completed: two-step Euler evolution, legacy conservative
snapshot continuation on two ranks, one-step Euler checkpoint, and Euler
continuation on two ranks. Both continuations reproduce the corresponding
continuous gas CSV byte-for-byte and particle rows exactly after combining IDs.
One expected rejection changes the scheme across restart: nonzero exit with
`EB spray context mismatch`, no output CSV. Artifacts:
`build/spray-euler-{old-restart,checkpoint,restart}/` and
`build/spray-evaporation-euler/restart-checks.json`.

## Reference comparison and remaining gap

The unchanged pinned C++ reference from
`build/spray-evaporation-sdc-run2/` is the same periodic 4³, two-step 400 K gas /
310 K liquid adaptation. C++ source identities remain those recorded in the
[full reference diagnostic](evaporation-full-reference-20260927.md).

| Particle value | Fortran minus reference |
| --- | ---: |
| Temperature | `2.6205e-11` K |
| Axial velocity | `-2.4469e-13` m/s |
| Single-drop mass | `-6.5431e-26` kg |
| Axial position | `4.5536e-18` m |

Composite absolute inventory residual is `2.7105e-20`.
This greatly reduces the preceding particle-update differences, but **full
gas/particle parity remains open**. A first all-cell gas comparison finds
relative-to-reference-field-maximum differences of `2.000e-4` for fuel density,
`1.899e-4` for axial momentum density and `1.890e-7` for total energy density.
Gas density differs by `1.130e-9` on the same normalization. Source deposition
and gas feedback need further diagnosis; these metrics alone establish no cause.
The stored plotfile `Temp` field is compared as stored, not recomputed from EOS.

```sh
build/pelec-plotfile-export/pelef_plotfile_csv build/spray-evaporation-sdc-run2/plt00002 build/spray-evaporation-euler/reference-gas.csv
python3 tools/compare_evaporation_particles.py --reference build/spray-evaporation-sdc-run2/spray00002.p3d --candidate build/spray-evaporation-euler/fortran.csv.particles.rank0.csv --reference-gas build/spray-evaporation-euler/reference-gas.csv --candidate-gas build/spray-evaporation-euler/fortran.csv --output build/spray-evaporation-euler/full-comparison.json
```

The diagnostic checks all 64 cells, converts CGS to SI, and records input hashes.
It has no parity acceptance gate; successful execution is not a parity pass.
Python syntax, project contract and diff checks passed. No test matrix,
installed-binary audit or physical validation was run.
