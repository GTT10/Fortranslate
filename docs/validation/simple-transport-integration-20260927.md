# Simple transport import and spray integration — 2026-09-27

Import/generation: `cc8ec1ef`; numerical/application integration: `4d0f2890`.

The optional heptane importer now verifies the exact pinned PelePhysics
`mechanism.H` SHA and Git blob before extracting viscosity, conductivity and
triangular binary-diffusion fits. The bundle retains those source bytes and
identities. Provenance checks reconstruct the coefficients; the Fortran
generator rejects changed coefficients or header bytes. CMake tracks the header
and importer as generation dependencies. Bundles without fits retain identical
generated source.

Selected serial/MPI EB tree applications accept
`simple_transport_enabled=.true.` only when built with a fitted bundle.
The transport records carry their own coefficients; gas transport and particle
film transport both use the [ported Simple kernel](simple-transport-kernel-20260927.md).
The spray path also applies pinned Drag.H's `MW_skin / MW_fuel` conversion to
the species diffusion coefficient. The existing analytical route remains the
default. This is currently an opt-in EB tree application interface, not a new
input option for every executable.

Checkpoint records and MPI configuration comparisons include the coefficients.
Coefficient-free checkpoint serialization retains its original layout.
Changing transport selection across restart is rejected.

## Completed checks

Fresh build: `build/heptane-simple-mpi-debug`, GNU Fortran 13.3.0 Debug,
OpenMPI 4.1.6, MPI/spray ON, tests/SUNDIALS OFF. Configured with explicit
`/usr/bin/gfortran`, `/usr/bin/mpifort`, `/usr/bin/mpiexec`, and
`PELEF_MECHANISM_BUNDLE`/`PELEF_MECHANISM_SOURCE` pointing to the new
`build/heptane-3sp-simple-reference/heptane_3sp.{json,yaml}`.
The selected mechanism and MPI EB tree spray executable built successfully.
The `pelef_mpi_support` library also built successfully to compile the touched
general MPI transport comparisons and checkpoint paths; log
`build/heptane-simple-mpi-support-build.log`. No MPI library test matrix was run.
Logs: `build/heptane-simple-mpi-configure.log`,
`build/heptane-simple-mechanism-build.log`, and
`build/heptane-simple-app{,-final}-build.log`.

Build commands used `env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build
build/heptane-simple-mpi-debug --target TARGET --parallel 2`, with targets
`pelef_selected_mechanism`, `pelef_mpi_reactive_eb_tree_3d_selected_spray`, and
`pelef_mpi_support` respectively.

```sh
build/cantera-venv/bin/python tools/prepare_spray_reference.py --fuel heptane-3sp --output-dir build/heptane-3sp-simple-reference --yaml build/heptane-3sp-reference/heptane_3sp.yaml --license build/heptane-3sp-reference/PelePhysics-license.txt --simple-transport-header build/pelec-pinned-hit-source/Submodules/PelePhysics/Mechanisms/heptane_3sp/mechanism.H
build/cantera-venv/bin/python tools/check_mechanism_provenance.py --bundle build/heptane-3sp-simple-reference/heptane_3sp.json --regenerate
python3 tools/check_simple_transport_reference.py --reference-build build/pelec-pinned-spray --candidate-build build/heptane-simple-mpi-debug --bundle build/heptane-3sp-simple-reference/heptane_3sp.json --output build/simple-transport-record-reference
```

Cantera 3.2.0 import/regeneration passes. Header tampering, coefficient tampering,
generator tampering and missing coefficient rows are rejected. Fourteen existing
Python provenance tests pass. Generated-source compilation and the six-state
compiled C++ comparison pass, including the production transport-record route:
30 values, maximum relative error `8.475e-16` against the existing `2e-12` gate.
One negative-composition kernel rejection also passes. No tolerance changed.

The same two-step 4³ evaporation adaptation was run with the new executable:

```sh
python3 tools/prepare_evaporation_reference.py --bundle build/heptane-3sp-simple-reference/heptane_3sp.json --output build/spray-evaporation-simple-transport
# From the prepared directory:
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none /usr/bin/mpiexec -n 1 /home/gtt13/projects/pelec/build/heptane-simple-mpi-debug/pelef_mpi_reactive_eb_tree_3d_selected_spray tree.nml spray.nml
```

| Particle quantity | Fortran minus pinned reference | Relative to reference change |
| --- | ---: | ---: |
| Temperature | `4.6943e-8` K | `9.9307e-6` |
| Single-drop mass | `-3.2393e-25` kg | `-2.8583e-9` |
| Axial speed | `3.0447e-8` m/s | `4.3621e-5` |

The preceding differences were -32.1% heating, +13.8% evaporation and -11.9%
deceleration. This joint transport intervention greatly reduces all three;
it does not isolate viscosity, conductivity and diffusion contributions.
Particle integration still differs: the existing Fortran frozen-gas drag
update uses an exponential while pinned Drag.H uses explicit Euler. No remaining
gap is claimed resolved, and full gas-field parity was not assessed.

Four successful application executions: full fitted run, disabled control,
one-step checkpoint and 1-to-2-rank continuation. Disabled gas, particle and
budget CSVs are byte-identical to the preceding analytical baseline. Restarted
gas CSV is byte-identical to the continuous fitted run; particles match exactly
after combining ranks by ID. One mismatched-transport restart returns nonzero
with `Gas snapshot restore failed`, publishing no result CSV.

Artifacts: `build/spray-evaporation-simple-{transport,disabled,checkpoint,restart}/`.
`transport/comparison.json` identifies the unchanged C++ reference output by
hash; `transport/integration.json` records restart equality and discrepancies.
Importer rejection results are in
`build/heptane-3sp-simple-reference/rejection-checks.json`; compiled transport
checks are in `build/simple-transport-record-reference/` (the earlier generated
loader check is `build/simple-transport-generated-reference/`).

Python syntax, project contract and diff checks pass. Existing real-comparison,
array-temporary and executable-stack warnings remain. No full test matrix,
installed-binary audit, long-time spray comparison or physical validation was run.
