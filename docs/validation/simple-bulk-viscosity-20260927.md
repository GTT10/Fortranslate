# Simple bulk viscosity — 2026-09-27

Implementation: `e5e6cf23`.

The fitted Simple transport path now computes bulk viscosity from the pinned
PelePhysics rotational-relaxation formula: species heat capacities, molecular
geometry, Lennard-Jones well depth, rotational relaxation and fitted shear
viscosity. The mixture uses the reference three-quarter-power mixing rule.
The 3D stress adds bulk viscosity times velocity divergence on its normal
component; bulk viscous work enters the energy flux. Face bulk viscosity uses
the reference default harmonic mean. The viscosity timestep limit includes
`(4*mu/3 + bulk)/rho`.

This applies to the fitted transport backend in the 3D face/timestep APIs.
Prescribed constant transport and the existing analytical backend retain their
previous behavior. It is not a claim of bulk-viscosity integration in 1D/2D.
Fitted transport checkpoint tags advance to `PELEF_SIMPLE_TRANSPORT_V2`, rejecting
snapshots from the prior zero-bulk fitted model. Coefficient-free layouts stay
unchanged.

## Diagnosis and verification

Disabling gas viscosity, conduction and species diffusion while retaining
particle Simple transport reduced the preceding axial-momentum discrepancy
from `5.156e-6` to `3.859e-7` (relative to the reference field maximum). This
motivated checking the viscous stress; it did not by itself prove a bulk-only
cause. Source inspection then found the missing bulk term.

The valid intervention is `build/spray-no-gas-transport2/`, containing matched
C++ and Fortran inputs, logs, outputs and `full-comparison.json`. The earlier
`build/spray-no-gas-transport/` experiment is **excluded**: its input edit also
disabled particle Simple transport. That error was corrected in a fresh run.

GNU Fortran 13.3.0 Debug / OpenMPI 4.1.6, MPI/spray ON, tests/SUNDIALS OFF:

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/heptane-simple-mpi-debug --target pelef_mpi_reactive_eb_tree_3d_selected_spray --parallel 2
python3 tools/check_simple_transport_reference.py --reference-build build/pelec-pinned-spray --candidate-build build/heptane-simple-mpi-debug --bundle build/heptane-3sp-simple-reference/heptane_3sp.json --bulk --output build/simple-bulk-flux-reference
```

Build passed (`build/simple-bulk-{,final-}build.log`). Six states / 36 transport
values, now including bulk viscosity, agree with compiled pinned C++ to
`1.103e-15` maximum relative error against the unchanged `2e-12` gate. Isotropic
compression stress, bulk viscous work and the longitudinal-viscosity timestep
bound also pass at `1e-12` relative tolerance. The existing negative-composition
kernel rejection passes. Artifacts, source identities and exact compiler commands:
`build/simple-bulk-flux-reference/`; the earlier bulk-coefficient-only check
is `build/simple-bulk-reference/`.

## Evaporation application

The full two-step 4³ case completes with the same pinned reference output as
the [preceding cell-source comparison](spray-cell-source-20260927.md).

| Gas field | Maximum difference / maximum reference magnitude |
| --- | ---: |
| Axial momentum density | `2.654e-6` |
| Transverse momentum density | `2.670e-6` |
| Total energy density | `2.360e-10` |
| Fuel density | `1.975e-12` |
| Total density | `7.405e-13` |
| Stored temperature | `2.218e-9` |

Bulk viscosity reduces the momentum discrepancy but does not remove it.
Stored-temperature and gas-source differences remain open; no full parity or
physical validation is claimed.

```sh
python3 tools/prepare_evaporation_reference.py --bundle build/heptane-3sp-simple-reference/heptane_3sp.json --output build/spray-bulk-transport
# From each prepared run directory, use -n 1 (restart: -n 2):
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none /usr/bin/mpiexec -n 1 /home/gtt13/projects/pelec/build/heptane-simple-mpi-debug/pelef_mpi_reactive_eb_tree_3d_selected_spray tree.nml spray.nml
```

Four post-change application runs complete: full evolution, checkpoint,
one-to-two-rank continuation and gas-transport-disabled control. Restart gas CSV
is byte-identical and particle rows match continuous output exactly. Disabled
transport reproduces its pre-change control exactly. One old-model restart
is rejected with `Gas snapshot restore failed`, publishing no CSV.
Artifacts: `build/spray-bulk-{transport,checkpoint,restart,disabled-control}/`;
comparison hashes and restart results are under the `transport` run directory.

Provenance preflight, Python syntax, project contract and diff checks pass.
No CTest matrix or installed-binary audit was run.
