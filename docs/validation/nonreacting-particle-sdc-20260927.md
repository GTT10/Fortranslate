# Reaction-free selected particle SDC — 2026-09-27

Implementation: `be25365087344799b8798f8ccd0aadddd9b0bee7`.
The shared selected EB tree front end now permits explicitly configured
zero-reaction mechanisms when chemistry is disabled. Chemistry ON is rejected
before state evolution/output. This connects the existing MPI particle SDC
implementation to the pinned reaction-free NC7H16/O2/N2 bundle.

The evaporation preparation tool also writes Fortran gas/particle inputs and
an all-positive nodal level set. The candidate uses a single fully regular 4³
patch, periodic boundaries, PLM CTU, two SDC iterations, and two `1e-7` s steps.
Gas density is derived from the source molecular weights, 79/21 mole fractions,
400 K and 101325 Pa. Four moving 310 K parcels carry 1e-10 kg initially.
Reference input now explicitly selects PLM order 4 to align reconstruction.

## Completed checks

Fresh `build/heptane-spray-mpi-debug`: GNU Fortran 13.3.0 Debug, OpenMPI 4.1.6,
MPI/spray ON, tests/SUNDIALS OFF; pinned bundle/source paths as in the serial
[integration record](nonreacting-spray-app-20260927.md). Configured with explicit
`/usr/bin/gfortran`, `/usr/bin/mpifort`, `/usr/bin/mpiexec`; built only
`pelef_mpi_reactive_eb_tree_3d_selected_spray` successfully.

```sh
python3 tools/prepare_evaporation_reference.py --bundle build/heptane-3sp-reference/heptane_3sp.json --output build/spray-evaporation-sdc-run2
# From the prepared directory:
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none /usr/bin/mpiexec -n 1 /home/gtt13/projects/pelec/build/heptane-spray-mpi-debug/pelef_mpi_reactive_eb_tree_3d_selected_spray tree.nml spray.nml
```

One successful two-step SDC run, plus one expected chemistry-ON rejection.
Composite absolute inventory residual: `1.3552527156068805e-19`.
Gas NC7H16 gain: `4.516830370340473e-15` kg; agrees with liquid loss within
relative `1e-9`. Four parcels survive, with active heating, drag and evaporation.
No comparison tolerance was changed. Chemistry-ON rejection published no gas CSV.
Provenance preflight, Python syntax, project contract and diff checks passed.

Artifacts: `build/spray-evaporation-sdc-run2/`, including `exchange.json`, gas,
particle and inventory CSVs and logs; rejection: `build/heptane-sdc-reject/`.
Build/configure logs: `build/heptane-spray-mpi-{configure,build}.log`.
The initial `build/spray-evaporation-sdc-run/` failed at the existing explicit
boundary guard because `particle_outflow` was absent. Adding that required flag
enables the explicit boundary path; periodic faces still remain periodic.
Headless MPI emitted X authorization warnings, but the completed solver run
returned zero and produced the stated outputs.

No CTest matrix, multiple-rank continuation or full reference spray trajectory
comparison is claimed. The pinned reference build remains in progress in
`build/pelec-pinned-spray`; its log is `build/pelec-pinned-spray-build.log`.
