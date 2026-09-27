# Embedded no-slip application connection — 2026-09-27

Implementation: `c74570be`, based on `4766b4c0`. The MPI gas and particle
embedded-boundary tree applications now accept `embedded_noslip=.true.` in
`eb_tree`, default false. It requires SDC, enabled transport and viscosity.
Both particle SDC wrappers forward it, including their wall-impact branches.
Heat remains independently selectable. Serial applications reject enabled
no-slip. The existing geometry cache and unconditional wall timestep cap now
serve viscosity as well as conduction. The new optional checkpoint prefix
`PELEF_EB_TREE_EMBEDDED_NOSLIP_V1` binds the setting; disabled runs retain the
previous header format.

## Completed checks

GNU Fortran 13.3.0 / OpenMPI 4.1.6, existing Debug build
`build/heptane-simple-mpi-debug`, MPI/spray ON, tests/SUNDIALS OFF.

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/heptane-simple-mpi-debug --target pelef_mpi_reactive_eb_tree_3d_selected_spray pelef_mpi_reactive_eb_tree_3d_selected pelef_reactive_eb_tree_3d_selected --parallel 2 > build/eb-noslip-app-build.log 2>&1
python3 -m py_compile build/eb-noslip-app/run_checks.py
python3 build/eb-noslip-app/run_checks.py > build/eb-noslip-app/driver.log 2>&1
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

All commands passed. The application checks use the checked-in
`cases/spray_eb_wall_heat_3d` fixture, enabling no-slip in generated inputs.
Eight bounded launches cover combined heat/traction, checkpoint stop, restart
with 1→2 ranks, changed-setting rejection, adiabatic common and recursive
particle steps, a deliberately small viscous CFL cap, and the prior heat-only
setting. The changed restart fails as intended before output. Gas CSV bytes
match uninterrupted execution; a separate comparison of all particle CSVs
across both ranks confirms the same four particle rows exactly.

For combined heat and traction, wall heat is `9.561100114981405e-9 J` and
energy residual is `-2.710505431213761e-19 J`. For both adiabatic paths,
wall heat is exactly zero, energy residual is `-1.8973538018496328e-19 J`, and
the tangential z wall impulse is `-4.565986383147719e-19 N s`. This tiny
particle-induced impulse demonstrates routing; the separate moving-gas unit
establishes the larger analytic drag check. With `transport_cfl=1e-8` and
automatic timestepping disabled, the first interval is capped from `1e-7 s`
to `1.273132029911567e-11 s`. The disabled-no-slip heat-only gas CSV is byte
identical to the preceding application artifact.

Artifacts: `build/eb-noslip-app/{result.json,driver.log,run_checks.py}` and
the named per-run directories. The selected spray executable SHA-256 is
`097b999d92c63aab1fab49482e258fce6478089651f410ce5ef4ad249e3ef88a`.
The report records the pre-commit base; the final implementation commit above
identifies the tested source changes. Linker executable-stack and MPI X11
messages remain nonfatal. No repository Python source changed.

This is application integration evidence, not full pinned Spray-EB parity or
physical validation. The recursive particle fixture has one level; multilevel
gas traction was checked separately. No full CTest matrix was run. The existing
full decane jet process remained live and unchanged.
