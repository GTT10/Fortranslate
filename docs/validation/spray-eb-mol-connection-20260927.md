# Pinned Spray-EB case selects MOL — 2026-09-27

Implementation `4843dea5` changes the preparer's generated `eb_tree` input to
`mol_mode=.true.` and records `hydrodynamic_update="mol"` in provenance. Pinned
`Exec/RegTests/Spray-EB/example.inp` has `pelec.do_mol = 1`. The source SHA checks,
64³ default grid, prescribed final time and five initial particle records are
unchanged. The adaptation list now records the additional positive-face slope
guard instead of the superseded CTU/SDC discrepancy.

## Completed execution and checks

The existing selected heptane Simple-transport MPI/spray Debug executable
(GNU Fortran 13.3.0 / OpenMPI 4.1.6, tests/SUNDIALS disabled) has SHA256
`160bc2ea44c7a045d1a1004b314e57cc3c86fc783b4bb8a6fed4449f53b4175f`.
It was built for source `19d167b8` in the preceding symmetry repair; this change
only changes the input preparer and documentation.

```sh
python3 tools/prepare_spray_eb_reference.py --bundle build/heptane-3sp-simple-reference/heptane_3sp.json --reference-source build/pelec-pinned-hit-source --output build/spray-eb-mol-connection --cells 8 --final-time 2e-7
cd build/spray-eb-mol-connection
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none /usr/bin/mpiexec -n 2 /home/gtt13/projects/pelec/build/heptane-simple-mpi-debug/pelef_mpi_reactive_eb_tree_3d_selected_spray tree.nml spray.nml > fortran.log 2>&1
cd ../..
python3 build/spray-eb-mol-connection/verify.py
python3 tools/prepare_spray_eb_reference.py --bundle build/heptane-3sp-simple-reference/heptane_3sp.json --reference-source build/pelec-pinned-hit-source --output build/spray-eb-mol-reference-inputs
python3 -m py_compile tools/prepare_spray_eb_reference.py build/spray-eb-mol-connection/verify.py
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

Final checks passed. The short case accepts two 100 ns steps, processes five
initial events and retains five particles. Gas and particle combined inventory
has maximum component residual divided by `max(1, inventory/exchange scales)`
of `3.349016974711027e-17`. Injected mass is
`9.32534655328326612e-11 kg`. The initial checker used exact equality against
the Python-prepared total and failed by one binary64 ULP
(`1.2924697071141057e-26 kg`); the final checker permits exactly one ULP for
that cross-language sum. Solver acceptance and conservation tolerances were
not changed. `result.json` records the final budgets and binary identity.

The default 64³ inputs were generated in `build/spray-eb-mol-reference-inputs`
but not executed. The 8³ result is a short connection check; it does not
establish full-duration Spray-EB parity, particle wall-impact fidelity or
arbitrary-cut-cell positivity. Full configuration matrices and upstream
time-history comparison were not run. Earlier CTU artifacts remain intact.
