# MOL interface flux connection — 2026-09-27

`mol_riemann_flux_3d` connects the independent MOL reconstruction to the existing
mixture PeleC Riemann solver. Inputs are four cell primitives, two central
sound speeds, direction, EB neighbor connectivity and reconstruction order.
First order skips outer donors; second order uses characteristic MOL slopes.
The output is intensive conservative flux in global momentum order, published
only after complete success.

Pinned `Source/MOL.cpp::pc_compute_hyp_mol_flux` passes the average of cell-center
sound speeds to the Riemann solver. The wrapper follows this convention and
the exact transverse-direction ordering. Unlike the CTU wrapper, it does not
replace the Riemann species flux with passive upwinding: the pinned MOL routine
retains the species flux returned by `riemann`. Reference revision remains
`bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

## Verification

GNU Fortran 13.3.0, existing heptane Debug core archive:

```sh
gfortran -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -I build/heptane-simple-mpi-debug/modules -J build/mol-riemann-unit src/hydro/mol_riemann_3d_mod.F90 tests/unit/test_mol_riemann_3d.F90 build/heptane-simple-mpi-debug/libpelef_core.a -o build/mol-riemann-unit/check
build/mol-riemann-unit/check
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/heptane-simple-mpi-debug --target pelef_core --parallel 2 > build/mol-riemann-build.log 2>&1
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

Final commands pass. The initial test compile failed because an ALLOCATE
statement referenced an array being allocated by that same statement; splitting
the allocations fixed the test. One executable passes 31 assertions covering
uniform first-order and linear second-order exact physical fluxes in all three
directions, species/mass closure, inactive NaN donors and active-donor rollback.
The fixed relative gate is `2e-13`. CTest registration was added; no CTest
matrix or compiled external-reference comparison ran. Only the core build was
executed; the selected source list also contains the new module.

Open-face selection, EB centroid interpolation, face area multiplication,
embedded pressure forces, rate assembly, time integration and application
routing remain separate work. The optional Lax-Friedrichs flux and rotating
frame branches are not part of this wrapper. This does not establish a usable
MOL application or full Spray-EB parity.

The independent jet continuation was still live at 7 minutes 50 seconds,
99.9% CPU, and had published `state.chk.step2100.gas` in
`build/decane-jet-retry-run`. Full completion remains unverified.
