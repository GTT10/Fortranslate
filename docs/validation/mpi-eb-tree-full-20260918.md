# MPI EB tree transport and reacting split steps — 2026-09-18

Implementation: `3a3dcfb6dd8010cd38bdacbc70a8c3b7624eb727`, based on
`69fc049`. New modules `mpi_eb_tree_transport_3d_mod` and
`mpi_eb_tree_full_3d_mod` are part of `pelef_mpi_support`.

## Implemented behavior

`advance_mpi_eb_tree_transport_3d` connects owner-local molecular transport to
private volume-weighted averaging, sparse face donors and conservative flux
synchronization. Physical domain faces are closed to diffusion; patch edges
inside the domain use donors. Viscosity, conduction, diffusion and barodiffusion
flags must agree across ranks, as must the interval. The effective root flux
gives the boundary ledger, and MPI minimum reduction gives the accepted limiter.

`react_mpi_eb_tree_3d` advances only exposed fluid cells, excluding regions
covered by children. It uses the existing explicit/implicit chemistry kernel,
then averages down. Composite inventories before/after give the chemical
exchange vector. Interval, tolerances and integrator controls are checked
collectively. Species/reaction/transport tables remain the caller's common
mechanism contract, as in the other MPI library APIs.

`advance_mpi_eb_tree_full_3d` performs R/2 T/2 H T/2 R/2 on private owner-local
fields. Chemistry and transport can be enabled independently; all stage fields,
boundary/wall/chemical inventories and minimum limiter publish together only
on success. Rank-inconsistent stage switches are rejected before branching.
No full state gather is used. Geometry remains replicated. Temporal halos and
subcycling, selected-target wiring, MPI checkpoint/application integration and
physical qualification remain open.

## Completed focused verification

Fresh GNU Fortran 13.3.0 Debug/Ninja build; compiler `/usr/bin/gfortran`, MPI
wrapper `/usr/bin/mpifort`, launcher `/usr/bin/mpiexec`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake -S . -B build/mpi-eb-tree-transport-debug -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran -DMPI_Fortran_COMPILER=/usr/bin/mpifort -DMPIEXEC_EXECUTABLE=/usr/bin/mpiexec -DPELEF_ENABLE_MPI=ON -DPELEF_ENABLE_TESTS=ON
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/mpi-eb-tree-transport-debug --target test_mpi_eb_tree_transport_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/mpi-eb-tree-transport-debug --target test_mpi_eb_tree_transport_3d test_mpi_eb_tree_full_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/mpi-eb-tree-transport-debug -R '^unit_mpi_eb_tree_(transport|full)_3d_[23]$' --output-on-failure --output-log build/mpi-eb-tree-transport-checks/ctest.log
python3 -m py_compile tools/check_project_contract.py tools/check_mechanism_provenance.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**3/3 tests passed**, no failures/skips, total **38.58 s**:

- Transport at two/three ranks: 1.82 s each. Three cut-plane orientations,
  four/five patches and distributed/single-owner layouts give 24 inventory
  checks, including empty ranks. U/T, boundary ledger and minimum limiter
  match serial exactly. Maximum normalized inventory residual is
  **6.798092530336388e-15**, below the unchanged 3e-14 bound. Invalid interval,
  rank-inconsistent interval/flags and deep nonphysical state preserve fields,
  boundary ledger and limiter.
- Coupled step at two ranks: 34.92 s. Five cases cover all four chemistry/
  transport switches with explicit chemistry and fully enabled implicit
  chemistry; the latter uses one owner and one empty rank. U/T, boundary
  ledger and limiter match serial exactly. Wall and chemical ledgers agree
  within 3e-14 normalized by initial inventory. Maximum full-step normalized
  budget residual is **1.8961134661160636e-15**. An invalid hydro solver after
  successful reaction/transport halves preserves all fields and diagnostic
  sentinels. Different chemistry switches are rejected collectively.

Configure/build succeeded. Compiler warnings remain; no warning-free claim is
made. Provenance, Python syntax, project contract and whitespace checks pass.
No full matrix, selected-mechanism, installed-binary, application/restart,
performance or physical-validation run was added.

Logs: `build/mpi-eb-tree-transport-checks/{configure,build,ctest}.log`.
Detailed output: `build/mpi-eb-tree-transport-debug/Testing/Temporary/LastTest.log`.
