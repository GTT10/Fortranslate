# MPI EB recursive hydro — 2026-09-18

Implementation: `1f810e3d2e61716a52aff336824a1482705bb413`, based on
`5dd142f`. Adds `advance_mpi_eb_tree_subcycle_hydro_3d` to `pelef_mpi_support`.

## Implemented behavior

Each level advances an Euler hydro interval from an owner-local snapshot. Its
provisional endpoint first receives same-level peer flux corrections. Children
then advance refinement-ratio substeps using independently interpolated ancestor
time endpoints. Effective child fluxes accumulate over the parent interval;
level-scoped synchronization restricts, mirrors, refluxes and averages them.
The final child interval uses the remaining parent time to avoid endpoint drift.

The MPI flux synchronizer now has a collectively agreed `peers_only` mode. It
updates provisional peer states and fluxes without restricting or averaging
descendants at another time. The default synchronization behavior is retained.
The recursive operator returns an effective root boundary ledger and an MPI-summed
exposed wall ledger. U/T and both ledgers publish together only after success.
Geometry stays replicated; persistent state and fluxes remain owner-local.

This is a hydro library API. MPI molecular-transport subcycling, coupled/application
dispatch and subcycle-scaled automatic intervals remain separate work.

## Completed focused verification

Fresh GNU Fortran 13.3.0 Debug/Ninja build, `/usr/bin/gfortran`, Open MPI wrapper
`/usr/bin/mpifort` and launcher `/usr/bin/mpiexec`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake -S . -B build/mpi-eb-tree-subcycle-hydro-debug -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran -DMPI_Fortran_COMPILER=/usr/bin/mpifort -DMPIEXEC_EXECUTABLE=/usr/bin/mpiexec -DPELEF_ENABLE_MPI=ON -DPELEF_ENABLE_TESTS=ON
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/mpi-eb-tree-subcycle-hydro-debug --target test_mpi_eb_tree_subcycle_hydro_3d test_mpi_eb_tree_flux_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/mpi-eb-tree-subcycle-hydro-debug -R '^unit_mpi_eb_tree_(subcycle_hydro|flux)_3d_[23]$' --output-on-failure --output-log build/mpi-eb-tree-subcycle-hydro-checks/ctest.log
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**4/4 tests passed**, no failures/skips, total **9.90 s**. Existing two/three-rank
flux checks took 0.82/0.84 s; recursive hydro checks took 4.06/4.14 s.

The hydro fixture covers x/y/z plane cuts, four/five-node three-level trees,
one-sided and paired deep refinement, distributed and single-owner layouts,
and empty ranks. MPI state, temperature and boundary ledger match the serial
recursive operator exactly. Wall ledgers agree within the unchanged 3e-14 bound.
The **24** MPI inventory comparisons have maximum normalized residual
**6.508800774330761e-15**, below the unchanged 3e-14 conservation bound.
Invalid solver, rank-disagreeing timestep and invalid input density reject
collectively while preserving U/T and sentinel ledgers. Existing common and
level-scoped flux synchronization checks also pass.

Configure/build, provenance, Python syntax, project contract and whitespace
checks pass. No failed build/test attempt occurred; compiler warnings remain.
No Python source changed. This does not qualify coupled/application subcycling,
temporal convergence, selected MPI, installed binaries, scaling, a full matrix
or physical validation.

Logs: `build/mpi-eb-tree-subcycle-hydro-checks/{configure,build,ctest}.log`.
