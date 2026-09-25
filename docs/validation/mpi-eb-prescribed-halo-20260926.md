# MPI EB tree prescribed inflow halos — 2026-09-26

## Source and behavior

Implementation dda6354, parent 21acfd5, local main. The owner-local MPI EB
halo API accepts boundary_faces and inflow_state with inflow_temperature,
matching the serial prescribed inflow/outflow interface. One-layer face halos
and wide CTU halos apply the same root boundary labels. Wide corners select
the first inflow axis, as in the serial implementation.

Before donor exchange, collective checks validate optional-data presence,
shapes, finite values, supported labels and conflict with Mach controls.
The min/max agreement vector includes label presence, six face codes,
inflow-state presence and all conserved inlet values, as well as existing
temperature controls. State and temperature recovery remain transactional;
invalid inlet EOS on an owning rank rejects publication on all ranks.
The exchange continues to fetch only requested donor cells.

## Completed checks

GNU Fortran 13.3.0, Open MPI 4.1.6 via /usr/bin/mpiexec, Ninja Debug,
Python 3.12.3; existing build/spray-amr-reactive-debug. No backend changes.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_eb_tree_halo_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_mpi_eb_tree_halo_3d_3$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**1/1 passed, 0 failed, 0 skipped**, 3.61 s test / 3.66 s total.
The three-rank fixture uses five patches and both distributed ownership and
all patches on rank zero, leaving two ranks idle. New prescribed-inlet checks
cover exact serial state/temperature agreement for one-layer and six-layer
halos, exact volume fractions, ancestor temporal interpolation, rank-disagreeing
inlet state rollback and collective invalid-EOS rejection. Existing outflow,
geometry, temporal-context and bad-donor checks also pass. No failed build or
test attempts. No tolerances or golden data changed.

Provenance, syntax, project contract and whitespace checks passed. Mechanism
SHA: 0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4.
Artifacts: build/mpi-eb-prescribed-halo-checks/build.log and LastTest.log.
Existing compiler warnings remain.

## Remaining scope

MPI hydro/transport/full-step forwarding and application/checkpoint input
remain pending. This check does not establish those paths or periodic/wall
mapping. No full matrix, installed audit, external solver comparison or
physical validation was run. Complete PeleC parity remains unproven.
