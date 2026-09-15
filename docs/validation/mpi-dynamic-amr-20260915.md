# MPI dynamic 3D AMR — 2026-09-15

Source: parent `3604754706bec38f231f66f811c996e4ab19bd4e` plus this
commit's MPI driver, application labels, test runner/registration and docs.
GNU Fortran 13.3 Debug with `/usr/bin/mpifort` and `/usr/bin/mpiexec`;
`build/port-regrid-mpi` uses only OpenMPI `libmpi.so.40` (verified with `ldd`).
This is focused local evidence, not full MPI/selected-chemistry qualification.

## Implemented behavior

At tagging steps the application gathers the synchronized hierarchy, runs the
same conservative serial remap, broadcasts the accepted geometry/counter and
rebuilds sparse coarse/fine slab ownership before redistribution. Existing
rank-local hydro/chemistry/transport operators advance the resulting hierarchy.
An initially empty fine owner can acquire cells, and empty fine owners are
retained. Tagging interval, buffer and threshold must agree across ranks.

Dynamic checkpoint reading restores the patch before distribution; the format
contains no rank ownership. Fixed chemistry-enabled restart remains rejected,
and the selected chemistry/transport context continues through the shared
reader/writer. Regrid and checkpoint operations aggregate on rank zero, so root
memory and communication remain scalability limitations. Multiple patches,
more than two levels, physical coarse boundaries and general 3D EB remain open.

## Completed focused checks

The Python integration checker ran five launches: a one-rank reference,
a four-rank checkpoint stop, a three-rank restart, and two three-rank rejected
restarts (changed tagging threshold and truncated file). All expected outcomes
passed. On the 8-cubed case, the patch moves at step zero and expands at step
four; the checkpoint is after step five. Both completed runs reach 2 microseconds
and report two regrids. The four-rank case includes an initially empty fine owner
that acquires cells and an owner that remains empty.

Coarse and fine CSVs match byte-for-byte between the one-rank reference and
four-to-three-rank restart. A separate comparison against the previous serial
Debug reference also passes:

- Coarse SHA256: `b5cd4361685a1d4cea13331d095d1a7f92bdd494ccb72695384c8907732b8808`
- Fine SHA256: `aa36eb4ea6ee2bf980b0defe42b75ddff2242155834b9b0e26d3af924d677687`

Artifacts: `build/port-regrid-mpi/dynamic-restart/` contains all inputs,
logs, checkpoint and CSVs. The executable was rebuilt after the final label and
indentation edits. Project contract, Python syntax, provenance preflight and
`git diff --check` pass. FFCM1_Red source SHA256 remains
`ddbba8c97d929530cfb213da45485cca4f979dccf97c1063721bd855d7b3f259`.

```sh
env -i PATH=/usr/bin:/bin LC_ALL=C cmake --build build/port-regrid-mpi --target pelef_mpi_amr_reactive_3d --parallel 8
python3 tools/check_dynamic_amr_3d_restart.py --executable build/port-regrid-mpi/pelef_mpi_amr_reactive_3d --mpiexec /usr/bin/mpiexec --work-dir build/port-regrid-mpi/dynamic-restart
python3 tools/compare_amr_reactive_3d_restart.py --reference-coarse build/port-spray-debug/tests/amr_dynamic_restart_3d/reference_coarse.csv --reference-fine build/port-spray-debug/tests/amr_dynamic_restart_3d/reference_fine.csv --restart-coarse build/port-regrid-mpi/dynamic-restart/resumed_coarse.csv --restart-fine build/port-regrid-mpi/dynamic-restart/resumed_fine.csv
python3 -m py_compile tools/check_dynamic_amr_3d_restart.py
python3 tools/check_project_contract.py --source .
python3 tools/check_mechanism_provenance.py --bundle build/port-methanol-reference/ffcm1_red.json --source build/port-methanol-reference/ffcm1_red.yaml
```

CTest registration `regression_mpi_dynamic_amr_3d_restart` invokes this checker
with configured MPI flags. The runtime check here was invoked directly because
this build has tests disabled; no CTest suite is claimed. Selected chemistry and
transport with dynamic MPI regridding were not runtime-tested in this increment.
