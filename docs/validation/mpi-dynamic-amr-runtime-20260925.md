# Dynamic MPI AMR runtime checks after launcher recovery — 2026-09-25

Tested source: `7981fc0b93ef93e3d67c2d1c46fe3ebf23bc7072`, clean tree before
this documentation update. This supplies runtime evidence missing from the
earlier Dynamic MPI AMR, wide-halo and application build-only records.

## Launcher observation and recovery

The prior CTest PID 2018710 and launcher PID 2018713 remained live. The launcher
had affinity only to CPU 12, approximately 27.5 ms accumulated CPU time, no child
ranks and pending termination signals. CPUs 0--39 were online; its cgroup was
not frozen. After verifying its command line, only this launcher had its affinity
restored to the current tool's allowed CPU set. It exited, allowing the old CTest
to report its MPI test as **timed out (1327.67 s)**. Its two subsequent serial
tests passed; that old run is **2/3 passed, one failure**, not a successful run.
No unrelated process was signaled or rebound.

New launches used `HWLOC_COMPONENTS=-x86` and
`OMPI_MCA_hwloc_base_binding_policy=none`. A two-rank `/bin/true` launch exited 0,
then the focused numerical suite completed. These observations locate the
stall in launcher/CPU-affinity behavior; they do not establish the underlying
kernel or hardware cause. Repository launch defaults were not changed.

## Completed verification

GNU Fortran 13.3.0, Open MPI 4.1.6 (`/usr/bin/mpiexec`), Ninja Debug. Fixed test
and application targets below rebuilt successfully. New CTest run: **4/4 passed**,
zero failures/skips, 56.22 seconds total:

| Test | Seconds | Evidence |
| --- | ---: | --- |
| `unit_mpi_amr_patch_tree_ghost_3d_2` | 3.05 | Seven-layer PCM/MC/minmod serial/MPI donor parity and width-disagreement rollback, plus existing temporal/boundary cases |
| `spray_mpi_amr_dynamic_application` | 14.96 | Dynamic Gaussian filter enabled with zero static coefficients, continuous vs 2-to-4-rank restart and changed-filter context rejection |
| `spray_les_mpi_dynamic_2` | 0.65 | Existing periodic regular-grid dynamic coefficient/flux/step parity across filter/interpolation cases |
| `spray_les_mpi_amr_tree_2` | 37.51 | Root and nested-tree dynamic/static serial/MPI state, temperature, timestep and ledger parity; empty/distributed ownership; interval/filter disagreement rollback |

Existing tolerances were unchanged. The application checker uses relative
`2e-12` and absolute `1e-25` field comparison; MPI AMR state comparison is scaled
`2e-12`, temperature `1e-8 K`, boundary ledger `1e-14`. The tests are bounded
implementation checks, not physical validation or production scaling. Selected
mechanism execution, full matrix and external PeleC parity were not run.

Provenance preflight, Python syntax, project contract and whitespace checks pass.

## Commands and artifacts

Logs: `build/mpi-launch-recovery-checks/` (`smoke.log`, `build.log`, `ctest.log`,
`provenance.log`). The failed prior run remains in
`build/les-dynamic-faces-checks/ctest.log`. Application inputs, outputs and
checkpoints are under
`build/spray-amr-reactive-debug/tests/spray_mpi_amr_dynamic_application/restart_esk5absq/`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_les_mpi_dynamic test_mpi_amr_patch_tree_ghost_3d test_les_mpi_amr_tree pelef_mpi_amr_spray_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^(spray_les_mpi_dynamic_2|unit_mpi_amr_patch_tree_ghost_3d_2|spray_les_mpi_amr_tree_2|spray_mpi_amr_dynamic_application)$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py tests/application/check_mpi_amr_spray_restart.py
python3 tools/check_project_contract.py --source .
git diff --check
```
