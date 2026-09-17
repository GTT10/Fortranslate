# Temperature-adapted MPI AMR spray runs — 2026-09-18

Implementation: `488b0607f36d1ab38779c198ab2b0bf27b79b8a0`, based on
`04ac177`. The shared fixed/selected AMR spray driver now invokes the existing
composite-temperature tagger at configured completed-step intervals. Changed
topology uses transactional gas/particle remapping and round-robin owners.
The input tree initializes the run; adaptation is processed before injection
and saving, including at final time. Prescribed and automatic regridding are
mutually exclusive. Levels, ratios, thresholds, tiles, buffer and node capacity
bind the checkpoint context.

Restart restores the adaptive tree and checks its depth/node/ratio bounds.
The saved step has already completed adaptation, so checks resume at the next
due completed step. No extra mutable adaptation counter is needed. This is
temperature-based tagging, not particle-density tagging or load balancing.

## Completed focused verification

Reused GNU Fortran 13.3.0 Debug/Ninja fixed tests-enabled and selected H2/O2
tests-disabled builds, `/usr/bin/gfortran`, Open MPI `/usr/bin/mpifort` and
`/usr/bin/mpiexec`. Compiler/MPI/backend configuration did not change.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_amr_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_amr_spray_3d_selected --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_amr_(adaptation|ports)$' --output-on-failure --output-log build/spray-amr-adaptation-checks/ctest.log
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tests/application/check_mpi_amr_spray_restart.py --exe build/spray-amr-reactive-selected-debug/pelef_mpi_amr_spray_3d_selected --mpiexec /usr/bin/mpiexec --source /home/gtt13/projects/pelec --work /home/gtt13/projects/pelec/build/spray-amr-adaptation-checks/selected --selected --multiport --adaptive
python3 -m py_compile tests/application/check_mpi_amr_spray_restart.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Fixed CTest **2/2 passed**, no skips: adaptive restart **7.94 s**, prescribed
multiport regression **4.15 s**, **12.14 s** total. Selected adaptive checker
also passed. A localized 100 K hotspot with a 1005 K threshold replaces the
initial corner-refined tree with **9 nodes / 120 exposed cells**, containing
both levels zero and one. Logs show checks at steps 1, 2 and 3 in the full run,
and only steps 2 and 3 after restarting the step-1 snapshot on four ranks.

Both frontends compare two-rank continuous and two-to-four-rank continuation
gas/particle CSVs at relative 2e-12 / absolute 1e-25, preserving 16 parcels,
port IDs, differing compositions and four injection events. Context changes,
truncated particle files and immutable-prefix overwrite are expected rejection
checks. Prescribed multiport continuation remains covered by the regression.

No build/test attempt failed. Provenance, Python syntax, project contract and
whitespace checks passed; compiler warnings remain. Gas chemistry/transport
are disabled in these CLI fixtures. No independent adaptive conservation/order
oracle, installed audit, load-balance study, scaling or physical validation was
run. Existing tagger/remap tests provide separate lower-level evidence.

Logs: `build/spray-amr-adaptation-checks/{build,build-final,build-selected,ctest,test-selected}.log`.
Fixed artifacts: `build/spray-amr-reactive-debug/tests/spray_mpi_amr_adaptation/restart_*`.
Selected artifacts: `build/spray-amr-adaptation-checks/selected/restart_*`.
