# Jet failure: isolated symbol-enabled replay — 2026-09-27

The original and retry-enabled Jet runs remain failed; no solver repair is
claimed here. The previous same-binary GDB replay localized terminal failure
to time `0.0004287942691061733 s`, accepted step 2436 and next injection 2144.
To identify the rejecting operation, source
`c4f6defc1e8e8212c9675d9d76bca073d1234306` was exported into
`build/jet-failure-source` and rebuilt separately in
`build/jet-failure-symbols`. The original Release binary was not modified.

## Completed preparation

GNU Fortran 13.3.0 / OpenMPI 4.1.6, selected decane Simple transport, MPI/spray
enabled, tests/Cantera-reference/SUNDIALS disabled. A clean configure initially
found the environment's Intel MPI launcher with the OpenMPI compiler; the
project's MPI installation check correctly rejected it. Explicit
`MPIEXEC_EXECUTABLE=/usr/bin/mpiexec` fixed configuration. The final selected
spray executable built and its MPI runtime audit passed (ABI 40).

```sh
mkdir -p build/jet-failure-source
git archive c4f6defc | tar -x -C build/jet-failure-source
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake -S build/jet-failure-source -B build/jet-failure-symbols -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran -DMPI_Fortran_COMPILER=/usr/bin/mpifort -DMPIEXEC_EXECUTABLE=/usr/bin/mpiexec -DPELEF_ENABLE_MPI=ON -DPELEF_ENABLE_SPRAY=ON -DPELEF_ENABLE_TESTS=OFF -DPELEF_ENABLE_CANTERA_REFERENCE=OFF -DPELEF_ENABLE_SUNDIALS=OFF -DPELEF_MECHANISM_BUNDLE=/home/gtt13/projects/pelec/build/decane-3sp-simple-reference/decane_3sp.json -DPELEF_MECHANISM_SOURCE=/home/gtt13/projects/pelec/build/decane-3sp-simple-reference/decane_3sp.yaml
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/jet-failure-symbols --target pelef_mpi_reactive_eb_tree_3d_selected_spray --parallel 2 > build/jet-failure-symbols-build.log 2>&1
python3 -m py_compile build/jet-symbol-replay/run.py
python3 tools/check_mechanism_provenance.py --bundle build/decane-3sp-simple-reference/decane_3sp.json
python3 tools/check_project_contract.py --source .
git diff --check
```

The final configuration, build, syntax, provenance and contract checks passed.
The initial rejected configure is retained in `build/jet-failure-configure.log`;
the corrected configure is in `build/jet-failure-configure-final.log`.
Compiler warnings and an executable-stack linker warning remain.

## Replay and evidence boundary

`python3 build/jet-symbol-replay/run.py` launched GDB at
`2026-09-27T10:29:26.111355+00:00`. The live driver, GDB and solver were checked
after launch (PIDs 642903, 642904, 642914 at that observation). Consult
`build/jet-symbol-replay/execution.json` and the actual process for later
status; this record does not report a completed replay.

The replay copies the earlier diagnostic input and restarts retry-run
`state.chk.step2400`. `provenance.json` records input and checkpoint hashes,
source commit and both binary identities. `trace.gdb` stops at the first
observed rejected common gas or particle operation and records its stack and
bounded local values. A first rejection need not be the terminal retry-limit
failure; further evidence must establish that relationship.

Diagnostic binary SHA256:
`d7048d385e463dd4a7f20d4d97644b6c299b67a7cae294ee91b524f3594f866c`.
Original binary SHA256:
`31b374ef2bdc0385bb464e0313a50baaa1548cd0a6692d01c7d3abc0c63bcdda`.
The source is the same but compiler options differ: original Release `-O3`,
diagnostic RelWithDebInfo `-O2 -g`. Reproduction is not assumed from source
identity alone. The original failed artifacts and current main checkout are
preserved. No full configuration matrix or new physical validation was run.
