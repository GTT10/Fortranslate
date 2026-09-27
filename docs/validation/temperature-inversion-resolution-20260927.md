# Temperature inversion at binary64 resolution — 2026-09-27

Implementation `d1488573` adds a representable-temperature stopping condition
to internal-energy inversion. The existing energy residual tolerance remains
`1e-11 * max(1,abs(target))`. If the evaluated lower/upper energies enclose
the target and the temperatures are adjacent representable values, the solver
returns the endpoint with smaller energy residual. No intervening binary64
temperature exists. The new termination is disabled across a NASA polynomial
branch boundary, where a coefficient discontinuity could mimic a root bracket.

## Failure evidence

The symbol-enabled Jet replay reached common SDC's final nonreacting state
recovery at exported-source line 652, then GDB itself aborted while printing
locals (`exit_code=-6`). Its captured stack remains in
`build/jet-symbol-replay/trace.log`; this was not a successful full replay.

A separate failure-only logging replay of the same source and `-O2 -g` build
completed with solver exit 1 after `263.12831960618496 s`, ending at
`2026-09-27T10:36:44.072119+00:00`. The first reported failure was root cell
`(2,1,2)`, iteration 2, interval `2e-7 s`. Density was
`46.643112388154066 kg/m³`, species densities were
`[0.46297573307930956,0,46.180136655074932] kg/m³`, and total energy density
was `139.69283779491434 J/m³`. Density/species closure was intact.

The extracted single-cell thermodynamic probe reproduces inversion failure
at specific internal energy `1.0351682798876758 J/kg`. By iteration 97 the
temperature bracket was `[436.85902346503832,436.85902346503838] K`, adjacent
binary64 values. Residual `2.6029400856941720e-11 J/kg` exceeded tolerance
`1.0351682798876757e-11 J/kg`; the next midpoint rounded to the same upper
endpoint through iteration 100. This establishes a resolution-induced failure
in the replay, rather than a low-temperature/out-of-range thermodynamic state.
The fixed probe returns `436.85902346503838 K` in four iterations.

## Completed checks and artifacts

GNU Fortran 13.3.0 / OpenMPI 4.1.6. The current fixed-core build and manually
linked `test_mixture_thermo` pass, including existing property/inversion/invalid
input checks and four new pure-N2 near-zero targets. For each new target,
energies evaluated at the returned temperature's immediate neighbors enclose
the target. The exact decane replay state also passes with the current source
compiled into `fixed-mixture.o` ahead of the selected runtime archive.

Key commands:

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/heptane-simple-mpi-debug --target pelef_core --parallel 2 > build/jet-cell-replay/current-core-build.log 2>&1
env -u LD_LIBRARY_PATH -u LIBRARY_PATH /usr/bin/gfortran -O0 -g -fcheck=all -ffree-line-length-none -I build/heptane-simple-mpi-debug/modules tests/unit/test_mixture_thermo.F90 build/heptane-simple-mpi-debug/libpelef_core.a -o build/jet-cell-replay/thermo-check
build/jet-cell-replay/thermo-check > build/jet-cell-replay/thermo-check.log 2>&1
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
python3 -m py_compile build/jet-cell-replay/run.py build/jet-energy-repaired/run.py
git diff --check
```

All listed checks passed. `build/jet-cell-replay/probe.F90`, `probe.log` and
`probe-fixed-final.log` retain the exact-state comparison. Its label `ORIGINAL`
names the library call; the later `ITER` lines intentionally emulate the old
algorithm and still stagnate, even in the fixed executable. Failure-only
instrumentation is confined to the exported old source, not current main.

## Bounded intervention replay

`build/jet-energy-repaired` restarts the same step-2400 checkpoint with only
the inversion change plus diagnostic logging, ending at `0.00043 s` instead
of the full case's `0.0015 s`. The selected application rebuilt and its MPI
runtime audit passed. Binary SHA256:
`7f9998254ba67171cfa9cd3aca96d00bff565e89cba2f936d65f00ba2619aa1e`.
Replay launched at `2026-09-27T10:37:58.697858+00:00`; it was running at the
last observation. Consult `execution.json` and the actual process for status.
This record does not yet claim that it crosses the failed interval, that the
original `-O3` behavior matches `-O2`, or that the full Jet case completes.
