# Fortranslate / PeleF

PeleF is an independent Modern Fortran reimplementation of selected PeleC
algorithms. It is not a mechanical C++ translation or an official Pele Suite
project. The pinned PeleC reference is
`bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`; recursive dependencies are recorded
in [the reference manifest](references/pelec_baseline.json).

## Current scope

The `0.245.0` milestone contains reacting-flow, molecular-transport, MPI,
selected-mechanism, and bounded AMR/embedded-boundary implementations. The
latest increment persists transport-aware fixed/selected static 3D AMR
checkpoints. This is an integration candidate, not a declaration of complete
PeleC parity or production diesel-spray capability.

Start with [current status and evidence](docs/current_status.md). In
particular, static periodic two-level 3D AMR is not dynamic general 3D AMR,
planar single-level 3D embedded boundaries are not arbitrary 3D geometry, and
H2/O2 chemistry tests are not detailed-hydrocarbon ignition validation.

## Build and test

Use a Fortran 2018 compiler, CMake 3.23 or newer, Python 3.11 or newer, and
Ninja or Make. The optional preset commands require CMake 3.25 or newer
because `CMakePresets.json` uses schema 6.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
cmake -S . -B build/debug -DCMAKE_BUILD_TYPE=Debug
cmake --build build/debug --parallel 2
ctest --test-dir build/debug --output-on-failure
```

MPI requires a matching compiler/wrapper with `mpi_f08`. Live reference tests
use Cantera 3.2.0; optional CVODE uses the official double-precision static
Fortran interface from SUNDIALS 7.2.0. See [contributor setup](CONTRIBUTING.md)
for separate configurations. Do not reuse build caches across MPI/compiler
families.

## Run a first case

Run from the repository root:

```sh
./build/debug/pelef cases/sod/sod.nml
python3 tools/compare_sod.py --input sod.csv
./build/debug/pelef0d_h2o2_full cases/zero_d_h2o2_full/reactor.nml
python3 tools/check_zero_d_h2o2_full.py --input zero_d_h2o2_full.csv
```

Further examples live with their inputs:
[Sod](cases/sod/README.md),
[full H2/O2 reactor](cases/zero_d_h2o2_full/README.md),
[3D reacting flow](cases/reactive_hotspot_3d/README.md), and
[static 3D AMR](cases/amr_reactive_3d/README.md).

## Project navigation

- [Documentation index](docs/README.md): design, numerical methods, validation,
  historical milestones, and the source-tree map.
- [Completion roadmap](docs/completion_roadmap.md): remaining work and exit gates.
- [Contributing](CONTRIBUTING.md) and [agent rules](AGENTS.md): one-purpose
  changes, test evidence, and branch discipline.
- [Integration PR #102](https://github.com/GTT10/Fortranslate/pull/102) and
  [tracking issue #101](https://github.com/GTT10/Fortranslate/issues/101).

The former long README is preserved in the
[immutable 0.245.0 checkpoint](https://github.com/GTT10/Fortranslate/blob/509d07929d52d9e570a6af1f6ec7e5836bd4fe5d/README.md).
Its milestone statements describe their original scopes, not current overall
qualification. Keep new implementation history in design/validation records
rather than appending it to this entry page. License selection remains an
owner decision; this change does not choose a license.
