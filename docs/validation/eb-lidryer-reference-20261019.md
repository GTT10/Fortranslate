# Uniform reacting EB trajectory against pinned PeleC, 2026-10-19

The new [matched case](../../cases/pelec_eb_reacting/README.md) connects the
PeleF intermediate-output path to an independent pinned PeleC application.
Both use LiDryer (9 species, 21 reactions), rather than comparing the previous
Cantera H2/O2 mechanism against a different reference chemistry. The YAML
is copied byte-for-byte from the pinned PelePhysics source and imported by
the existing Cantera 3.2.0 importer. Regeneration is byte-identical.
No numerical solver implementation or acceptance tolerance was loosened.

## Completed comparison

The accepted case has 2000 K, 101325 Pa, H2/O2/N2 mass fractions
0.06/0.50/0.44, an 8-cubed mesh in a 1 cm cube, and a wall at z = 0.31 cm.
All 384 exposed cells are compared, including 64 cut cells at volume
fraction 0.52. Four samples span 0, 100, 200 and 300 ns. Both solvers use
two-iteration MOL and CVODE, with no molecular transport. The uniform
initial state makes this a uniform-reacting limit on an EB mesh; gradient,
redistribution, AMR/reflux, regrid and physical ignition validation remain open.

| Quantity | Maximum observed difference | Predeclared acceptance |
| --- | ---: | ---: |
| Temperature | 7.7761e-8 K | 1e-3 K |
| Species mass fraction, all nine species | 4.8162e-10 | 1e-7 |
| Total energy density, relative to field maximum | 1.2329e-16 | 1e-8 |
| Hydrogen consumption, relative to consumed amount | 2.2235e-5 (0.0022235%) | 1e-3 (0.1%) |

Reference and Fortran hydrogen consumption fractions were respectively
`4.4788038653e-5` and `4.4787042806e-5` of initial H2. Both exceed the
`1e-6` signal gate. Density, all momenta, retained internal energy, initial
conditions, physical times, fluid volumes and independently integrated
Fortran budgets also pass. The checker converts CGS reference quantities to
SI before comparison. It excludes covered reference cells by volume fraction
before checking fluid-state finiteness.

The accepted physical executions took 2.62 s for PeleC and 16.01 s for
PeleF. These are small-case wall times with different Release/Debug builds,
not a production speed comparison. The final condition-contract audit reused
those completed outputs after verifying executable, input, mechanism and
contract hashes; neither solver was rerun for that audit.

## Identity and reproduction

Dirty-main base: `0a0e00ebe9138145fbb0e0e745665d0d768df009`.
The source reference remains PeleC
`bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8` with the dependencies in
`references/pelec_baseline.json`. Ninety-one PeleC source/problem/CMake files
were compared against a freshly downloaded archive of that commit and
matched. Three pinned PelePhysics mechanism/reactor files also matched
their recorded hashes. The 94-file contract is
`references/pelec_eb_reacting.json`; the complete dependency build matrix
was not independently requalified here. Source archives lack Git metadata,
so native startup reports `unknown` Git hashes; archive/source and binary
hashes are the identity evidence, not those startup strings.

Artifacts are under `build/eb-lidryer-reference-20261019/`. `result.json`
points to the accepted comparison and original physical execution reports.
`runtime-source-sha256.json`, `source-sha256.json`, `source/`, input copies,
binary hashes, native plotfiles, CSVs, source-audit manifest and build logs
retain the exact evidence. The original accepted physical run is
`runs/matched_eb_ts8rbqqs`; the final checker/condition audit is
`runs/matched_eb_hsxvfotu`.

The existing GNU 13.3.0 configurations were reused without compiler/MPI/backend
changes. PeleC EB-C3 was built in the Release reference directory. PeleF's
Debug MPI/spray/SUNDIALS configuration added the LiDryer selected bundle and
built only `pelef_reactive_eb_tree_3d_selected`; shared numerical sources
were checked unchanged from the preceding snapshot build. Both use SUNDIALS
7.2.0. The existing compiler warnings are retained in logs.

```sh
env -i PATH=/usr/bin:/bin LC_ALL=C \
  python3 tools/check_eb_reacting_reference.py \
  --pelec build/pelec-pinned-spray/Exec/RegTests/EB-C3/PeleC-EB-C3 \
  --fortran build/reacting-hot-diagnosis-20261019/registration/pelef_reactive_eb_tree_3d_selected \
  --exporter build/pelec-plotfile-export/pelef_plotfile_csv \
  --reference-source build/pelec-pinned-hit-source \
  --work build/eb-lidryer-reference-20261019/runs
```

This was one accepted direct comparison, not a CTest/configuration matrix.
Mechanism provenance/regeneration, Python syntax, project contract and
whitespace checks passed. No installed-binary audit was run.

## Retained unsuccessful preparation attempts

PeleC first refused a 4-cubed input because EB requires a blocking factor
of at least 8; both inputs were resized to 8 cubed. It then refused the
abbreviated CVODE selector; the documented `denseAJ_direct` selector was
used. At 1500 K both solvers completed, but an initial checker rejected
undefined derived quantities in zero-volume covered cells. Recomparison
corrected the fluid-cell selection without repeating either solver. The
1500 K trajectory met the field limits but failed the minimum reaction signal
(`2.65e-7` H2 consumption fraction). Both inputs were raised to 2000 K,
retaining the field/signal limits and adding the consumed-amount relative
gate before that run. Every failed attempt remains in its own directory.
