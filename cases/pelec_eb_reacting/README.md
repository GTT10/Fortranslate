# Matched PeleC planar-EB reacting trajectory

This case compares the unchanged pinned PeleC `EB-C3` problem against the
selected PeleF EB-tree gas application. Both use the pinned PelePhysics
LiDryer mechanism (9 species, 21 reactions), 2000 K, 101325 Pa, and H2/O2/N2
mass fractions 0.06/0.50/0.44. The domain is a 1 cm cube with an 8-cubed
single-level mesh and a planar wall at z = 0.31 cm. Fluid is above the wall.
There are 384 exposed cells, including 64 cut cells with volume fraction 0.52.

Both applications use MOL with two iterations and CVODE, no molecular
transport, and three fixed 100 ns steps. The fields are initially uniform,
so this exercises the uniform-reacting limit on an EB mesh. It does not
establish gradient, coarse/fine, redistribution, regrid or ignition parity.
`pelec.inp` uses CGS; `pelef.nml` uses SI. The latter's density and mole
fractions are derived from the same mass fractions and molecular weights.

## Build

The source YAML is byte-pinned to PelePhysics
`6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b/Mechanisms/LiDryer/mechanism.yaml`.
The normalized bundle is generated with the repository importer and Cantera
3.2.0; do not hand-edit the bundle or generated rates.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/lidryer.json
cmake -S . -B build/lidryer -DCMAKE_BUILD_TYPE=Debug \
  -DPELEF_ENABLE_SUNDIALS=ON \
  -DCMAKE_PREFIX_PATH=/path/to/sundials-7.2.0-fortran-static \
  -DPELEF_MECHANISM_BUNDLE="$PWD/mechanisms/lidryer.json" \
  -DPELEF_MECHANISM_SOURCE="$PWD/mechanisms/lidryer.yaml"
cmake --build build/lidryer --target pelef_reactive_eb_tree_3d_selected --parallel 4
```

Build `PeleC-EB-C3` and the existing `tools/plotfile_export` exporter against
the source/dependencies pinned in `references/pelec_baseline.json`. The
checker expects the prepared source archive's
`pelef-reference-provenance.json` to match that manifest. It independently
verifies the source hashes in `references/pelec_eb_reacting.json` before
running either solver; its report records executable hashes as well.

```sh
python3 tools/check_eb_reacting_reference.py \
  --pelec /path/to/PeleC-EB-C3 \
  --fortran build/lidryer/pelef_reactive_eb_tree_3d_selected \
  --exporter /path/to/pelef_plotfile_csv \
  --reference-source /path/to/prepared-pinned-PeleC \
  --work build/eb-reacting-comparison
```

## Comparison contract

`acceptance.json` owns the physical conditions and tolerances. The checker
compares all four saved physical times, matches exposed cells and fluid
volumes, converts CGS fields to SI, and checks density, three momenta, total
and retained internal energy, temperature and all nine species. It also
integrates the Fortran CSVs independently against cumulative budgets.
Covered reference cells are excluded by their zero volume fraction before
checking derived-field finiteness; undefined covered-cell quantities are
not fluid states.

Hydrogen consumption must exceed `1e-6` of its initial inventory and agree
within 0.1% of the consumed amount, in addition to field tolerances. This
prevents a nearly inert trajectory from passing only because its changes
fall below an absolute species tolerance. Reactor error controls differ
because PeleC integrates temperature/species with component-scaled
tolerances while PeleF uses energy-constrained reduced species; comparison
tolerances remain separate from those integration controls.

Every execution creates a unique directory with logs, native plotfiles,
Fortran snapshots, converted reference fields, hashes and `checks.json`.
`--reuse-run /path/to/completed-run` copies existing completed simulations
after verifying their input, executable, mechanism and contract identities,
then reruns only conversion/comparison. It preserves the original failed
report and avoids repeating integration when repairing a checker.
