# Reaction-free spray reference import — 2026-09-27

Implementation: `048c81201b8fef5a9ec723552e8f83fbb3ab2f4e`.
The pinned Spray-Evaporation executable uses PelePhysics `heptane_3sp`:
NC7H16/O2/N2, three species, zero reactions. The previous importer rejected
both its `kinetics_model='none'` phase and its empty reaction list.

`--nonreacting` now explicitly enables reaction-free import. It requires zero
source reactions and records `nonreacting: true` in the normalized bundle.
Default import still rejects an empty reaction mechanism; marking a reacting
source nonreacting also fails. The generator requires source provenance,
thermodynamic data and transport data for this mode and emits typed empty
integer/logical metadata arrays. No dummy reaction is synthesized.

The pinned reference manifest is `references/heptane_3sp.json`, from
PelePhysics `6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`.
YAML SHA-256: `5df1e27a9551d3c54e17945e12e6e1e6db366f90df01ea548b0f6f19f30284a4`.
The upstream description mentions decane, but the species composition is
NC7H16; the source bytes and chemical identity were preserved.

## Completed checks

Pinned Cantera 3.2.0 through `build/cantera-venv/bin/python`:

```sh
build/cantera-venv/bin/python tests/unit/test_ingest_cantera_mechanism.py PinnedCanteraBundleTest.test_explicit_nonreacting_bundle PinnedCanteraBundleTest.test_bundle_is_complete_ordered_and_generator_inputtable
build/cantera-venv/bin/python tools/prepare_spray_reference.py --fuel heptane-3sp --output-dir build/heptane-3sp-reference --yaml build/pelec-pinned-hit-source/Submodules/PelePhysics/Mechanisms/heptane_3sp/mechanism.yaml --license build/pelec-pinned-hit-source/Submodules/PelePhysics/license.txt
python3 tools/check_mechanism_provenance.py --bundle build/heptane-3sp-reference/heptane_3sp.json
```

Two tests passed, zero failures/skips (0.140 s). They check explicit opt-in,
source hash preservation, typed empty arrays, default rejection, rejection of
a reacting source marked nonreacting, and ordinary H2/O2 bundle generation.
Reference preparation verified YAML SHA/Git blob and license SHA and produced
three species and zero reactions.

An earlier direct import produced `build/heptane-nonreacting/heptane_3sp.json`.
The generator emitted `heptane_3sp.F90` there; GNU Fortran 13.3.0 compiled it
with existing `build/amr-cvode-app-debug/modules` and an isolated module/object
output directory. This was a generated-module compilation, not a linked
application run. Python syntax for the three tools and touched test, the
existing H2/O2 provenance preflight, project contract and diff checks passed.

No generated source was hand-edited. The selected application and full
Spray-Evaporation comparison remain pending. Upstream's example starts gas at
273 K, below the imported NC7H16 NASA7 lower bound of 300 K. That range has not
been widened; usable comparison conditions need an explicit temperature
decision. The supplied example is also two-dimensional with legacy `Hard`
boundary labels, so it cannot be presented as a verbatim current 3D case.
