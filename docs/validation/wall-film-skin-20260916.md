# Wall-film skin interpolation correction — 2026-09-16

## Source and correction

Tested parent: `84bcded024993b3988b9f2df9b73c040b3bf4af3` plus this commit's
source/tests. Exact file SHA-256 values and build/test logs are retained under
`build/wall-film-skin-checks/`.

Pinned PelePhysics `6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`
`Source/Spray/BreakupSplash/WallFilm.H` declares `rule=1./2.` and passes it to
`calcVaporState` in `Source/Spray/Drag.H`. Its nearby comment says one-third,
but the executed temperature and composition interpolation uses one-half.
The source hashes and license are recorded in `references/spray_wall_impact.json`.

The prior implementation incorrectly used one-third for the wall-film skin.
This corrects both temperature and vapor mass fraction to one-half, retaining
the nonnegative gas-minus-film temperature difference and normalized carrier
composition. Droplet interpolation retains its existing expression.
For a 350 K liquid film and 500 K gas, the wall skin is now 425 K rather than
400 K. Gas transport and density therefore change as well.

The new `wall-film-evolution-v2-half-skin` checkpoint context rejects evolving
film restarts from the previous rule. Static-film and droplet-only contexts
are unchanged. This supersedes the one-third-rule description in the earlier
wall-film evolution evidence; that dated evidence remains a record of the
previous implementation.

## Focused verification

GNU Fortran 13.3.0, existing Debug spray build, MPI disabled:

```sh
cmake --build build/port-spray-debug --target pelef_spray_3d test_spray_core --parallel 8
ctest --test-dir build/port-spray-debug -R '^spray_(film_binary|application_cli)$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 -m py_compile tests/spray/test_spray_application.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Two CTest entries passed in 26.38 s, including 22 CLI tests, with no failed or
skipped tests. Provenance, syntax and contract checks passed. An initial build
failed because the new test omitted the gas-constant helper's required `ok`
argument; this was corrected before the successful build and tests.

The focused kernel checks independently evaluate binary diffusion at 425 K,
ideal-gas density with the midpoint composition, and the cold-gas temperature
clamp. The existing droplet binary limit remains checked at 400 K. CLI tests
include exact split-run CSV/history agreement and rejection of a checkpoint
whose context is changed back to the previous film rule.

No full-suite, external-solver parity or physical qualification was performed.
Mixture films, lateral transport and the previously documented film-model
limitations remain open.
