# Checkpoint source-integrity repair (2026-09-08)

## Scope and cause

The integration checkpoint `509d07929d52d9e570a6af1f6ec7e5836bd4fe5d`
failed Configure because its committed `mechanisms/h2o2_cantera.yaml`
was not byte-identical to the source recorded in `h2o2_full.json`.
The difference is exactly two removed spaces on line 4, inside the YAML
`description` block. No reaction, thermochemistry, or transport data differ.

| Source | Bytes | SHA-256 |
| --- | ---: | --- |
| Checkpoint YAML | 10031 | `fa5076030f7ccea6c0b42d034a906d8c0c942a3b0f5cd5ba4cf3f197976e6974` |
| Cantera 3.2.0 packaged `h2o2.yaml` | 10033 | `0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4` |

The repair restores the packaged source bytes, including those two spaces.
It does not change the bundle hash, generated Fortran, rate constants, solver,
acceptance tolerances, or Configure-time provenance enforcement.
The path-specific Git/editor settings preserve this intentional upstream
whitespace; do not run trailing-whitespace cleanup on this vendored source.

## Independent diagnostic evidence

[Mechanism integrity run 34207681798](https://github.com/GTT10/Fortranslate/actions/runs/34207681798)
ran on commit `47d99fdec3077d54b4192414de3f173342d7913a` with
Python 3.12 and Cantera 3.2.0. All 14 new dependency-free checker tests passed.
Regenerating the complete bundle from the checkpoint YAML changed **only**
`source.sha256`; all other normalized fields matched exactly. The packaged
Cantera YAML had the originally recorded hash. This is why restoring the
original source is justified rather than weakening the check.

The run deliberately failed its provenance gate and retained
[artifact 10048438553](https://github.com/GTT10/Fortranslate/actions/runs/34207681798/artifacts/10048438553):
source commit/tree, repository archive, dependency versions, packaged source,
regenerated bundle and unified diff. Artifact ZIP SHA-256:
`8fae0f85d059d6213275360fd8f0f7cc42aec509dcbe34654b3df47fb63f7196`.
Artifacts expire; the source identities and diagnosis above are retained here.

## Reproduce

```sh
python3 -m unittest discover -s tests/python -p test_mechanism_provenance.py -v
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
# With the pinned optional Cantera 3.2.0 dependency installed:
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json --regenerate
```

The byte-only check requires no Cantera. `--regenerate` compares the full
normalized bundle, including runtime provenance, and does not update inputs.
`--output /path/to/candidate.json` saves a separate review candidate, not an
in-place repair. Any mismatch still returns a nonzero exit status.

The source repair alone is not numerical qualification. The historical
4570-test matrix in the handoff and any new CI results must be reported with
their own tested source/configuration identities; neither replaces the other.
