# Pinned decane spray mechanism — 2026-09-27

Implementation: `88aec3c5`. This supplies the reaction-free NC10H22/O2/N2
mechanism required by pinned PeleC's Spray-Jet case. It does not implement or
validate that case's injection trajectory.

`references/decane_3sp.json` pins the YAML and generated transport header to
PelePhysics `6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`. Both were checked against
that commit's local Git objects and the existing pinned source export.
YAML SHA-256: `cd27f0bc24ad9ebd31464022a7ac0fcb54501cb979cbcb6a2cb1c51468f7f9d1`.
The upstream source and license bytes are preserved; no rates or coefficients
were hand-edited. The importer emits three species, zero reactions and optional
Simple fits. Regeneration now selects either the heptane or decane manifest by
the complete pinned header identity, then checks species ordering and bytes.
CMake tracks both manifest dependencies.

## Reproduction

With Cantera 3.2.0 in the existing environment:

```sh
build/cantera-venv/bin/python tools/prepare_spray_reference.py --fuel decane-3sp --output-dir build/decane-3sp-simple-reference --yaml build/pelec-pinned-hit-source/Submodules/PelePhysics/Mechanisms/decane_3sp/mechanism.yaml --license build/pelec-pinned-hit-source/Submodules/PelePhysics/license.txt --simple-transport-header build/pelec-pinned-hit-source/Submodules/PelePhysics/Mechanisms/decane_3sp/mechanism.H
python3 tools/check_mechanism_provenance.py --bundle build/decane-3sp-simple-reference/decane_3sp.json
python3 tools/generate_elementary_mechanism.py --input build/decane-3sp-simple-reference/decane_3sp.json --output build/decane-3sp-simple-reference/decane_3sp.F90
python3 tools/check_simple_transport_reference.py --reference-build build/pelec-pinned-spray --reference-case Spray-Jet --candidate-build build/heptane-simple-mpi-debug --bundle build/decane-3sp-simple-reference/decane_3sp.json --bulk --output build/decane-simple-transport-reference
python3 -m unittest discover -s tests/python -p test_mechanism_provenance.py
```

The compiled checker uses the configured Spray-Jet C++ include flags and the
actual pinned decane header; it does not require a built Spray-Jet application.
Its Fortran harness compiles the newly generated decane module and links the
existing mechanism-independent transport libraries. The candidate build's
heptane name does not select the harness species.

GNU C++/Fortran 13.3.0: six states / 36 transport coefficients pass the existing
`2e-12` gate, maximum relative difference `1.1242436427969756e-15`. The harness
also passes invalid-composition rejection, bulk compression/viscous-work and
timestep checks, three face-coefficient comparisons and fitted checkpoint
record checks. Exact commands, pinned identities, source hashes and results
are in `build/decane-simple-transport-reference/`.

Fourteen existing provenance tests pass. A focused negative audit rejects a
heptane filename under the decane header identity and an NC7H16 species under
the decane header (`build/decane-simple-import-check2/result.json`). Its first
attempt tried to execute a compiler diagnostic from a mixed commands/output
log; that orchestration failure remains in `build/decane-simple-import-check/`
and is not counted as a pass. The reusable checker above succeeds directly.
Heptane generated source remains byte-identical to the preceding successful
face-check source. Modified Python syntax, project contract and diff checks pass.

No decane application build/run, gas-particle parity, nozzle validation or full
configuration matrix is claimed. Selected applications must keep chemistry OFF
for this zero-reaction mechanism and stay inside its source thermo ranges.
