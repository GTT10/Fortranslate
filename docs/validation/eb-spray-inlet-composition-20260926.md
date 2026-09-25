# Independent EB spray inlet composition — 2026-09-26

Implementation 6f8c4c3, parent f37cf3d, local main, version 0.246.0.
Fixed/selected single-level spray applications now enable the shared per-face
composition inputs. Inlet construction resolves named mole fractions against
the active mechanism, converts to mass fractions and recovers temperature from
density and pressure. Each default face resets to initial gas composition.
Coupled checkpoint context includes all six counts and active names/fractions;
changed mixtures reject continuation. Particle injection remains separate.

## Completed focused evidence

GNU Fortran 13.3.0, Python 3.12.3, Ninja Debug, existing OpenMPI wrappers.
Fixed tests-enabled build and selected tests-disabled H2O2 full build remain
separate. No backend configuration changes.

Each final checker passed **3 successful executions and 13 expected rejections**,
with no skips. Continuous/save/resume gas, particle and film CSVs agree exactly,
as do printed cumulative inventories. Reflection and particle thermal/mass
transfer checks pass. Films are empty; this does not validate film evolution.
The new fixture overrides lower-y inlet with O2/N2 mole fractions .2/.8.
New rejections cover changed restart composition, unknown species and negative
fractions, in addition to existing boundary, immutable checkpoint, corrupted
payload, geometry/context and output-alias checks.

The initial selected build used a nonexistent target and stopped before
compilation; the corrected target below succeeded. The initial fixed checker
expected the word "context" in a rejection message, but the actual reader
reports "incompatible spray checkpoint". The assertion was corrected and the
complete fixed checker rerun successfully. No numerical corrections, tolerance
changes or golden-data updates were needed. Compiler warnings remain.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_eb_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_reactive_eb_3d_selected_spray --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tests/application/check_eb_spray.py --exe build/spray-amr-reactive-debug/pelef_eb_spray_3d --source . --work build/eb-spray-inlet-composition-checks/fixed --inflow --restart --inlet-composition
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tests/application/check_eb_spray.py --exe build/spray-amr-reactive-selected-debug/pelef_reactive_eb_3d_selected_spray --source . --work build/eb-spray-inlet-composition-checks/selected --inflow --restart --inlet-composition --selected
python3 -m py_compile tests/application/check_eb_spray.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Provenance, syntax, project contract and whitespace passed.
Mechanism SHA: 0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4.
Pinned PeleC: bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8.
Pinned PelePhysics: 6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b.

Artifacts: build/eb-spray-inlet-composition-checks/. Final runs:
fixed/eb_spray_ctc28zpl and selected/eb_spray_x31y04vq, with fixed-final.log,
selected.log, build.log and build-selected-final.log. Initial build/checker
failures are retained in build-selected.log and fixed.log.

## Limits

One inlet mixture and the H2O2 mechanism were exercised. Multiple simultaneous
mixtures, other mechanisms, MPI spray, full matrix, installed binaries,
external solver parity and physical nozzle/ignition validation were not rerun.
This is a coupled application input/restart check; complete PeleC parity remains
unproven. Domain-wall tree boundaries remain a separate implementation item.
