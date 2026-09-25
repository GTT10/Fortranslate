# Serial EB spray periodic application — 2026-09-25

## Source and scope

Implementation f3c95dc, parent 4ec6211, local main. The shared fixed/selected
serial EB spray application now accepts paired periodic faces with CTU/pelec.
Injection and restart particle checks use the periodic field stencil. Existing
boundary checkpoint identity and wrapped extended geometry are reused.
Explicit boundaries still require particle_outflow=true, which enables outer
face processing; periodic faces wrap particles without removing inventory.
Periodic impact/film options reject at startup. Specular embedded-wall drift
remains supported.

## Completed checks

GNU Fortran 13.3.0, Ninja Debug, Python 3.12.3. Existing fixed tests-enabled
build/spray-amr-reactive-debug and selected tests-disabled
build/spray-amr-reactive-selected-debug configurations; selected bundle
mechanisms/h2o2_full.json. No compiler/backend changes.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_eb_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_reactive_eb_3d_selected_spray --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tests/application/check_eb_spray.py --exe build/spray-amr-reactive-debug/pelef_eb_spray_3d --source . --work build/eb-spray-periodic-application-checks/fixed --periodic --restart
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tests/application/check_eb_spray.py --exe build/spray-amr-reactive-selected-debug/pelef_reactive_eb_3d_selected_spray --source . --work build/eb-spray-periodic-application-checks/selected --periodic --restart --selected
python3 -m py_compile tests/application/check_eb_spray.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Each application checker passed **3 successful executions and 10 expected
rejections**, with no unexpected failures or skips. A parcel pair injected
1 nm below the high-y face at 20 m/s wraps to low y with positive mass and
heated liquid. Coupled chemistry and transport are enabled. Continuous and
resumed gas, particle and film CSVs, plus cumulative printed inventory ledgers,
match byte-for-byte. The film arrays are empty in this fixture.

Expected rejections cover changed hydro, changed periodic pairs, an unpaired
face, a periodic direction incompatible with the EB geometry, periodic impact
films, checkpoint overwrite, changed spray context, changed geometry context,
truncated payload and output alias. The boundary-change rejection checks also
verify existing CSVs are untouched.

Provenance, syntax, project contract and whitespace checks passed. Mechanism
SHA remains 0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4.
No tolerances or golden data changed. Existing compiler warnings remain.

Artifacts under build/eb-spray-periodic-application-checks/:

- build-fixed.log and build-selected.log.
- fixed/eb_spray_omy5dgg8: input, output, restart and rejection logs.
- selected/eb_spray_uojr5kbu: input, output, restart and rejection logs.

## Limits

This qualifies bounded serial application dispatch and continuation. MPI/tree
periodic particle boundaries, periodic impact/film coupling and domain-wall
particle application dispatch remain open. No full configuration matrix,
installed audit, external solver comparison or physical spray validation was
run. Full PeleC parity remains unproven.
