# Chemically activated rates — 2026-09-15

Source: parent `5a34567` plus this commit's activated-rate branches,
importer/generator, checkpoint handling, tests and documentation.
GNU Fortran 13.3 Debug in `build/port-spray-debug`; Cantera 3.2.0 oracle.

The pressure law follows
[Cantera 3.2 chemical activation](https://cantera.org/3.2/reference/kinetics/rate-constants.html#chemically-activated-reactions):
`k = k0 F/(1+Pr)`, where `Pr=k0[M]/kinfinity`. Lindemann, Troe and SRI
broadening are supported. The analytic concentration derivative includes both
the rational pressure term and the broadening derivative. Below the existing
reduced-pressure floor, the latter is frozen consistently with clipped F.
The normalized bundle distinguishes activated reaction types, and runtime
kinds 5/6 distinguish their law in packed MPI/checkpoint metadata. SRI records
are included for kind 6 as well as ordinary SRI kind 4.

## Focused evidence

- 10 portability tests passed, including activated Lindemann/SRI import and
  generation. The 12-test generator suite also passed.
- Four Fortran tests passed: ordinary kinetics, existing AMR checkpoint,
  SRI AMR checkpoint and activated-SRI AMR checkpoint. The latter verifies
  coefficient persistence and transactional rejection when D changes.
- Generated Fortran rates and analytic collider derivatives were compared
  with Cantera over 15 states for each activated broadening model:
  T=500/1000/2000 K, Pr=1e-6/1e-3/1/1e3/1e6. Maximum relative errors:

| Model | Rate | Concentration derivative |
| --- | ---: | ---: |
| Lindemann activated | 2.665e-15 | 1.019e-7 |
| Troe activated | 2.554e-15 | 1.250e-11 |
| SRI activated | 2.665e-15 | 8.710e-10 |
| Ordinary SRI regression | 2.665e-15 | 5.227e-10 |

The derivative oracle uses a fourth-order centered stencil with relative
interval 1e-3. The original second-order 1e-5 interval suffered cancellation
at the almost-flat Lindemann low-pressure limit (relative error 1.756e-6).
The oracle stencil was improved; the 1e-6 derivative tolerance was not relaxed.
Rate tolerance remains 1e-12. The synthetic reaction and coefficients are a
numerical fixture, not a physical validation dataset.

Artifacts under `build/activated-check/`:
- `sri-stkk20ps/`: activated Lindemann
- `sri-ncv8cpxc/`: activated Troe
- `sri-5qr6k1zm/`: activated SRI
- `sri-8mvxst1l/`: ordinary SRI regression

Each contains source YAML/JSON, generated Fortran, probe and `result.json`.
Activated-SRI source SHA256:
`9808bab362866e32bcc988c02d4733ea681fd7cefa48dc2162ca9614a45fee96`.
Source provenance, project contract, Python syntax and whitespace pass.

```sh
cmake --build build/port-spray-debug --target test_elementary_kinetics test_amr_reactive_3d_checkpoint --parallel 8
ctest --test-dir build/port-spray-debug --output-on-failure -R '^unit_(elementary_kinetics|amr_reactive_3d.*checkpoint)$'
python3 -m unittest discover -s tests/python -p test_mechanism_portability.py
python3 -m unittest discover -s tests/python -p test_generate_mechanism_bundle.py
build/cantera-venv/bin/python tools/check_sri_falloff.py --activated --model lindemann --build-dir build/port-spray-debug --work-dir build/activated-check
build/cantera-venv/bin/python tools/check_sri_falloff.py --activated --model troe --build-dir build/port-spray-debug --work-dir build/activated-check
build/cantera-venv/bin/python tools/check_sri_falloff.py --activated --model sri --build-dir build/port-spray-debug --work-dir build/activated-check
build/cantera-venv/bin/python tools/check_sri_falloff.py --model sri --build-dir build/port-spray-debug --work-dir build/activated-check
python3 tools/check_mechanism_provenance.py --bundle build/activated-check/sri-5qr6k1zm/sri.json --source build/activated-check/sri-5qr6k1zm/sri.yaml
python3 tools/check_project_contract.py --source .
python3 -m py_compile tools/generate_elementary_mechanism.py tools/ingest_cantera_mechanism.py tools/check_sri_falloff.py
git diff --check
```

No new activated-reaction MPI/EB run, ignition trajectory or physical benchmark
was performed. PLOG, Chebyshev and Tsang rates remain open.
