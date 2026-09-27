# Periodic SDC constant transport input, 2026-09-27

Implementation `b55b2250698e08361a23e5656da2063791afb5be`, based on
`39d7d67c`, local `main`.

## Integration

Regular fixed/selected gas applications now read optional
`reactive_constant_transport` with `enabled`, `viscosity` and `conductivity`.
Enabled input requires finite nonnegative SI coefficients, periodic boundaries,
SDC and enabled transport. The allocated coefficient pair is passed to both
transport timestep evaluation and old/new SDC transport fluxes. Omitting it
retains mixture transport. Optional LES rates continue to contribute alongside
the configured transport source.

The regular gas SDC and LES SDC library APIs additionally accept the pair.
The core rejects coefficients without transport, with caller-supplied halos,
or with nonperiodic boundaries rather than silently using molecular values.
Species diffusion remains the mixture model and can be disabled through the
existing flags. No generated database or checkpoint record changed.

See [input example](../../cases/reactive_hit_3d/README.md). This completes the
periodic regular application connection of the preceding kernel API. Automatic
Taylor-Reynolds coefficient construction, nonperiodic/AMR support and full HIT
reference parity remain separate. The standalone LES library forwarding is
compiled; the application test uses the shared application diffusion callback.

## Completed checks

Existing `build/amr-cvode-app-debug`, GNU Fortran 13.3.0 Debug, Open MPI 4.1.6,
SUNDIALS 7.2.0, MPI/spray/SUNDIALS enabled and selected full H2/O2. No compiler
or backend configuration changed. Existing real-comparison, conversion and
executable-stack warnings remain.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/amr-cvode-app-debug \
  --target pelef_reactive_3d pelef_reactive_3d_selected --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest \
  --test-dir build/amr-cvode-app-debug -R '^application_constant_transport$' --output-on-failure
python3 -m py_compile tests/application/check_hit_input_app.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final CTest **1/1 passed, 0 failures, 0 skips**, 3.57 s (total 3.62 s).
The checker completes seven successful solver runs and four expected
rejections. A manufactured HIT field supplies nonuniform velocities on 8×4×4
cells. The constant values are 0.001 Pa s and 0.5 W/(m K).
ASCII/binary/fixed/selected initialization CSVs are identical; periodic SDC to
1e-8 s with WALE gives exact fixed/selected CSVs. Constant transport without
LES differs from the zero-viscosity/conductivity result, establishing active
application routing. Maximum normalized Euler inventory residual across
successful runs is `2.3414515441381875e-15`, below `3e-12`.
Negative coefficients and nonperiodic constant transport reject without a
final CSV; the reused checker also rejects truncated HIT binary data and
nonascending coordinates. Timestep-limit dominance and reacting constant
transport were not separately tested.

Initial compilation failed on a line exceeding Fortran's 132-column limit.
An accidentally subsequent test invocation used old binaries and failed the
active-coefficient assertion; it is preserved as a stale-binary failure, not
counted as validation. After line continuation repair, both targets rebuilt
successfully before the final passing CTest. No tolerance or golden data changed.

Provenance, builds, Python syntax, project contract and whitespace checks pass.
Pinned YAML SHA-256 remains
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
Logs are in `build/constant-transport-app-checks/`, including `LastTest.log` and
initial failure logs. Final inputs/CSVs/per-case logs are under
`build/amr-cvode-app-debug/constant-transport-app/hit_e6wb212c/`.
