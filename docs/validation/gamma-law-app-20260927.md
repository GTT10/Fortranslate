# Calorically perfect HIT/SDC/LES gas, 2026-09-27

Implementation `ab7be8dc7a9eed46cd6f4bfc8f04b01f61e8c0c0`, based on
`f93e2cc5`, local `main`.

## Implemented model

The fixed regular executable `pelef_reactive_3d` now explicitly permits
`thermo_model='gamma_law'` with optional `gamma_law` namelist `gamma` (default
1.4). Other configuration readers retain their existing capability boundary;
configure-time selected mechanisms are not overwritten. The model is one
pseudo-species `AIR`, molecular weight 28.97 kg/kmol, with zero reactions.

The existing NASA7 thermodynamic interface exactly represents this calorically
perfect gas with the sole nonzero coefficient `cp/R=gamma/(gamma-1)` and zero
formation energy. Both temperature intervals use the same polynomial. The
supported range is 1 to 1e6 K, with a mathematically immaterial join at 1000 K.
The heat-capacity ratio must be finite, greater than one and representable with
positive cv. This is an analytic model constructor, not modified generated
mechanism data.

The formulas and AIR molecular weight match pinned PelePhysics
`6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`, inspected in
`Source/Eos/GammaLaw.H` and `Source/PhysicsConstants.H` with `git show`.
The model requires nonreacting SDC and no named composition. Enabled transport
requires prescribed constant coefficients and disabled species diffusion;
there is no molecular transport database for pseudo-air. The shared application
now permits an empty reaction array when chemistry is disabled. LES's redundant
two-species minimum was removed: one-species stress/heat equations and zero
composition gradients are valid. All other LES state checks remain.

This closes the regular HIT thermodynamic-model connection. It does not establish
external evolved-field parity, MPI/AMR gamma-law support or long-time turbulence
statistics. The new 16³ case template was not run in this increment.

## Completed checks

Existing `build/amr-cvode-app-debug`: GNU Fortran 13.3.0 Debug, MPI/spray/
SUNDIALS enabled, Open MPI 4.1.6 and SUNDIALS 7.2.0. No compiler/backend
configuration changed. Existing real-comparison, conversion, static-storage
and executable-stack warnings remain.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/amr-cvode-app-debug \
  --target pelef_reactive_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest \
  --test-dir build/amr-cvode-app-debug -R '^application_gamma_law$' --output-on-failure
python3 -m py_compile tests/application/check_gamma_law_app.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final CTest **1/1 passed**, 0 unexpected failures and 0 skips, 0.92 s
(total 0.97 s). The checker executes six successful applications and four
expected rejections. It uses periodic 4³ cells, a manufactured HIT velocity
cube, initial 300 K/101300 Pa, Mach 0.1, Taylor-based constant transport and SDC
to 1e-8 s. Gamma 1.4 is exercised with Smagorinsky, WALE, Vreman and Dynamic;
gamma 1.67 with WALE. Every final cell satisfies independent
`p=rho*R*T` and `rhoE-kinetic=p/(gamma-1)` formulas, with maximum relative
residual `5.7411274570087557e-16` against 3e-13. Pseudo-air mass fraction stays
exactly one. Derived transport coefficients match independent ideal-gas
cp/sound-speed formulas within 3e-14. Default gamma output is byte-identical
to explicit gamma 1.4. Every successful run passes the unchanged 3e-12 Euler
inventory gate.

Gamma=1, enabled chemistry, species diffusion and unconfigured molecular
transport are rejected without a final CSV. The first run exposed the
unnecessary two-species LES guard; its failure is preserved separately.
No tolerances or golden data changed.

Provenance, build, Python syntax, project contract and whitespace checks passed.
Pinned YAML SHA-256 remains
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
Artifacts: `build/gamma-law-app-checks/{build.log,ctest.log,LastTest.log,single-species-failure.log}`;
final per-case records are in
`build/amr-cvode-app-debug/gamma-law-app/gamma_4jb5fvug/`.
