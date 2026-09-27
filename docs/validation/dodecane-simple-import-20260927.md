# Dodecane Simple transport import — 2026-09-27

Implementation `6f0cfe5685bee6f806a28a7da7a9c902e2e4014d` extends the existing
source-verified Simple transport importer to Lu n-dodecane (53 species,
268 reactions). This is a prerequisite for pinned Spray-A-Wbreakup, whose
GNUmakefile selects `dodecane_lu` and `Simple`. Application input preparation
and execution of that case remain open.

## Provenance

PelePhysics `6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`, header
`Mechanisms/dodecane_lu/mechanism.H`:

- Git blob `ff4a703dc2a6e603e915fc5e0175b8df7df6f791`.
- SHA256 `5d943543c3f705867de594d1c8c9481ef7c4db1b22f403eeaeca515532cbc90c`.
- Local export bytes were compared directly with `git show` from the pinned
  commit in `/mnt/sn850/Projects/gtt2/external/PeleC/Submodules/PelePhysics`.

The header calls YAML species `CH2*` macro `CH2D_ID` at index 10;
`mechanism.cpp::CKSYMS_STR` confirms the chemical name. The reference manifest
records this one explicit alias. Header SHA/blob and species-order checks
remain required. No generated rates, thermochemistry or source YAML changed.

## Completed checks

Preparation used the existing Cantera environment:

```sh
build/cantera-venv/bin/python tools/prepare_spray_reference.py --fuel dodecane --output-dir build/dodecane-simple-reference --yaml build/pelec-pinned-hit-source/Submodules/PelePhysics/Mechanisms/dodecane_lu/mechanism.yaml --license build/pelec-pinned-hit-source/Submodules/PelePhysics/license.txt --simple-transport-header build/pelec-pinned-hit-source/Submodules/PelePhysics/Mechanisms/dodecane_lu/mechanism.H
python3 tools/generate_elementary_mechanism.py --input build/dodecane-simple-reference/dodecane_lu.json --output build/dodecane-simple-reference/generated.f90
python3 tools/check_mechanism_provenance.py --bundle build/dodecane-simple-reference/dodecane_lu.json
python3 build/check_dodecane_simple.py
```

The prepared bundle SHA256 is
`5aa088da8ed7da0d084fde0f5bb24b574fccceaba69a33fb91c6ab006f8913ee`.
The first attempt used system Python without Cantera and failed. The next
attempt exposed the previously unsupported `CH2*` macro alias. The final
preparation, regeneration and provenance checks pass.

GNU C++/Fortran 13.3.0, generic transport kernel from the existing Debug
`build/heptane-simple-mpi-debug/libpelef_core.a`. The C++ probe uses pinned
PelePhysics Simple.H and the dodecane header; compile definitions/includes
come from the configured pinned Spray-Evaporation build, with the dodecane
mechanism include prepended. Exact commands are retained in
`build/dodecane-simple-check/commands.log` and `result.json`.

- All 6,148 imported coefficients equal compiled C++ output exactly:
  53 viscosity fits, 53 conductivity fits, 1,431 binary-diffusion fits,
  each with four coefficients.
- Four 53-species states at 300, 900, 1,200 and 1,800 K compare viscosity,
  conductivity and all 53 rho-diffusion values. Maximum relative error
  `1.7286239934987665e-15`, below the existing transport comparison gate
  `5e-12`. States include pure N2, air, fuel/CO2/N2 and all species.
- Regenerated Fortran was compiled with `-g -fcheck=all`, fixed core modules
  and library. `build/dodecane-simple-check/generated_probe` compares all
  6,148 loaded coefficients exactly with C++ and rejects an undersized array.
  Log: `generated-probe.log`. Generated module/object output stays in that
  build directory (`-J build/dodecane-simple-check`).
- Species-order permutation and a changed header are rejected. The generic
  Fortran transport probe also rejects negative composition.

Python syntax for both modified tools and the probe script, fixed h2o2
provenance, `tools/check_project_contract.py --source .` and
`git diff --check` pass. No full selected application build, breakup case,
reacting CFD, physical validation or configuration matrix was run here.
