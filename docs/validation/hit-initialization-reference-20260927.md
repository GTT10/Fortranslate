# Pinned HIT initialization comparison — 2026-09-27

Checker source: `9ae3a59a6df1f86f540998d3db7400b407e55744`.
Fortran implementation: `90085c9b2908540b24a4ca25e7981b9fe6e2e2ea`
(documentation-only successor `e2866396`).

`tools/check_hit_initialization_reference.py` retrieves SHA-checked blobs from
PeleC `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8` and PelePhysics
`6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b` with `git show`, without changing
either checkout. It compiles the unchanged `pc_initdata` and `locate` functions
with minimal array/geometry adapters. The adapter supplies gamma-law CGS
thermodynamic properties and converts output density, momentum and energy to SI.
This does not compile the full PeleC application or its EOS implementation.

Both sides read the pinned 32³ velocity file. Comparison uses 4³ cells in a
2π cm cube, 300 K, 101300 Pa, Mach 0.1, gamma 1.4 and molecular weight
28.97 kg/kmol. Fortran uses coordinate factor 0.01 and a tiny SDC step of
`1e-20` s, with transport, chemistry and LES disabled. Cell coordinates and
output time are checked to prevent accidental row/trajectory mismatches.

## Completed checks

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_hit_initialization_reference.py --reference /mnt/sn850/Projects/gtt2/external/PeleC --exe build/amr-cvode-app-debug/pelef_reactive_3d --work build/hit-initialization-reference
python3 -m py_compile tools/check_hit_initialization_reference.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

All completed checks passed. One final reference comparison, 64 cells, zero
failures/skips. GNU C++ 13.3.0, C++17, `-O0`; existing GNU Fortran 13.3.0 Debug
application with MPI/spray/SUNDIALS enabled. No compiler/backend changes or
Fortran rebuild were needed. No CTest matrix was rerun.

| Field | Maximum normalized absolute difference |
| --- | ---: |
| Density | 1.8872860278019619e-16 |
| x momentum | 3.766091378023539e-15 |
| y momentum | 3.464804067781656e-15 |
| z momentum | 3.766091378023539e-15 |
| Total energy density | 2.2984268870075976e-16 |
| Temperature | 1.8947806286936006e-16 |

Gate: `1e-12`. Normalizers are initial density, density times target RMS
velocity, `p/(gamma-1)`, and initial temperature, respectively. These fixed
scales avoid relative errors at nearly zero local momentum.

Final artifacts: `build/hit-initialization-reference/hit_reference_iw14x_5r/`:
source blobs, C++ adapter/binary, inputs, CSVs, logs and `summary.json` with
source and Fortran-binary SHA-256 hashes. An earlier invocation used an absent
`app/pelef_reactive_3d` executable path and stopped before running Fortran; it
is excluded. The first successful comparison preceded added row/time checks;
the final result above includes those checks.

This establishes initialization-kernel agreement. It does not establish
evolved-field parity, forcing statistics, AMReX execution parity, or physical
turbulence validation.
