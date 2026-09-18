# MPI EB component property controls — 2026-09-18

Implementation `af0192f`, base `02d742f`.

The MPI EB namelist exposes existing component density/heat-capacity polynomial,
Antoine pressure and Clausius surface-cap fields. Components are initialized
against actual gas species and validated; unused entries must retain defaults.
Port templates and injection inventories use evaluated mixture density, heat
capacity and energy in a constant temporary liquid, preventing component-one
polynomials from being applied twice. All controls enter restart context.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6 (`/usr/bin/mpifort`, `/usr/bin/mpiexec`),
Debug/Ninja; reused fixed tests-enabled and selected H2O2 tests-disabled builds.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_eb_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tests/application/check_mpi_eb_nozzle.py --exe build/spray-amr-reactive-debug/pelef_mpi_eb_spray_3d --mpiexec /usr/bin/mpiexec --source . --work build/spray-mpi-eb-properties-checks/runs --rate --properties
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_reactive_eb_tree_3d_selected_spray --parallel 2
python3 -m py_compile tests/application/check_mpi_eb_nozzle.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

The final driver passed all five launches: full two/four-rank, two-rank stop,
four-rank resume and intentional invalid-density rejection. No final failures
or skips. Both executable builds, provenance, syntax, contract and whitespace
passed. Selected execution was not run. The driver is registered as
`spray_mpi_eb_properties_application`; this run called it directly, not via CTest.

Two earlier driver attempts failed. The first fixture requested H2O2, absent
from the fixed elementary application mechanism; initialization correctly
rejected it. A debugger inspection identified component two as the failure.
The fixture now uses supported H2 alongside H2O, with synthetic property fits.
The next attempt passed all positive calculations but expected the negative
case to fail at component validation; it actually rejected at injected-mixture
validation. Correcting that expected diagnostic completed the final driver.
No production tolerance or property-validity rule was loosened.

At 350 K, the component density polynomials give 1030 and 865 kg/m3. Two ports
use 0.25/0.75 and 0.75/0.25 mass fractions. Fresh-particle size reconstruction
and profile-derived speed use independently calculated harmonic mixture density.
Component fractions agree within 1e-13, speed within relative 1e-12, event dose
within relative 1e-13. Both components have polynomial heat capacity; water uses
Antoine pressure. Earlier parcels exchange heat/mass with the reacting gas.
The application enforces its composite mass/species/momentum/energy budget.
Gas and particle CSVs agree across ranks and after restart to relative 3e-12,
scale floor 1e-30. A negative density fit rejects before output.

These are synthetic interface/conservation checks, not calibrated liquid-H2
property data or physical mixture-boiling validation. Surface-cap behavior,
wall impact/film behavior with these fits, and other mechanisms were not newly
executed. No broad matrix, installed audit, scaling or external parity was run.

Artifacts: `build/spray-mpi-eb-properties-checks/` contains build/provenance/
contract logs and `application.log`, `application-final.log` (earlier failures),
`application-verified.log` (pass). Final inputs, outputs and checkpoints are under
`runs/nozzle_lmfgyel9`; earlier attempts remain under `runs/nozzle_bu7jufqp` and
`runs/nozzle_5tab4y_n`.
