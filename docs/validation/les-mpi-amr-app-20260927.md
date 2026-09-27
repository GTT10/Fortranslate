# MPI gas application LES SDC — 2026-09-27

Implementation `614cfc16`, based on `e18d0c7b`, local main, version 0.246.0.

## Change

With `PELEF_ENABLE_SPRAY=ON`, the MPI fixed/selected gas frontends route the
optional `amr_tree_les` namelist through the LES runtime. Enabling LES requires
SDC. Smagorinsky, WALE, Vreman and Dynamic options feed the instantaneous MPI
LES evaluator; common and recursive SDC retain their existing transport,
chemistry, forcing and synchronization paths. An absent/disabled LES group uses
the original gas path. The generic driver receives a diffusion callback rather
than depending on the optional LES runtime.

Before a step, the driver evaluates diffusion and limits the timestep by 0.1
using the largest diffusivity and the current smallest patch spacing. This is
conservative even for recursive stepping. Each stage receives current topology,
owners, boundary reservoirs and parent temporal endpoints. The canonical LES
configuration is appended to checkpoint context, including its required final
newline; changed or removed LES cannot continue the same checkpoint.

## Independently completed checks

Artifacts: `build/les-mpi-amr-app-checks/`, plus per-run input/log/output files in
`build/spray-amr-reactive-debug/les-mpi-amr-sdc-app-{common,recursive}/`.
Existing fixed-mechanism Ninja Debug build, GNU Fortran 13.3.0, Open MPI 4.1.6
(`/usr/bin/mpiexec`).

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_amr_tree_reactive_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^les_mpi_amr_sdc_app_(common|recursive)$' --output-on-failure
python3 -m py_compile tools/check_les_mpi_amr_sdc_app.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final **2/2 passed**, no failures/skips: common 19.07 s, recursive 18.97 s,
CTest total 38.09 s (`ctest-final.log`, `LastTest-final.log`). Build, MPI runtime
ABI audit, mechanism provenance, Python syntax, project contract and whitespace
checks passed. Mechanism source SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
Existing real-comparison and executable-stack linker warnings remain.

Both cases use a two-level 4-cubed-root hotspot with moving no-slip walls,
molecular transport, chemistry disabled and Smagorinsky LES. LES changes the
final gas CSV compared with disabled LES. Common stepping takes 8 steps;
recursive stepping takes 4. Saving after step 1 on 2 ranks and continuing on
4 ranks reproduces both final gas and budget CSVs byte-for-byte. Changing `cs`
or removing the LES group rejects restart due to context mismatch.

Development failures were corrected: the first configure placed target wiring
before target creation; the first runtime test passed an unassociated procedure
pointer; the next exposed the input-group whitelist; checkpoint writing then
exposed the missing final context newline. A long argument declaration also
failed compilation before being wrapped. Retained initial/context failure logs
are separate from the final passing evidence. Open MPI emitted display
 authorization messages in some captured subprocess logs without preventing the
final checks.

## Limits

Selected frontend wiring is present but was not built/run in this configuration.
The focused application checks cover Smagorinsky only; the preceding library
checks cover all four closures. Active chemistry, regridding, inflow/outlet
variants, dynamic application settings, installed binaries, external PeleC parity
and physical LES validation were not run. Serial gas and particle-SDC application
routing remain separate work.
