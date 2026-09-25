# MPI EB spray CTU and signed reconstructed species — 2026-09-25

Numerical correction `188fc05`; spray integration `3dc714c`, based on
`cb41723`. PeleC pin remains `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

## Change and reference

The MPI spray wrapper forwards an optional extended geometry cache into the
coupled gas operator. The shared tree application now accepts `plm_ctu` for
MPI EB spray, including both film/no-film dispatch branches. Particle updates
remain transactional around the gas update. The gas checkpoint identity binds
the method; restarting CTU as PCM rejects without writing result output.

Initial pure-N2 spray unit cases rejected CTU double-transverse predictions:
new passive weights reached `-6.6115534049436878e-24` and
`-2.7929069933295029e-32`. This was not fixed by changing tolerances or clipping
those weights. Pinned `Source/Godunov.H::pc_transdd_passive` keeps signed
reconstructed weights; `Source/Riemann.H::riemann` passes them into EOS before
its existing star-species nonnegative clamp. Pinned PelePhysics
`6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`, `Source/Eos/Fuego.H::RPY2Cs/RYP2E`,
evaluates weighted gas constant, Cv and internal energy without a species-sign
check. The Fortran physical-state converter had incorrectly been reused for
these transient states.

The private reconstructed-state EOS now evaluates signed weights for primitive
Riemann calls, with finite/positive density, pressure, gas constant, Cv and sound
speed and unchanged NASA temperature bounds. Nonnegative inputs retain their
previous path. CTU transverse prediction retains signed weights. Ordinary
physical-state conversion still rejects negative composition, and cell recovery,
conservation checks and rollback tolerances are unchanged.

## Completed evidence

GNU Fortran 13.3.0, OpenMPI, Python 3.12.3, existing Ninja Debug build
`build/spray-amr-reactive-debug`; reference adapter uses `c++ -std=c++17 -O0
-ffp-contract=off`. No compiler/MPI configuration change.

| Test | Result | Seconds |
| --- | --- | --- |
| `regression_reactive_eb_hydro_3d` | PASS | 13.58 |
| `unit_reactive_pelec_riemann` | PASS | 0.00 |
| `spray_mpi_eb_reactive_2` | PASS | 1.06 |
| `spray_mpi_eb_reactive_ctu` | PASS | 2.22 |
| `spray_mpi_eb_reactive_ctu-subcycle` | PASS | 4.28 |

5/5 passed, 21.18 s total; the added physical-negative-composition assertion was
then independently rerun in the Riemann unit test (1/1, 0.05 s total).
CTU spray cases use four ranks, two owner patches and two idle ranks. They check
nonzero chemical exchange, parcel heating/evaporation, reflection/drag, composite
gas/parcel/wall conservation at the unchanged `3e-12` criterion and complete
rollback after gas rejection. Both common and recursive gas stepping passed.

The bounded pinned-C++ Riemann oracle now has 64 cases/896 values, including
negative weights `-1e-24` and `-.001` in four velocity regimes. Maximum scaled
error `3.111830054382345e-14`, unchanged criterion `3e-12`. This uses an independent
NASA7 SI adapter with exported coefficients, not a full PeleC executable.

The fixed-mechanism application checker passed with `--ctu --restart` after the
numerical correction: 2/4-rank gas and particle parity, wall reflection, exact
2-to-4-rank continued gas/particle outputs, budget agreement, overwrite protection,
method/context mismatch and truncated-payload rejection. Final time is 30 ns;
checkpoint at step 2 leaves a subsequent step to execute. This application run
used common gas intervals. Selected application, film-specific CTU behavior,
nozzle/diesel/LES physics and whole-solver parity were not qualified here.

Initial unit suite was 1 pass/2 failures. Diagnostic prints were removed and
original sources restored before applying the numerical correction. Logs of
those failures remain in `build/mpi-eb-spray-ctu-checks/LastTest-initial.log` and
`debug*-test.log`; final evidence is `LastTest-signed.log`,
`LastTest-physical-rejection.log`, `signed-oracle/report.json`, build logs and
`application-final/eb_spray_4c8vv1iu/`. The earlier application run
`application/eb_spray_bpcjzh8_/` passed before the signed-state correction and
is not the final verification. Compiler warnings remain.

Commands (repository root; MPI checks used DISPLAY/LD_LIBRARY_PATH/LIBRARY_PATH
unset, HWLOC_COMPONENTS=-x86, OMPI_MCA_hwloc_base_binding_policy=none):

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_reactive test_reactive_pelec_riemann test_reactive_eb_hydro_3d pelef_mpi_eb_spray_3d --parallel 2
ctest --test-dir build/spray-amr-reactive-debug -R '^(spray_mpi_eb_reactive_(2|ctu|ctu-subcycle)|unit_reactive_pelec_riemann|regression_reactive_eb_hydro_3d)$' --output-on-failure
python3 tools/check_pinned_riemann_kernel.py --reference /mnt/sn850/Projects/gtt2/external/PeleC --executable build/spray-amr-reactive-debug/tests/test_reactive_pelec_riemann --work-dir build/mpi-eb-spray-ctu-checks/signed-oracle
python3 tests/application/check_mpi_eb_spray.py --exe build/spray-amr-reactive-debug/pelef_mpi_eb_spray_3d --mpiexec /usr/bin/mpiexec --source . --work build/mpi-eb-spray-ctu-checks/application-final --ctu --restart
python3 -m py_compile tools/check_pinned_riemann_kernel.py tests/application/check_mpi_eb_spray.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Provenance, Python syntax, contract and whitespace passed. Mechanism SHA remains
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
The new CTest `spray_mpi_eb_ctu_restart` registration was inspected with `ctest -N`;
the equivalent application command above was run directly, not again through CTest.
