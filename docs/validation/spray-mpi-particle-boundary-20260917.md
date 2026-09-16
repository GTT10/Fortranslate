# Nonperiodic MPI particle coupling — 2026-09-17

Source: parent `748f908` plus this change. Source hashes and logs:
`build/spray-mpi-particle-boundary-checks/`. Existing GNU Fortran 13.3.0 Debug /
Open MPI 4.1.6 configuration.

`advance_mpi_mixture_coupling` now accepts optional physical boundary labels and
a boundary-change output. Sampling and feedback use the corresponding physical
stencils. After transfer/breakup, the serial parcel-boundary routine wraps periodic
coordinates, removes open-face escapees and reflects stationary-wall crossings.
Mixture density, heat capacity and internal energy define the effective liquid
used by that routine. Component inventories are cleared on removal; mass, momentum,
energy and mapped gas-species inventories enter a replicated global ledger,
positive into the gas-plus-particle system. Surviving parcels migrate using the
actual periodic axes. Input arrays, optional ID watermark and boundary output
commit only when the whole interval succeeds; failure returns zero boundary change.

Boundary handling also visits newly generated KH children. This change does not
connect mixed deposition/splash, moving-wall impact or film evolution to the
coupled driver. Stationary specular reflection is the behavior of the existing
serial parcel-boundary routine, not a substitute claim of wet-wall modeling.
Nonperiodic gas/particle split orchestration and application inputs remain open.

Completed commands:

```sh
cmake --build build/spray-mpi-migration-gnu-debug --target test_spray_mpi --parallel 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^spray_mpi_migration_[24]$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

Final focused tests passed 2/2 in 15.93 s, with no failures or skips. Provenance,
contract and whitespace checks passed.

New ballistic cases put three two-component parcels on rank zero: one exits the
lower x face, one reflects from the lower z wall, and one crosses into rank one's
slab. Heat, drag and evaporation are disabled to isolate boundary motion. The
survivor IDs, owners and positions match the expected trajectories. The global
ledger accounts for escaped mass/components, reflection impulse and both mixture
internal and kinetic energy. A longer interval with a one-substep cap rejects
after partial candidate advancement and exactly preserves gas, temperature,
particle positions/velocities and component masses, returning zero ledger.
Existing periodic coupling, breakup, film, restart and gas-boundary checks remain
in the same focused entries. An initial compile caught a long migration call;
wrapping that line fixed the build. No Python source changed.

This is focused numerical evidence, not physical spray qualification or a full
configuration matrix.
