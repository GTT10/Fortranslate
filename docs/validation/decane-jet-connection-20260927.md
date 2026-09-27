# Decane jet application connection — 2026-09-27

Base: `219d7f8a`; implementation is the accompanying
`prepare_jet_reference.py`, pinned case manifest and case guide. Numerical
solver source is unchanged. This is a short adaptation of pinned PeleC
`bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8` Spray-Jet, using the
[pinned decane mechanism](decane-spray-import-20260927.md).

The preparation tool verifies the exact upstream `example.inp` and `roi.dat`
SHA-256 values and the mechanism/transport provenance. It preserves those
original files, converts rate units and writes gas, particle and level-set
inputs. Preparation refuses an existing output directory. A deliberately
altered ROI file is rejected before output creation
(`build/jet-preparation-rejection/rejection.log`).

## Completed execution

Exact configuration and execution commands are in the
[case guide](../../cases/spray_jet_3d/README.md). GNU Fortran 13.3.0 Debug,
OpenMPI 4.1.6, MPI/spray ON, tests/SUNDIALS OFF; fresh build
`build/decane-simple-openmpi-debug`. The selected MPI spray target builds and
one rank completes three `2e-7 s` SDC steps to `6e-7 s`, processing two
injection events. Logs:

- `build/decane-simple-openmpi-configure.log`
- `build/decane-simple-openmpi-build.log`
- `build/decane-jet-connection/fortran.log`

The first configuration attempt used automatic MPI detection, which selected
Intel MPI incompatible with GNU `mpi_f08`. It failed at configuration and is
retained in `build/decane-simple-mpi-configure.log`. The successful fresh
directory explicitly selects `/usr/bin/mpifort` and `/usr/bin/mpiexec`.
Existing compiler real-comparison and runtime array-temporary warnings remain.

## Independent inventory reconstruction

The exact first two active ROI intervals supply
`1.0425546181770002e-11` and `1.7375445019439998e-11 kg`. Dividing each interval's
mean rate by liquid density, nozzle area and discharge coefficient gives
`14.385466164356453` and `23.975134929127897 m/s`.

The output contains 128 gas cells and eight unique, positive-mass parcels.
Integrating fuel density using each CSV fluid volume and summing parcel
mass times multiplicity gives:

| Quantity | kg |
| --- | ---: |
| Expected injected fuel | `2.7800991201209998e-11` |
| Final gas fuel | `2.0307791120968503e-15` |
| Final liquid fuel | `2.7798960422097902e-11` |

Adding the reported fuel boundary/wall ledger reproduces the expected total
to the printed double precision (relative error zero; focused gate `1e-8`).
This checks the CSV inventories independently of the driver's final inventory
accumulator. The driver's energy residual is `-1.76368e-15 J`; its reported
composite inventory error is `1.80923e-16`. This is numerical conservation
evidence, not experimental spray validation.

Reproduction and results:

```sh
python3 build/decane-jet-connection/check_inventory.py
```

`preparation.json` records interval masses, velocities, input identities and
adaptations. `inventory-check.json` records the focused check. Provenance,
Python syntax, project contract and diff checks pass. No CTest matrix,
rank-changing restart or full C++ jet trajectory comparison was run.

## Scope

Nitrogen gas remains at the reference 440 K / 6 MPa initially, with 363 K
decane, a 90 micrometre nozzle and 4.5 micrometre drops. The adaptation uses
4x8x4 all-fluid cells without AMR, deterministic four-parcel events, averaged
event velocities, a nozzle center displaced 1e-8 m inside the wall, and 300 K
fuel reference temperature instead of the source's out-of-range 298.15 K.
The source thermo bounds and generated coefficients remain unchanged.
The 3.8 ms refined stochastic reference jet, penetration, breakup and physical
validation remain open; this result establishes the selected decane
gas/particle application connection only.
