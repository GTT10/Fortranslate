# Spray-A KHRT normalization and child release — 2026-09-27

Implementation `fb97202be6ef9f7a4c8e0e2385ec9460a8a2c304` corrects the
Spray-A preparer's KHRT release normalization. Pinned
PelePhysics `6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`,
`Source/Spray/SprayParticles.cpp`, computes the average injection reference
as `num_ppp * average_diameter**3` for this one-jet case. The pinned input sets
`inject_ppp=0.1`; `SprayInjection.H` replaces that value only when it is negative.
The preparer had supplied diameter cubed alone. It now supplies
`0.1*(89e-6)**3 = 7.04969e-14 m³` and records it in preparation metadata.
The algorithm source is added to the case's verified input manifest (SHA256
`0e376f77df434ae7eec9fba39ab67af5ad4c20fb4c46529a9ff35b6d978b4d05`).
This correction follows the pinned rule; no threshold was tuned to make a test
pass. Fixed four-parcel event sampling remains an explicit adaptation.

## Superseded diagnostic

The preceding 12-event Debug run used the missing-factor input. Its live
creation observations are historical diagnostics, not qualification of the
correct input. After identifying the defect, only that run's verified launcher
PID 676110 was sent SIGINT. Solver PID 676113 terminated; driver exit 1,
elapsed 1066.8943 s. Artifacts remain in `build/spray-a-children/`, including
`interruption.json`. It is classified as deliberately interrupted, not passed
or naturally failed. The separate frozen Jet process was left running.

## Completed corrected runs

The existing GNU 13.3.0/Open MPI 4.1.6 Debug build and a separate clean Release
build both use the pinned 53-species dodecane Simple bundle. Release configure
uses the same options as [the original build](spray-a-connection-20260927.md),
with `CMAKE_BUILD_TYPE=Release` and build directory
`build/dodecane-simple-mpi-release`. Configure and selected spray target build
pass; logs: `build/dodecane-release-{configure,build}.log`. Release binary
SHA256: `9a2db9506db5529879370d49a32773cdfced2d559931a6172a030a5c6b677756`.
The solver source is the previously recorded tree implementation plus the
particle-ID output change; this correction changes preparation, not the solver.

Three corrected one-rank application runs completed:

| Directory under build/ | Mode | Events | End time | Result |
| --- | --- | ---: | ---: | --- |
| spray-a-weighted-reference | Debug | 2 | 0.6 µs | exit 0, 8 parcels, fuel error 1.192e-16 |
| spray-a-release-short | Release | 2 | 0.6 µs | exit 0, same gas and particle CSV values exactly |
| spray-a-release-children | Release | 8 | 2.8 µs | exit 0, 1264 parcels, fuel error 3.130e-16 |

The short runs took 17.1665 and 5.4031 s respectively. These are individual
run timings, not a performance benchmark. All gas/particle numeric CSV fields
agree exactly in `build/spray-a-release-short/debug-comparison.json`.

The extended Release run took 103.1687 s, accepted 16 steps and injected
32 parcels / `6.607651862518499e-10` kg of fuel. It reports ID watermark 1264,
zero computational splits and 1264 unique surviving parcels. Thus 1232 new
KHRT particle IDs were created in this no-splash case. Fuel CSV and signed
boundary/wall accounting agree within `3.129631485449313e-16` relative against
the unchanged `1e-8` gate; injected mass matches exactly. Energy ledger residual
is `8.524973262036273e-14` J; application composite inventory error is
`5.291951755036430e-16`. The minimum accepted interval is `1.0588e-22` s at an
event boundary; this does not qualify the time-rounding policy.

Commands:

```sh
python3 tools/prepare_spray_a_reference.py --bundle build/dodecane-simple-reference/dodecane_lu.json --reference-source build/pelec-pinned-hit-source --output build/spray-a-release-children --events 8
python3 build/run_spray_a_release_children.py
python3 tools/check_jet_inventory.py --run-dir build/spray-a-release-children --fuel-species NC12H26 --allow-khrt-children --output build/spray-a-release-children/inventory.json
```

The driver records exact command, binary hash and exit status in each
`execution.json`, unsets DISPLAY/library-path variables and uses the previously
recorded MPI environment. Short-run drivers are `build/run_spray_a_weighted.py`
and `build/run_spray_a_release_short.py`. Their inventory reports pass as well.
No other new application runs are counted. Fixed and selected mechanism
provenance, Python syntax, project contract and diff checks pass.

This verifies child release and conservative coupling on a coarse short
adaptation. Full 3 ms, upstream mesh/AMR/stochastic sampling, reference
trajectory/penetration comparison and physical validation remain open.
