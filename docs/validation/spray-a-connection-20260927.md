# Spray-A connection — 2026-09-27

Implementation `c2df43e8c588ed0c4ebad3f76e08bdf10350433f`, based on `eb51f3d7`.
The selected solver was built from that base; the implementation adds input
preparation and a fuel label to the inventory checker, without solver changes.
Seven case files in `references/spray_a_case.json` were compared byte-for-byte
with PeleC `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8` using the local Git
object store. Preparation verifies their SHA256, mechanism and imported Simple
fits before writing a new directory. It retains copied upstream inputs.

## Clean build and execution

GNU Fortran 13.3.0/Open MPI 4.1.6, Debug, MPI/spray ON, tests/SUNDIALS OFF:

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake -S . -B build/dodecane-simple-mpi-debug -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran -DMPI_Fortran_COMPILER=/usr/bin/mpifort -DMPIEXEC_EXECUTABLE=/usr/bin/mpiexec -DPELEF_ENABLE_MPI=ON -DPELEF_ENABLE_SPRAY=ON -DPELEF_ENABLE_TESTS=OFF -DPELEF_ENABLE_SUNDIALS=OFF -DPELEF_MECHANISM_BUNDLE=/home/gtt13/projects/pelec/build/dodecane-simple-reference/dodecane_lu.json -DPELEF_MECHANISM_SOURCE=/home/gtt13/projects/pelec/build/dodecane-simple-reference/dodecane_lu.yaml
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/dodecane-simple-mpi-debug --target pelef_mpi_reactive_eb_tree_3d_selected_spray --parallel 2
python3 tools/prepare_spray_a_reference.py --bundle build/dodecane-simple-reference/dodecane_lu.json --reference-source build/pelec-pinned-hit-source --output build/spray-a-connected
python3 build/run_spray_a.py
python3 tools/check_jet_inventory.py --run-dir build/spray-a-connected --fuel-species NC12H26 --output build/spray-a-connected/inventory.json
```

The driver unsets DISPLAY and both library-path variables, sets
`HWLOC_COMPONENTS=-x86`, `OMPI_MCA_hwloc_base_binding_policy=none`, and runs
`/usr/bin/mpiexec -n 1 <selected-spray-binary> tree.nml spray.nml` in the case
directory. Configure/build logs: `build/dodecane-simple-configure.log`,
`build/spray-a-build.log`. Both pass, including MPI ABI audit. Existing
array-temporary and executable-stack warnings remain; no installed-binary
or full configuration-matrix qualification is claimed.

Binary SHA256:
`1756bfabb90029773a58e5de857e9441dda3ac61b5e8f95765c0db92f8a654c4`.
Execution exits 0 in 17.0585 s; `build/spray-a-connected/execution.json`
records the command and identity. It completes 3 steps to 0.6 µs, two events,
eight retained parcels and 128 gas cells. Expected injected fuel is
`2.7114158689350002e-11` kg. Gas + liquid − signed boundary/wall fuel agrees
within `1.1916926152145067e-16` relative, below the unchanged `1e-8` gate.
Injected mass matches exactly. Energy ledger residual is
`-2.6941190706293583e-14` J. Input/output hashes and inventories are in
`inventory.json`; application composite inventory error is `1.6724e-16`.

## KHRT intervention and limits

One additional run removes only breakup-related options:
`build/run_spray_a_no_breakup.py`, `build/spray-a-no-breakup/`. It also completes
and passes fuel accounting. With KHRT, parcel multiplicities grow by factors
1.0030951–1.0046102 relative to initial injection counts; without breakup,
counts remain unchanged to `6.66e-16` relative. This establishes active
representative-drop breakup evolution in this fixture, not child release or
upstream transient parity. Eight parcels remain in both runs. The recorded
comparison is `build/spray-a-connected/breakup-intervention.json`.

The [case guide](../../cases/spray_a_reference_3d/README.md) enumerates physical
inputs and adaptations. Initial preparer drafts under `build/spray-a-connection`
and `build/spray-a-connection-final` were not run; only `spray-a-connected`
uses the final pressure-dependent boiling configuration. Full 3 ms,
32×32×160 plus AMR, stochastic injection, child release, reference penetration
comparison and physical validation remain open.

The existing decane short output passes the updated inventory checker
(`build/spray-a-decane-checker-regression.json`); no decane calculation was
rerun. Fixed h2o2 and selected dodecane provenance, Python syntax, project
contract and `git diff --check` pass. Two application runs were completed.
