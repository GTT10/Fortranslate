# Same-level EB peer flux reconciliation — 2026-09-17

Implementation: `df9a695fc6ebfb71e9c77f89b15c449ea0655a43`, based on
`5889f25` on local `main`.

`reconcile_reactive_eb_peer_flux_3d` takes an authoritative owner flux and an
already advanced receiver patch. For each shared open face it replaces the
receiver's integrated flux with the owner's physical transfer, correcting the
receiver conserved state by the signed difference divided by fluid volume.
Open-area and Cartesian-area factors are applied once; the stored per-open-area
receiver flux accounts for any roundoff-sized area difference. Owner data is
read-only. Receiver flux, conserved state and EOS temperature commit together.

The pair operation checks geometry validity, same spacing, aligned origins,
nonoverlapping boxes, matching apertures and valid flux/state shapes. It does
not choose an owner, validate the whole collection or replace the peer halo
API's tangential-centroid checks. Callers must establish a consistent shared
geometry and deterministic ownership policy. EOS rejection is transactional;
reconciliation does not guarantee that the authoritative limited flux will
preserve a receiver's positivity. A hierarchy can reject and retry the interval.

This is a building block, not integrated multiple-patch or MPI stepping.

## Independently executed checks

GNU Fortran 13.3.0, Debug/Ninja, existing
`build/selected-tree-mpi-debug`; the unit is serial.

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug \
  --target test_eb_reactive_reflux_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug \
  -R '^unit_eb_reactive_reflux_3d$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

**1/1 passed**, zero failures/skips, 0.01 s (0.04 s total).
Provenance, contract and whitespace checks passed; no Python changed.
Mechanism source SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

The expanded test retains coarse/fine reflux cases and adds all six peer
orientations with partial shared faces and cut geometry. Receiver integrated
inventory changes equal the independent face-transfer difference for every
conserved component. Corrected temperatures remain positive. Reapplying an
already synchronized flux leaves conserved state exactly unchanged. An
unphysical correction preserves receiver state, temperature and all flux arrays;
overlapping boxes reject as well. This test uses prescribed integrated fluxes,
not a multiple-patch time integrator or a viscous-stencil parity experiment.

Logs: `build/eb-peer-flux-checks/{build,ctest}.log`.
No broader configuration matrix, installed audit, external solver parity,
convergence study or physical validation was rerun.
