# MPI particle-impact and film creation coupling — 2026-09-17

Source: parent `78bad52` plus this change. Source hashes and logs:
`build/spray-mpi-impact-coupled-checks/`. Existing GNU Fortran 13.3.0 Debug /
Open MPI 4.1.6 configuration.

`advance_mpi_mixture_coupling` accepts film state and wall settings together,
requiring physical boundaries and an ID watermark. Each trial performs existing
frozen-gas transfer and breakup, then detects planar crossings, queries distributed
film footprints and resolves mixed impacts before gas feedback commits. KH children
inherit their original parent's trajectory start and sampled pressure for that
collision batch. Splash children consume globally reserved IDs. Film/particle
candidates and the ID watermark are discarded together when feedback or a later
substep fails. Dead parents are compacted between accepted substeps; final films
and free particles migrate to their spatial owners.

The impact API now excludes periodic-face crossings when selecting the first
physical hit and wraps tangential periodic coordinates at the impact location.
Thus a periodic crossing earlier in the same trajectory cannot suppress a later
wall event. Open faces still participate in earliest-physical-face selection.

Global impact wall increments are accumulated separately from local open/reflecting
boundary increments, then subtracted exactly once after local reduction. The
reported boundary ledger is positive into gas, particles and films. The outer
reaction/flow splitter passes its private film candidates through both particle
half-steps, so created films reach the existing terminal film half-step and all
state stays inside the outer transaction.

This synchronous MPI ordering applies breakup before the collision batch, uses
pre-batch film heights, and starts splash descendants at the accepted substep end.
It differs from serial sequential impact/breakup and wet-film accumulation; no
blanket trajectory parity is claimed. Application film controls remain open.

Completed commands:

```sh
cmake --build build/spray-mpi-migration-gnu-debug --target test_spray_mpi --parallel 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^spray_mpi_migration_[24]$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

Final focused tests passed 2/2 in 17.05 s with no failures or skips. Provenance,
contract and whitespace checks passed.

New 2/4-rank cases include an empty rank and multiple source owners. Low-speed
impacts deposit whole-parcel component inventories, consume parents without new
IDs and redistribute films to rank zero. Their paths cross a periodic x face
before striking the lower z wall. Faster dry impacts reserve four splash IDs per
parent and retain lineage through owner migration. Total gas/particle/film
inventories agree with the wall ledger within 2e-12 of initial inventories;
for splash momentum, the initial vector norm sets the scale for all components.
A one-substep cap rejects an interval after candidate deposition and exactly
preserves original gas, temperature, parent masses/positions/composition, empty
film state and ID watermark; the output ledger and substep count are zero.
Existing film-split rollback, periodic breakup and nonperiodic flow tests remain
in the same entries. Combined KH shedding and collision was not separately exercised.
No Python source changed; no full matrix or physical impact validation is claimed.

The new oblique-splash check initially rejected a 4.4484e-23 transverse momentum
residual because that initially zero component used a 1e-12 inventory floor
(2e-24 absolute threshold). The initial momentum norm was about 1.2e-5. The new
check now scales all momentum components by that vector norm, retaining the
2e-12 relative coefficient; solver arithmetic and pre-existing tests were not
changed. The failure and measured residual remain in
`ctest-tangential-budget-failure.log` and `diagnostic.log`.
