# Owner-local MPI SDC checkpoint I/O — 2026-09-26

Implementation `618bf656`, base `2c355f0a`; PeleF 0.246.0.
Pinned PeleC `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

## Behavior

The MPI checkpoint library now supports the serial SDC V3 plain and V4 indexed
formats. The replicated external-source ledger resides in metadata; reaction
estimates in `metadata%reaction_source` are owner-local. Presence, layout and
finite values are checked collectively before writing. Source records follow
the gas records and preserve their signed values. Existing V1/V2 paths remain.

MPI I/O writes each rank's own source records at calculated byte offsets. On
read, rank zero scans source record positions one cell at a time, broadcasts
only offsets, and each new owner seeks to its assigned records. It does not
retain a full source/state tree. Candidate validation includes reaction-field
finiteness and layout before publishing snapshot or ownership. Gas indexed
record offsets retain their existing meaning.

## Completed focused checks

Artifacts: `build/sdc-mpi-checkpoint-checks/` (local).
GNU Fortran 13.3.0, Ninja Debug, existing tests-enabled
`build/spray-amr-reactive-debug`; Open MPI 4.1.6 through absolute
`/usr/bin/mpiexec`. No compiler/backend changes.

- Provenance passed: `python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json`;
  SHA `0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
- Builds passed: `env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build
  build/spray-amr-reactive-debug --target test_mpi_amr_patch_tree_checkpoint_3d
  --parallel 2`. Final test build log `build-test-final.log`; warnings remain.
- Direct four-rank unit executable passed, exit 0 (`run-final.log`). It includes
  existing checkpoint preflight/late-EOS rollback checks plus V3/V4 signed-source
  round trips, source-only records without reaction estimates, serial read/write
  interchange and **4-rank write to 2-rank read** using a split communicator.
  Reaction estimates and external inventories are exact; non-owned source arrays
  remain absent. A NaN in the final reaction record rejects the candidate and
  preserves prior source fields and owner mapping.
- Command, from the artifact directory:
  `env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86
  OMPI_MCA_hwloc_base_binding_policy=none /usr/bin/mpiexec -n 4
  /home/gtt13/projects/pelec/build/spray-amr-reactive-debug/tests/test_mpi_amr_patch_tree_checkpoint_3d`.
- An initial CTest filter for `_4` found **zero tests** because registered variants
  are `_2` and `_3`; this is not counted as a pass. The explicit four-rank run
  above completed. An earlier direct run passed before adding the NaN check;
  it is not counted as a separate evidence class. Logs retain display-environment
  authorization warnings from the launcher, followed by all test PASS markers.
- Python syntax for `tools/check_project_contract.py`, project contract,
  Fortran line lengths and `git diff --check` passed.

This is checkpoint interchange, not a continued MPI SDC application trajectory.
Application routing, regridding of MPI reaction estimates and recursive MPI SDC
remain open. No selected build, full matrix, install audit, full-solver pinned
comparison or physical validation was performed.
