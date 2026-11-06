# Recursive native MOL common-clock prefix retries, 2026-10-15

## Implemented scope

The serial/MPI recursive drivers now retry a refusing ordinary native MOL
common clock using its observed source prefix. Retries occur after the existing
entry funding payment/reset and before any descendant call or publication of
the accepted receipt. Each common trial retains the same physical entry,
ancestor endpoints and current funding stock; the common-step transaction
commits gas/T/H/I_R, faces, source receipt and ledgers only on acceptance.

The retry family descends through spatial weights 1/2, 1/4, ..., 1/4096, then 0,
with at most 13 retries after the original unit-weight attempt. Each failed
native trial can extend the observed source cache. A later trial consumes every
already observed packet and asks for an unobserved stage only when first
reached. Declined source packets and native-admitted caller reactor refusals
stop immediately. A successful prefix trial marks its accepted interval as
replayed. No conservation, reaction or EOS threshold was relaxed.

The scope requires chemistry, MOL, ordinary redistribution, independent H and
an energy-correction ledger, with built-in spatial physics. Routed/prepared
paths, caller spatial/predictor providers and MPI level-post callbacks retain
their existing behavior. Explicit complete-window replay disables these local
retries. The existing complete-window parent closure path is unchanged.

The optional `native_common_retries` int64 output is failure diagnostics,
initialized to zero and incremented for actual scheduled retry attempts.
It is distinct from accepted `level_steps` and physical ledgers. MPI agrees
its presence as part of the existing wrapper option profile.

## Fresh checks

Eight CTests passed, with zero failures and zero skips:

- `unit_mpi_eb_tree_held_natural_replay_3d_3`: 74.44 seconds.
- `unit_mpi_eb_native_mol_recursive_prefix_retry_3d_3`: 110.69 seconds.
- `unit_mpi_eb_native_mol_source_prefix_3d_3`: 10.70 seconds.
- `unit_mpi_eb_native_mol_prefix_consumer_3d_3`: 47.49 seconds.
- `unit_mpi_eb_native_mol_future_source_3d_3`: 232.34 seconds.
- `unit_mpi_eb_tree_held_scalar_replay_3d_3`: 76.72 seconds.
- `unit_mpi_eb_tree_held_native_funding_3d_3`: 637.30 seconds.
- `unit_mpi_eb_tree_held_native_material_3d_3`: 737.98 seconds.

The new test uses full native H2O2 chemistry, inert nitrogen at 1400 K,
two-level cut geometry and both `[0,2]` and `[0,0]` owner layouts on three ranks.
It verifies serial and MPI common-clock retry wiring and failure rollback:

- An accepted, deliberately unfundable old sink makes 13 native retries while
  calling its source only once; it never requests a future source.
- A newly accepted sink at the first new stage likewise makes 13 retries while
  retaining both observed packets, with only two source calls in total.
- Old/new source refusals (including an idle-rank refusal), a native-admitted
  caller reactor refusal and rank-mismatched diagnostic presence make no retry.
- Every failed outer transaction retains caller gas/T/H/I_R, gas/scalar faces,
  physical ledgers, accepted step counts, history and funding context.
- An ordinary successful pass makes no retry, advances root/two fine clocks and
  preserves exact serial/MPI gas/T/H/I_R, faces and ledger agreement.

These refusal fixtures deliberately have no admissible material balance. They
qualify retry/cache/rollback behavior; they do not demonstrate a successful
newly recovered partial trajectory. The separate bounded search in `probe.*`
looks for such a candidate using native physical gas/H inputs. A search hit is
not promoted evidence: it requires a reproducible serial/MPI fixture, accepted
source/flux/funding inventory checks and rollback qualification.

The consumer/capture, current-stock and four existing automatic complete-window
checks are regressions on the same linked source. Whole parent/descendant
partial-window replay, accepted recovery with material funding, application or
installed-binary audits, pinned solver parity and physical validation remain
separate open gates.

## Identity and commands

The dirty main HEAD remains `0a0e00ebe9138145fbb0e0e745665d0d768df009`.
GNU Fortran 13.3.0 / OpenMPI 4.1.6 Debug, fixed/selected h2o2_full, MPI/spray and
SUNDIALS OFF remain unchanged. Four source/test/CMake files change from the
preceding consumer delivery. Artifacts, exact source SHA-256, incremental diffs,
linked identities for two libraries/two executables, commands/results, native
logs and JUnit are in `build/native-mol-recursive-prefix-retry-20261015`.
The initial development build's overlong argument line failed compilation;
that log is retained, the line was wrapped, and the failed build is not counted
as passed. Provenance, 14 Python provenance tests, Python syntax, project
contract and diff whitespace are also required before promotion.

The first test candidate passed four prefix/funding checks but failed four MPI
parent-window regressions. Prefix capture had expanded the accepted interval
reaction-source field to inactive owned levels, violating its sparse selected
level schema. Whole-entry prefix I_R is now stored separately; accepted interval
records retain their original schema. The new retry unit also asserts that
schema directly on accepted recursive records. All failed-candidate source,
linked files and JUnit are retained in `failed-v1`; none are counted as passed.
