# Isolated parallel integration, 2026-11-06

Application update, 2026-11-06 14:29 JST: after user approval, the packaged changes
were applied to `/home/gtt13/projects/pelec` on main. All 30 patch outputs matched
their recorded hashes, and all 801 files in the final source manifest matched the
tested candidate. Existing staged changes were preserved; no commit or push was
performed. The following sections retain the pre-application qualification record.
The application receipt and preimages are in
`/mnt/sn850/Projects/gtt/docs/analysis/pelef_parallel_20261106/application_20261106_142935/`.

This candidate combines worker04 Simple transport reuse, worker07 ETAB CGS
conversion, worker08 static LES face-cp correction, and worker09 MPI optional
collective ordering/transactional AMR checkpoint publication. It additionally
repairs wall-born parcel continuation/checkpoint metadata and connects serial
species-bounded conservative reflux, including coupled held energy and effective
faces, using the same funding kernels and traversal as MPI.

The base is the supplied dirty-main snapshot, not a clean Git commit. Original
main and existing jobs were not changed. Mechanisms, reference limits and existing
tolerances were preserved. Numerical changes, corrected deposition fixture input,
and CMake registration proposals remain separate reviewable patches in the package.

## Executed integration checks

- Initial combined revision: 30 CTests executed, 27 passed, 3 known baseline
  failures, no skips. All five newly registered transport/LES/publication tests
  actually executed and passed. The failures remain unit_mixture_transport
  (NASA midpoint expectation), spray_les_dynamic_boundaries (temperature equality
  inconsistent with allowed state rounding), and spray_les_amr_level (root-level
  timestep controls the fixture). Their tolerances/expectations are unchanged.
- Final revision after spray/reflux followups: 23 affected CTests passed, no skips
  or failures. This includes active common/recursive gas reflux recovery/refusal
  at MPI2/3, real single/mixed wall-child continuation, MPI migration/restart,
  wall/transport MOL restart, AMR LES continuation and publication refusal.
- Four controlled I/O failures passed again against final executables, preserving
  the previous checkpoint bytes. Initial-revision records are preserved separately.
- MPI2 retained-internal-energy common/recursive lifecycle passed independently.
- Mechanism provenance, project contract, Python syntax, whitespace, source hashes
  and exact patch application are recorded in the delivery package.
- Independent numerical/state review, spray-format review, reflux review and
  validation-evidence review found no blocking issue within their recorded scope.

Worker05 remains separate. A paired same-library single-object Release comparison
of a 4-cubed reacting hotspot with transport reached 10 ns in all ten runs with
byte-identical CSVs. Five alternating pairs gave median CPU 6.74 versus 6.61 s,
about 1.93% shorter, not the 18--40% inverse-kernel figure. The known standalone
property regression and broader workload uncertainty remain; no worker05 source
is included in this candidate. This is a short fixed-mechanism application
measurement, not production or AMR throughput qualification.

## Bounded followup conclusions

Worker11 preserves valid wall-born N0 without enabling aerodynamic shedding.
Pending shed volume remains rejected. New binaries read old zero-count files,
whose bytes are unchanged, and use an explicit N0 tag for nonzero retained counts.
Older binaries intentionally reject that new tag. Negative/malformed payload
checks retain caller state, film, parcel, clocks and ledgers.

Worker12 closes the demonstrated serial/MPI low-level acceptance gap. Six topology
modes exercise full conserved-vector and held-energy transfers, effective faces,
budget closure and atomic refusal. New full-split active repair is gas-only;
stored-energy lifecycle and existing rank-changing restart are separate evidence.
An active repaired-state restart and a full held-energy trajectory with active
repair were not newly qualified. No causal claim about prior reference errors is
made, and no long trajectory/grid sweep was rerun.

Worker13 captured actual callback57 and all ten species/source inputs. At the
stall Y_OH approaches -5e-13 and the RHS points outward while dt repeats near
3e-22 s. An explicitly capped diagnostic rejects transactionally; this is not a
passing replacement for the original kinetic test. No generic progress guard,
mechanism edit, relaxed bound or clipping was adopted. That kinetic timeout,
EB short-step field mismatch/convergence, canonical C2, physical validation and
full release/install qualification remain open.

## Delivery

The durable package is
`/mnt/sn850/Projects/gtt/docs/analysis/pelef_parallel_20261106/integration_result/`.
It contains the combined patch, ordered original patches, before/after hashes,
command/test logs, independent reviews, the full actual-cell reproducer and a
source archive excluding historical build outputs. The integration source is
`/tmp/pelef-parallel-20261106/integration/source` with a synthetic baseline commit.
This record does not claim application to original main, publication, or overall
PeleF completion. Reproduction and application checks are in the package README.
