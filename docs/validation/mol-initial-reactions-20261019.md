# Fresh MOL reaction-history initialization, 2026-10-19

Regular serial/MPI tree applications now compute the initial reaction estimate
before their first reacting MOL predictor. Previously they used zero until the
first accepted reaction update. The production initializer also passes the
pinned nonuniform TGReact comparison without importing reference reaction
history. This connects the fresh-start requirement identified in the
[preceding operator comparison](tgreact-amr-reference-20261019.md).

## Numerical change

`initialize_amr_tree_reaction_source_3d` and its MPI counterpart run the existing
SDC reaction packing/backend on private state copies with zero external rates.
They publish only the reaction estimate. Gas, temperature and actual UEINT
remain unchanged; failure on any level or owner preserves the complete prior
estimate. Idle MPI ranks do not call the reactor.

The root uses the selected coarse timestep. Recursive finer levels use the
successively divided intervals; common clocks use the same interval on every
level. Ordering follows pinned `PeleC::post_init`: root reaction initialization,
restriction of the still-zero fine estimates, then initialization of finer
levels. The root's covered estimates consequently remain zero. The existing
reaction backend supplies CVODE or the native implicit integrator; rotation
and actual held internal energy are forwarded through its existing interface.

The application invokes this only for a fresh reacting `mol` start, after
first-step/event/final-time limits. If the first attempt fails, it recomputes
the estimate at the reduced timestep. After an accepted step, the normal
persisted history takes over. Restarts bypass this initializer regardless of
the restored step number. Checkpoint layout/context were not changed; existing
checkpoint source histories are neither discarded nor retroactively reseeded.
Fresh reacting MOL results can change relative to the former zero-history
startup. Non-MOL and reaction-disabled routes do not activate this change.

## Completed checks

- One three-rank unit CTest passed in 0.58 s. It covers common/recursive initial
  intervals, nonzero synthetic H2-to-H conversion, the pinned root ordering,
  exact serial/MPI estimates, unchanged gas/T/UEINT, an idle rank, and rollback
  after a fine-level failure following successful root work. The synthetic
  conversion is an interface/transaction check, not physical chemistry evidence.
- The retained native TGReact plotfiles were compared with one new Fortran
  trajectory using `--initialize-history`. All prior trajectory acceptance
  limits passed unchanged. Temperature difference is at most `1.210e-5 K`;
  hydrogen-consumption relative difference is `2.244e-5` (0.002244%). Initial
  gas/T/UEINT remain exactly imported. Initial reaction estimates are now
  computed independently; their differences from native estimates are recorded
  and are not claimed to be bitwise identical.
- Five selected-LiDryer application runs passed: serial continuous, two-rank
  continuous, two-rank stop at step one, three-rank resume and serial resume.
  Final field and budget files are byte-identical. Each fresh run reports one
  initialization; both resumed runs report none. The checkpoint remains unchanged.

The application case uses a nonuniform 2000–2050 K hotspot, a 4³ root and one
2³ fine patch, two 100 ns coarse steps, recursive clocks and CVODE. Hydrogen
consumption is `2.329e-5` of the initial inventory. Independently integrated
CSV fields agree with native inventory rows; maximum scaled gas/energy budget
residual is `8.147e-12`. H/O/N relative drifts are respectively
`2.887e-15`, `7.772e-16`, `2.221e-16`; species positivity and sum closure pass.
Signed momentum residuals use absolute cell inventories for normalization.
The existing `source_change` column contains external forcing, so chemical
conversion is audited by element inventories rather than requiring each
species inventory to be constant.

Two earlier application checks stopped after a successful serial calculation
because the new checker incorrectly normalized zero net transverse momentum,
then treated chemical species conversion as an external-source budget failure.
Those checker mistakes were corrected without changing solver output or any
reference trajectory tolerance. A preceding PPM input was refused before
stepping and was corrected to the required characteristic PLM. All attempts
remain in the artifact directory.

## Identity and reproduction

Branch `main`, base `0a0e00ebe9138145fbb0e0e745665d0d768df009`, with substantial
pre-existing uncommitted work retained. GNU Fortran 13.3.0 Debug, Open MPI 4.1.6,
CMake 3.28.3 and SUNDIALS 7.2.0. No compiler/MPI/backend configuration changed.
Both fixed and selected frontends use the changed shared drivers; application
execution in this increment was selected LiDryer only.

Artifacts under `build/mol-initial-reactions-20261019/`:

- `source-identity.json`, `source/`, `before/`, `working-tree.patch`,
  `git-status.txt`, `CMakeCache.txt`: exact source, changes and toolchain.
- `reference-runs/matched_tgreact_oz6vuu1n/checks.json`: computed-history
  reference comparison; the native run was reused, not rerun.
- `app-runs/initial_mol_ffwdug9w/checks.json`: completed five-run lifecycle,
  exact commands, binary/input hashes, durations and inventory results.
- `unit.log`, `build.log`, `unit-build.log`, `app-build.log`: build/unit results.
- `app-runs/initial_mol_ry4ulcv7/`, `initial_mol_7lebzeri/`,
  `initial_mol_jaxx28qv/`: startup refusal and the two checker failures.

Build targets: `test_amr_tree_reacting_reference`,
`test_amr_tree_initial_reaction`, `pelef_amr_tree_reactive_3d_selected`, and
`pelef_mpi_amr_tree_reactive_3d_selected`, in
`build/reacting-hot-diagnosis-20261019/registration`, using `--parallel 2`.
Focused checks:

```sh
HWLOC_COMPONENTS=-gl,-x86 OMP_NUM_THREADS=1 ctest --test-dir build/reacting-hot-diagnosis-20261019/registration -R '^unit_amr_tree_initial_reaction_3$' --output-on-failure
```

The exact orchestration commands are recorded in `commands.json`. The reference
check uses `tools/check_tgreact_amr_reference.py --initialize-history`
with the previously captured native reference. The application checker is
`tools/check_amr_tree_initial_reaction.py --serial <selected serial> --mpi
<selected MPI> --mpiexec /usr/bin/mpiexec --work <artifact directory>`.
It uses repository case inputs `mol_initial_reactions.nml/.dat` under
`cases/amr_tree_reactive_3d`.

Mechanism provenance plus regeneration, Python syntax for both checkers,
project contract and `git diff --check` passed; their logs are retained. There
were no skips in the selected CTest or final comparisons. No full matrix,
installed-binary audit, fixed-mechanism application, rotating application,
three-level initialization, or earlier-checkpoint compatibility run was made.
Ordinary source initialization is connected; full-duration/convergence,
EB initialization and the nonreacting MOL2 policy difference remain separate
completion work. The result does not qualify physical ignition.


## 2026-10-28 current-source application rerun

The selected LiDryer serial/MPI applications and initial-reaction unit were
rebuilt with the existing GNU Debug/OpenMPI/SUNDIALS configuration. The
three-rank unit passes1/1. All five original application runs pass: serial,
MPI2, checkpoint stop on MPI2, resume on MPI3, and resume in serial. Fresh
runs initialize reaction history once; both restored runs do not reseed it.
CSV fields and budgets match exactly, and the input checkpoint is unchanged.
Maximum scaled budget residual is2.314365e-12; H/O/N relative inventory drift
is2.220447e-16. Source, executable and retained output/checkpoint hashes were
independently checked before extending the checker. Commands and results are
in `build/amr-reacting-lifecycle-current-20261028/verified.json` and
`record.json`. This confirms application equivalence, separately from the
fixed-PeleC operator comparisons and physical validation.

The checker now accepts `--subcycling 0` for common level timesteps and
`--transport` to enable the selected application's molecular transport.
The default remains recursive and transport-off. It keeps per-run CSV and
budget snapshots with hashes, allowing independent reinspection of continuous
and resumed results. All four combinations pass their five-run lifecycle (20 completed runs).
Retained per-run CSV and budget hashes independently reproduce exact serial,
MPI2, resumed MPI3 and resumed serial identity in each combination. Fresh
initialization markers, absence of restart reseeding, unchanged checkpoint
and source/executable identities were rechecked. Maximum scaled budget
residual across the matrix is4.085348e-12; maximum element drift is4.440893e-16.
The unchanged limits remain1e-9. Evidence is in
`build/amr-reacting-lifecycle-matrix-20261028/verified.json`. Application
transport uses the selected bundle's loader; this equivalence matrix does
not substitute for the separate pinned-Simple fixed-reference probe.
