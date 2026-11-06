# Density-coupled momentum interpolation, 2026-10-22

## Cause and bounded implementation

The smooth 900 K H2/N2 contact exposed a velocity defect in the existing
`separate_momentum` halo policy. Its scalar limiter can flatten the density
slope (including in response to roundoff-sized stored-temperature extrema),
while independently limited momentum retains a nonzero slope. The ratio
momentum/density then varies even though every donor has the same velocity.
The incoming fine halos show maximum x-velocity error 0.0887952547 m/s in
stage 1 and 0.0888144921 m/s in stage 2; incoming halo pressure and recovered
temperature remain at roundoff. These logs describe the original halos,
before the diagnostic candidate repairs them.

`interpolate_eb_linear_halo_3d` now accepts a further explicit optional
`consistent_velocity` flag, which requires `separate_momentum`. With donor
velocity v0=m0/rho0, the candidate limits the neighboring residual
m-v0*rho, then forms each directional momentum slope as
v0 times the limited density slope plus the limited residual slope.
Constant velocity is therefore retained even when the scalar limiter clips
density. The interpolation remains linear in child offset, preserving the
average momentum of symmetric children. Density/species/scalar interpolation
and the default and existing separate-momentum paths are unchanged.

The reference probe opts in only with the seventh argument
`--consistent-velocity`, after `--smooth-contact`. It replaces fine ghost
momentum in the paired-pressure diagnostic, updates the ghost kinetic energy
consistently with held energy, and retains paired pressure reconstruction.
Valid leaf states are not overwritten by halo repair. This is a candidate
core interpolation option plus diagnostic integration; ordinary serial/MPI
applications do not select the new option. General velocity gradients,
nonuniform temperature, cut geometry, kinetic/thermodynamic positivity and
full held-MOL closure still need qualification before application exposure.

## Completed comparison

Both modes use the same smooth periodic contact, 900 K, 101325 Pa, x velocity
10 m/s, four 50 ns SSPRK2 steps and the existing accepted-step reset.
The baseline mode retains the original halo policy. The paired mode adds
the density-coupled momentum candidate to the earlier pressure candidate.

| Metric | Earlier pressure-only candidate | Pressure plus velocity candidate |
| --- | ---: | ---: |
| Pressure offset at 200 ns (Pa) | 0.001550295812 | 2.91038e-11 |
| Stage 2 fine pressure offset before sync (Pa) | 1.93274347e-4 | 1.45519e-11 |
| Stage 2 coarse pressure offset after sync (Pa) | 9.66369844e-5 | 2.91038e-11 |
| Cumulative absolute reset energy (J) | 4.77372e-11 | 5.88898e-17 |

The maximum over all accepted candidate clocks is 4.36558e-11 Pa, below the
1e-8 Pa equilibrium diagnostic bound. The worst independently checked
correction-accounted inventory residual is 2.23299e-14. E-H-KE is zero after
each accepted reset. The baseline accepted history is byte-identical to the
previous smooth900 baseline, confirming unchanged behavior when the candidate
is not selected. This is not a native-PeleC comparison or physical validation.

The earlier apparent second-stage thermodynamic problem is thus explained,
for this contact, by a velocity-inconsistent halo feeding the spatial update.
Reflux must remain conservative; removing synchronization was not needed.
The next integration step is to carry the candidate through actual held-MOL
closure and qualify nonuniform-temperature/cut controls before serial/MPI
application adoption.

## Verification and artifacts

Artifacts: `build/eb-contact-stage-20261022/`. The one-step `probe-before`
run records the original halo defect; `probe-after` and four.dat freeze the
completed candidate run. `record.json` records commands, return codes, binary
and input hashes. `checks.json`, `location.json`, `sync.json` and
`contact-gates.json` hold independent reductions; source-identity.json records
the dirty source hashes and both binaries. Build/rebuild and CTest logs are
retained. The completed tool event was consumed before report analysis.

Base HEAD: `0a0e00ebe9138145fbb0e0e745665d0d768df009`, dirty main.
GNU Fortran 13.3.0 Debug, existing isolated OpenMPI registration build.
Build warnings include existing conversion/comparison diagnostics and the
probe executable-stack linker warning; no compiler/MPI configuration changed.

Completed checks:

- `cmake --build build/eb-eos-application-20261021/isolated-openmpi/registration --target test_eb_tree_reacting_reference test_eb_reactive_transfer_3d --parallel 2`
- `ctest --test-dir build/eb-eos-application-20261021/isolated-openmpi/registration -R '^unit_eb_reactive_transfer_3d$' --output-on-failure`: 1 passed. New assertions reproduce the old constant-velocity failure, preserve disabled-policy bits, check constant positive/negative velocity, child-average conservation, resolved velocity slope and transactional refusal.
- `python3 build/eb-contact-stage-20261022/run.py`: 2 policy trajectories complete; all inventory, cell and stage reductions pass.
- `python3 -m unittest discover -s tests/python -p 'test_eb_pressure*.py'`: 19 passed, zero failures/skips.
- Mechanism provenance, Python syntax, `python3 tools/check_project_contract.py --source .` and `git diff --check`: passed.

No MPI, installed-binary, full configuration matrix or native-reference rerun
is included in these counts.
