# Static two-level reactive 3D AMR

Run the public ratio-two, full-H2O2 entropy-wave case from a build directory:

```sh
./pelef_amr_reactive_3d \
  ../cases/amr_reactive_3d/amr_entropy_wave.nml
```

The 8-by-8-by-8 periodic coarse grid refines coarse cells 3 through 6 in all
three directions. The fine patch takes two SSPRK2 substeps per coarse step,
uses linearly time-interpolated coarse boundary states, replaces all six
coarse/fine interface fluxes through reflux, and averages the covered cells
back to the coarse level. The application writes separate coarse and fine CSV
files and rejects a result that violates composite conservation, species
closure, physicality, or average-down synchronization.

The checkpoint/restart parity sequence uses the same physical case:

```sh
./pelef_amr_reactive_3d \
  ../cases/amr_reactive_3d/amr_entropy_wave_checkpoint.nml
./pelef_amr_reactive_3d \
  ../cases/amr_reactive_3d/amr_entropy_wave_restart.nml
python3 ../tools/compare_amr_reactive_3d_restart.py \
  --reference-coarse amr_reactive_3d_coarse.csv \
  --reference-fine amr_reactive_3d_fine.csv \
  --restart-coarse amr_reactive_3d_restart_coarse.csv \
  --restart-fine amr_reactive_3d_restart_fine.csv
```

The first command stops after coarse step four and writes a versioned
checkpoint. The second restores both levels and the original conservation and
reflux history, completes the requested final time, and must reproduce both
uninterrupted CSV files byte for byte. A checkpoint is accepted only when its
schema, full NASA7 species data, solver, boundary mode, mesh, patch, CFL,
physics switches, array extents, physical state, temperature recovery, and
coarse/fine synchronization match the requested run.

The rank-local sparse application accepts the same namelists and an output
prefix. PCM rank-count and changed-rank restart parity can be reproduced with:

```sh
mpiexec -n 1 ./pelef_mpi_amr_reactive_3d \
  ../cases/amr_reactive_3d/amr_entropy_wave.nml mpi_amr_np1
mpiexec -n 8 ./pelef_mpi_amr_reactive_3d \
  ../cases/amr_reactive_3d/amr_entropy_wave.nml mpi_amr_np8
mpiexec -n 2 ./pelef_mpi_amr_reactive_3d \
  ../cases/amr_reactive_3d/amr_entropy_wave_checkpoint.nml checkpoint_np2
mpiexec -n 4 ./pelef_mpi_amr_reactive_3d \
  ../cases/amr_reactive_3d/amr_entropy_wave_restart.nml restart_np4
python3 ../tools/compare_mpi_amr_reactive_3d.py \
  --reference-prefix mpi_amr_np1 \
  --candidate-prefix mpi_amr_np8 \
  --candidate-prefix restart_np4 \
  --checkpoint amr_reactive_3d.chk
```

The characteristic-PLM path uses the corresponding `*_plm*.nml` inputs:

```sh
./pelef_amr_reactive_3d \
  ../cases/amr_reactive_3d/amr_entropy_wave_plm.nml
mpiexec -n 8 ./pelef_mpi_amr_reactive_3d \
  ../cases/amr_reactive_3d/amr_entropy_wave_plm.nml mpi_amr_plm_np8
mpiexec -n 2 ./pelef_mpi_amr_reactive_3d \
  ../cases/amr_reactive_3d/amr_entropy_wave_plm_checkpoint.nml \
  checkpoint_plm_np2
mpiexec -n 8 ./pelef_mpi_amr_reactive_3d \
  ../cases/amr_reactive_3d/amr_entropy_wave_plm_restart.nml restart_plm_np8
python3 ../tools/compare_mpi_amr_reactive_3d.py \
  --reference-prefix amr_reactive_plm_3d \
  --candidate-prefix mpi_amr_plm_np8 \
  --candidate-prefix restart_plm_np8 \
  --checkpoint amr_reactive_plm_3d.chk \
  --reconstruction characteristic_plm --limiter mc
```

Coarse and fine state and temperature are rank-local throughout stepping.
Fine ownership follows the parent coarse slab, so some ranks may own no fine
cells. Checkpoint and CSV access is selected-root-only; the checkpoint payload
contains no owner map and may therefore be resumed with a different valid
rank count. Root gather/scatter is compatibility I/O rather than scalable
distributed output.

The fixed full-H2/O2 chemistry boundary uses paired inert and reacting
hotspots:

```sh
./pelef_amr_reactive_3d \
  ../cases/amr_reactive_3d/amr_hotspot_control.nml
./pelef_amr_reactive_3d \
  ../cases/amr_reactive_3d/amr_hotspot_chemistry.nml
mpiexec -n 8 ./pelef_mpi_amr_reactive_3d \
  ../cases/amr_reactive_3d/amr_hotspot_chemistry.nml chemistry_np8
python3 ../tools/check_amr_reactive_chemistry_3d.py \
  --control-prefix amr_hotspot_control_3d \
  --reactive-prefix amr_hotspot_chemistry_3d \
  --nx 8 --ny 8 --nz 8 \
  --i-lower 3 --i-upper 6 --j-lower 3 --j-upper 6 \
  --k-lower 3 --k-upper 6 --ratio 2
python3 ../tools/compare_mpi_amr_reactive_3d.py \
  --reference-prefix amr_hotspot_chemistry_3d \
  --candidate-prefix chemistry_np8 --final-time 1.0e-9 \
  --reconstruction characteristic_plm --limiter mc
```

Each chemistry phase advances both levels privately, averages down, and
recovers the covered coarse temperature. The application gates Euler and
H/O/N conservation separately from reacting species activity. A rank with no
fine cells still participates in reaction/integrator and success consensus.
Fixed-mechanism reacting checkpoint/restart now stores its reaction and
source-integrator context; see the fixed chemistry restart section below.

This case qualifies serial and rank-local sparse, static, strictly interior,
single-patch PCM or characteristic-PLM hydrodynamic AMR with optional fixed
elementary/full-H2/O2 chemistry and periodic coarse boundaries. Selected
mechanisms, chemistry restart, molecular transport, CTU/PPM, dynamic
regridding, scalable I/O, embedded boundaries, and boundary-touching patches
remain outside this case.


## Temperature-driven regridding (serial and MPI)

Set `regrid_interval_steps>0` in `&amr_reactive_3d` to rebuild the single
ratio-preserving interior patch at synchronized coarse steps, including step
zero. Cells at or above `regrid_temperature_threshold` (kelvin) define a
bounding box, enlarged by `regrid_buffer_cells` and clipped to the strictly
interior coarse cells. Empty tags retain the current patch; complete removal
of the fine level is not implemented.

The previous fine solution is averaged into the coarse grid before remapping.
Overlapping fine cells are copied exactly. Newly refined cells use constant
prolongation for PCM, or conservative component-wise limited-linear prolongation
for characteristic PLM. Density slopes equal the sum of species slopes;
a common per-parent slope factor is reduced until every child is nonnegative
and thermodynamically admissible. The symmetric child offsets preserve parent
averages. This transfer is not a new characteristic reconstruction scheme.
Coarse/fine temperatures are recovered from conserved state and the whole
remap commits only after the composite conservation check passes.

`dynamic_hotspot.nml` shows regridding during advection. Both serial front ends
use this path; chemistry/transport operators continue to use the current patch.
Dynamic checkpoint/restart uses schemas 6--9 and preserves the current patch,
completed regrid count, tagging interval/threshold/buffer and initial patch
configuration, alongside the existing physics context. Restart commits only
after geometry, state, EOS and synchronization checks pass. Checkpoints are
written after an accepted coarse step; restart resumes tagging at that stored
step. Static schemas 2--5 remain unchanged. Fixed-mechanism reacting checkpoints
use separate schemas; selected-mechanism chemistry retains its existing context.
MPI uses the same tagging/remap operator. At tagging steps, it gathers the
hierarchy on rank zero; after a patch change it rebuilds sparse slab ownership
and redistributes the remapped state. Ranks with no fine cells remain supported.
Checkpoint geometry is broadcast before restart distribution, so rank counts
can change. Fluid/chemistry/transport stepping remains rank-local. Regridding
and checkpoint I/O still require root memory for the complete hierarchy; this
is not a scalable distributed-remap implementation.
See [serial restart evidence](../../docs/validation/dynamic-amr-restart-20260915.md)
and [MPI regrid/restart evidence](../../docs/validation/mpi-dynamic-amr-20260915.md).
Multipatch/multilevel regridding and 3D EB coupling
remain open. The default interval zero preserves static runs.

## Physical boundaries

`physical_wall.nml` runs a strictly interior two-level patch with a 600 K
isothermal wall moving tangentially at 20 m/s and an opposite outflow.
The serial and sparse MPI fixed/selected frontends accept the regular `reactive_3d`
`boundary_faces`, reservoir inflow composition/temperature/pressure/velocity,
and wall temperature/velocity settings. For example:

```sh
mpiexec -n 4 ./pelef_mpi_amr_reactive_3d ../cases/amr_reactive_3d/physical_wall.nml
```

Physical runs report `Conserved boundary change:` in conserved-variable order,
positive into the composite hierarchy. Conservation diagnostics subtract this
cumulative surface transfer. Inlet wave speeds, inlet diffusion and isothermal
wall heat transfer also constrain the timestep. Static fine patches may touch
physical faces. Temperature-tagged movement/resizing can include physical faces with the
explicit setting below.

Checkpoint schemas 10–13 (static) and 14–17 (dynamic) extend the existing schema
with a boundary-settings signature and cumulative transfer. Restart rejects
changed faces, wall settings or reservoir inputs. Periodic schemas 2–9 retain
their layout. Fixed-mechanism reacting restart is also supported;
selected-mechanism context handling is shared, but the dated physical-boundary
application checks exercised only the fixed H2/O2 frontends. MPI checkpoint
restart is rank-neutral; four-to-two-rank physical restart is compared to serial
CSV output and retains the global boundary ledger.

See [serial evidence](../../docs/validation/amr-physical-app-20260917.md) and
[MPI evidence](../../docs/validation/amr-mpi-physical-app-20260917.md).

## Fixed chemistry restart

Serial and sparse MPI AMR support checkpoint/restart with `chemistry_enabled =
.true.` for both `elementary` and `full_h2o2`. The fixed mechanism remains fixed:
restart compares the saved reactions with the compiled model, alongside initial
composition, chemistry tolerances, the full-H2/O2 source SHA and the default integrator (elementary explicit,
full H2/O2 implicit). Species thermodynamics and transport retain their existing
checks. A selected bundle cannot be substituted through this path.

Reacting fixed schemas are the corresponding nonreacting fixed schema plus 16:
18/20 for periodic static hydro/transport, 22/24 for periodic dynamic cases,
26/28 for physical static cases and 30/32 for physical dynamic cases. The old
nonreacting and selected layouts remain unchanged. Reads commit only after all
context, state and geometry checks pass, including the physical boundary ledger.
See [focused restart evidence](../../docs/validation/amr-fixed-chemistry-restart-20260917.md).

## Boundary-touching static patches

The serial and sparse MPI fixed/selected AMR frontends accept fine patches that reach
physical domain faces. For example:

```sh
./build/debug/pelef_amr_reactive_3d cases/amr_reactive_3d/physical_touching_wall.nml
```

This short reacting full-H2/O2 case refines the cold moving x-low wall. Fine
wall heat/work fluxes enter the composite boundary ledger; wall thermal timestep
limits account for both coarse and fine states. Checkpoint/restart retains the
physical patch geometry and boundary/chemistry context with the existing static
schemas. Periodic seams remain unsupported. Boundary-touching dynamic regridding is
available in serial and sparse MPI with the explicit setting below.

[Application restart evidence](../../docs/validation/amr-boundary-touching-app-20260917.md)
covers a full-H2/O2 reacting wall case and an elementary reacting corner case.
The [MPI restart check](../../docs/validation/amr-mpi-boundary-touching-app-20260917.md)
compares four-to-two-rank restart with serial output for these cases.
Use `mpiexec -n 4 ./pelef_mpi_amr_reactive_3d` with the same static input.
The selected frontend shares this implementation; it has not been independently
qualified for boundary-touching application inputs in that check.

## Physical-boundary regridding

Set `regrid_physical_boundaries = .true.` together with a positive
`regrid_interval_steps` to let temperature tags and their buffer reach physical
faces. Periodic axes retain an interior-cell margin. The example
`dynamic_touching_wall.nml` expands the wall-touching patch at step zero and can
be run with the same serial or MPI executable as the static example.

Remapping preserves overlapping fine states, restricts derefined cells, and
uses physical coarse ghosts for conservative limited prolongation of new cells.
Checkpoint context includes this tagging mode, the moved patch and regrid count.
The default remains interior tagging, preserving existing serial/MPI cases.
MPI compares the tagging policy across ranks, performs the conservative remap
on rank zero, then rebuilds ownership and redistributes the new hierarchy.
Checkpoint/restart retains the tagging policy and supports changed rank counts.
This still requires root memory for the complete hierarchy.
See [serial evidence](../../docs/validation/amr-boundary-regrid-20260917.md) and
[MPI dynamic-boundary evidence](../../docs/validation/amr-mpi-boundary-regrid-20260917.md).
