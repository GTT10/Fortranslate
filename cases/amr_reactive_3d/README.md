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
cells. CSV rows are written directly from local slabs with collective MPI file I/O,
using the same byte layout as the serial writer. Final inventory, synchronization
and composition diagnostics use reductions without gathering fields. Checkpoint
writes also use local slabs, preserving the existing file layout after collective
metadata, temperature and level-synchronization checks. The payload contains no
owner map and may therefore be resumed with a different valid rank count. Fresh
runs initialize analytic conditions directly on each local slab and reduce
initial inventories. Restart reads the canonical fixed-width cell rows
directly into local candidates and commits only after collective EOS, synchronization
and metadata checks. Metadata is replicated; no full hierarchy is gathered for I/O.
Reformatted list-directed files use a bounded-memory streaming fallback on each
rank; that compatibility path scans all rows and is not a scalable I/O path. See the
[checkpoint input checks](../../docs/validation/amr-mpi-checkpoint-input-20260917.md),
[checkpoint output checks](../../docs/validation/amr-mpi-checkpoint-output-20260917.md) and
[distributed CSV checks](../../docs/validation/amr-mpi-distributed-csv-20260917.md).

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
for characteristic PLM and characteristic PPM. Density slopes equal the sum of species slopes;
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
MPI scans temperature tags in local slabs and reduces the global tag bounds.
When the patch changes, each rank restricts its old fine cells, exchanges
coarse halos, prolongs new children and reuses local fine overlap. Coarse slab
ownership stays fixed; fine ownership changes with the patch. No full-level
field gathering is needed for remapping. Ranks with no fine cells remain supported.
Checkpoint geometry is broadcast before restart distribution, so rank counts
can change. Fluid/chemistry/transport stepping and remapping remain rank-local.
Checkpoint reads and writes now operate on local slabs.
The library also provides `repartition_mpi_amr_sparse_3d` for explicit contiguous
coarse-slab counts (at least one coarse column per rank). It transfers both
levels and temperatures directly with MPI all-to-all exchange; subsequent
remaps preserve those counts. By default the application chooses equal coarse
column counts. Optional automatic balancing is described below; production
scaling remains unqualified.
See [repartition evidence](../../docs/validation/amr-mpi-repartition-20260917.md).
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
MPI compares the tagging policy across ranks and computes global tag bounds
from local temperature slabs. Changed patches are conservatively remapped
within each fixed coarse slab, with two-layer neighbor exchange for slopes.
Fine overlap stays on the same rank because its coarse parent ownership is fixed.
Checkpoint/restart retains the tagging policy and supports changed rank counts.
Checkpoint reads and writes now operate on local slabs.
See [serial evidence](../../docs/validation/amr-boundary-regrid-20260917.md) and
[MPI dynamic-boundary evidence](../../docs/validation/amr-mpi-boundary-regrid-20260917.md),
with the newer [distributed-remap check](../../docs/validation/amr-mpi-distributed-remap-20260917.md).

## MPI cell-count balancing

Set `load_balance = .true.` in `&amr_reactive_3d` to rebalance after
initialization, restart and a patch geometry change. The default is `.false.`.
The example `balanced_touching_wall.nml` can be run with four MPI ranks.

The partitioner counts coarse cells plus fine children in each coarse x-column,
minimizes the maximum count over contiguous rank slabs, and retains at least
one coarse column per rank. It transfers data only when that maximum decreases;
ties keep the current ownership. `Balanced coarse columns:` reports each rank's
new coarse width. Chemistry stiffness, subcycle counts and measured runtime
are not included in this cell-count estimate.

The serial executable accepts the same input and performs no repartitioning.
The option does not change checkpoint physics/context records; the MPI layout
is recomputed for the current rank count on restart. Existing global state
and boundary-ledger comparisons remain applicable. See
[balancing and restart evidence](../../docs/validation/amr-mpi-load-balance-20260917.md).

## Pressure outlet

The serial and MPI fixed/selected two-level applications accept `outlet_pressure(6)`
in `reactive_3d`, using the regular six-face order and Pa units. Positive finite
values up to 1e15 require `outflow`; zero keeps extrapolation. Targets feed hydro,
molecular transport and coarse interpolation ghosts during regridding. Existing
`regrid_physical_boundaries` policy still controls boundary-touching regrids.

Checkpoint boundary identity includes all six targets when any is positive.
Restart rejects changed or removed targets. MPI restart can change rank count;
see [MPI evidence](../../docs/validation/pressure-outlet-mpi-application-20260918.md).
Patch-tree pressure input remains unconnected and is rejected by its readers. See
[application/restart evidence](../../docs/validation/pressure-outlet-amr-application-20260918.md).

## Two-level PPM applications

The fixed and selected serial/MPI two-level front ends accept
`reconstruction="characteristic_ppm"`. This connects conservative CTU hydro,
reaction/transport splitting, dynamic limited-linear remapping and checkpoint
restart. The stored reconstruction must match on restart. Sparse MPI uses four-layer halos and the same fine CTU update with rank-local
storage. With `riemann_solver="pelec"`, CTU evolves transverse primitive pressure
separately and uses multidirectional shock flattening. Other solvers retain the
conservative/EOS correction. The stored solver must also match on restart.
Multi-patch tree MPI PPM is a separate path.
See [the focused application check](../../docs/validation/ppm-two-level-application-20260925.md)
for the reacting, transporting physical-boundary/regrid fixture and its limits.
The [MPI application check](../../docs/validation/ppm-mpi-two-level-application-20260925.md)
records four-to-two-rank restart against serial output.

### Artificial-viscosity coefficient

Fixed/selected serial and MPI two-level applications accept `difmag` in
`&reactive_3d`. Default `-1` selects 0.1 for PPM/pelec and zero otherwise.
Explicit nonnegative input requires PPM; zero disables the correction.
The resolved value reaches coarse and fine substeps and is included in the
checkpoint reconstruction identity. Automatic and equivalent explicit values
can resume each other; changed coefficients reject restart. Older PPM
checkpoints without a coefficient identity are incompatible. PCM/PLM records
retain their existing reconstruction identity. See
[input/restart evidence](../../docs/validation/difmag-two-level-input-20260925.md).
