# Rotating-energy serial transport — 2026-09-27

Source `3aa1b3615e95f95f1f3dad3ce1429de65ab65fb7`, main, base `2462089b`.

## Change

Serial MOL SDC now accepts viscosity, conduction, species diffusion and
barodiffusion with rotation, including embedded-wall heat and no-slip traction.
Conserved halos are filled in rotating energy first. Private state/halo copies
then restore ordinary relative total energy at each receiving cell before
entering existing transport kernels. Covered cells and unfilled halo slots are
excluded from conversion. The public conserved tree remains in rotating energy.

Transport has zero total-mass flux, so its energy flux needs no `Phi` correction.
The same private-input conversion applies to the embedded-wall stencil. The
previous transport-with-frame rejection has been removed; CTU with rotation
remains unsupported, consistent with the pinned rotating MOL path.

## Completed evidence

GNU Fortran 13.3.0 Debug fixed core; isolated heptaneSimple selected serial build.
Existing configuration: MPI/spray ON, CTest targets/SUNDIALS OFF. Manual tests
use `-g -fcheck=all -ffree-line-length-none`, compile from the artifact directory
(to avoid old root-level module files), and unset both library-path environment
variables. Selected serial application was built, not run.

- New `test_rotating_transport_3d` passed four separate viscosity/conduction/
  diffusion/barodiffusion comparisons on a 4³ all-fluid single patch. Initial
  density, velocity, pressure and H2/N2 fractions vary across cells, O2=.2;
  pressure is about 300 kPa. Rotation is 100 rad/s about z through `(10,10,0) m`
  and rotating energy is negative. Each operator produces a nonzero rate.
  Physical-domain diffusive faces are closed in this comparison, so the same
  physical interior states determine both ordinary and rotating results.
- Rates and time-integrated face fluxes matched exactly (gate `3e-11`);
  limiters matched exactly and total-mass rates were exactly zero.
- A cut plane at x=.31 m with 400 K embedded wall exercised rotating heat and
  tangential traction through a two-cell halo. Integrated energy rate was
  `-180.74991110138967 W`, y momentum rate `-8.9529447823897813e-5 N`;
  wall frequency was positive and source/ledger integration residual was zero.
  This is a source evaluation, not an evolved rotating wall trajectory.
- The existing rotating tree test now passes eight modes. Two added modes
  enable reaction and all transport flags with two SDC iterations, one root
  step of `1e-9 s`, two levels and three patches. Common/acoustic and
  recursive/LF composite relative residuals are `3.5379380752486782e-15` and
  `4.1698297793315234e-15`, below the unchanged `3e-14` gate. This tree fixture
  is pure N2; chemical conversion and nonzero species diffusion are established
  by their separate focused cell/operator tests, not by these coupled cases.
  Combined tree executable reports 6,227 assertions, mostly prior cell checks.
- Existing `test_serial_sdc_species_3d` passed four ordinary CTU/MOL acceptance
  and rollback cases. Fixed/selected builds, provenance, artifact Python syntax,
  project contract and diff whitespace checks passed. No test failed in this
  increment; existing executable-stack linker warnings remain.

## Artifacts and remaining scope

Exact commands and logs: `build/rotating-transport/commands.txt` and neighboring
files. No new compiled external-solver comparison, full CTest matrix, installed
binary audit, rotating wall trajectory or physical validation was performed.
MPI transfer/stepping, conservative remapping, timestep selection and application/
checkpoint frame identity remain open, as does the earlier 1000 rad/s LF
absolute patch-conservation condition.

Executable SHA-256:

- check: `a92a1343d6e610c728f1657ed85760d99f82c5e74d5231e7ad5d674b4801b7e6`
- tree: `c90a6aa67780a0b7e6e3fdb032f9e8297daec0ec1b791dcf3b969ccde7a20b87`
- default: `ed5f18308b48b58c303019c281c8ed7cb6cb0130e22477d87cc1e614d78433d6`
- selected serial: `8db52325c18155cc1f6079405b03f91b2624e097991139d7b56a43f363b80945`
