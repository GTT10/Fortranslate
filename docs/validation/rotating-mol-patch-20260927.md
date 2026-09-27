# Rotating-energy EOS and single EB MOL patch — 2026-09-27

## Scope and source

Implementation `d1d0bf967c01f98c6e76dbc95eb03f6eec25c1bb`, on `main` after
`443eb6c6`. This connects the earlier [rotation kernels](rotating-frame-kernels-20260927.md)
to primitive/conserved conversion, EB slip-wall pressure and a single MOL
spatial patch. The energy density is `rho*e + rho*|v_relative|^2/2 - rho*Phi`.
EOS recovery adds `rho*Phi` before subtracting kinetic energy. Negative rotating
total energy is accepted when the recovered thermal state is valid.

The optional frame supplies cell-center potentials for both adjacent Riemann
states and the physical cut cells. Covered gas cells may retain NaNs.
Absent rotation keeps the existing ordinary-energy path. Checked EOS arithmetic
rejects invalid potentials; patch outputs are published only after success.
Rotational body force is **not** added by this spatial operator: old/new source
evaluation and time integration remain separate work.

## Completed checks

GNU Fortran 13.3.0, Debug fixed core, runtime checks on the manually linked tests.
The existing development build enables MPI and spray, disables CTest targets and
SUNDIALS; the selected serial target uses the three-species heptaneSimple bundle.
Builds and links unset `LD_LIBRARY_PATH` and `LIBRARY_PATH`.

- Fixed core and isolated selected serial application builds passed. Existing
  executable-stack linker warning remains. Selected application was not run.
- Extended `test_eb_mol_hydro_3d`: **16 cases passed** at default angular speed
  400 rad/s: ordinary/rotating × stationary/nonuniform × acoustic/LF × order 1/2.
  Inclined-wall geometry has 10³ cells including three halo layers, a central
  4³ physical patch, 0.25 m spacing and integration weight `dt=1e-7 s`.
  Rotation uses z axis through `(0.13,-0.07,0.2) m`; states use pure N2,
  density near 1 kg/m³ and pressure near 150 kPa.
- EOS checks cover exact zero-potential identity, energy subtraction, unchanged
  other components, primitive/temperature/sound-speed round trips at `1e-12`,
  negative rotating total energies, and invalid-potential rejection.
- Maximum absolute integrated conservation residual was `2.2204460492503131e-16`,
  below the existing `1e-15` gate. Stationary acoustic wall balance and ordinary
  LF balance passed `1e-7`; rotating LF stationary balance is not asserted,
  because its diffusion acts on spatially varying rotating energy.
- Active-cell NaN and invalid-frame failures preserved patch rates, face fluxes
  and both ledgers exactly. Existing `test_reactive_pelec_riemann` passed.
- Mechanism provenance, artifact Python syntax, project contract and
  `git diff --check` passed. No repository Python source changed.

## Failed additional condition

The initial 1000 rad/s fixture exceeded the unchanged absolute conservation gate.
The reproducible command `build/rotating-mol-patch/check 1000` exits 1 for the
rotating, nonuniform, second-order LF case: residual
`2.6645352591003757e-15` versus `1e-15`. The committed default fixture uses
400 rad/s and still exercises negative rotating total energies. This is a
restricted passing condition, **not** resolution of the 1000 rad/s failure.
Neither solver tolerances nor the original absolute test gate were loosened.
The failed log is retained separately; broader rotation qualification is open.

## Artifacts and limits

Exact commands: `build/rotating-mol-patch/commands.txt`. Build, default regression,
patch success and larger-speed failure logs live in that directory.
Executable SHA-256:

| Artifact | SHA-256 |
| --- | --- |
| `check` | `4cf3bf02010ea44b5549dcf37a090c69927493b5f94f4902d3d63040b51d467e` |
| `default` | `3727aac62d73f29ac0848ccef58a2fdcc7e4475f6355b710f360392cfccfdceb` |
| selected serial application | `b75d633662d04fdf9190deb258f73931ae0cc73970a1542179618c8c07882746` |

No new full reference-solver run, installed-binary audit, CTest matrix or physical
validation was performed. Tree/halo transfer, reaction and transport, AMR
restriction/reflux/remapping, source time integration, initialization, boundary
conditions and checkpoint frame identity still need integration before rotating
CFD can be enabled. The unrelated frozen Jet continuation remained running.
