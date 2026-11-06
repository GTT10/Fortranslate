# Velocity-consistent halo controls with planar cut cells, 2026-10-22

The pressure/velocity diagnostic now accepts `--smooth-cut-contact` and
`--imported-state` as its sixth argument. The former initializes the same
900 K H2/N2 composition contact with x velocity 10 m/s on the imported planar
z geometry. The variation and velocity are tangential to the wall, so the
analytic x cell averages also give exact fluid-volume averages in cut cells.
Solid states are excluded and the same restriction-consistency gate applies.
The latter preserves the complete imported field, including its temperature
and composition gradients. The optional seventh `--consistent-velocity`
argument remains diagnostic-only.

The candidate now applies planner-provided reflection signs to momentum
stencil samples and reflects the reconstructed normal momentum at physical
walls. Valid leaf fields, density and all species are checked for exact
immutability during the ghost-only pressure/velocity repair. Incoming smooth
halo diagnostics skip solid cells. No production application policy changed.

`tools/run_eb_velocity_controls.py` freezes its input and executable into a
new output directory, runs pressure-only and pressure-plus-velocity policies
for both profiles, and checks accepted inventories, cell coverage and all
intermediate synchronization stages. It also requires identical baseline
histories with/without the candidate flag and bounds the smooth cut contact
pressure error by 1e-8 Pa. Failures are retained in record.json; output
directories cannot be reused by this runner.

Artifacts: `build/eb-velocity-cut-20261022/`, including input.dat, frozen-probe,
record.json, source-identity.json, four run logs, CSVs and checker reports.
The baseline mode in each run uses the existing separate-momentum policy;
the paired mode uses pressure reconstruction with/without the velocity option.

## Completed results

All four runs (eight baseline/paired histories) complete four 50 ns steps
through 200 ns. Every inventory, cell and stage checker passes. For each
profile the baseline histories are byte-identical with/without the velocity
flag. The smooth cut initializer restriction residual is 4.40626e-16.

| Paired diagnostic | Pressure only at 200 ns (Pa) | Plus velocity consistency (Pa) |
| --- | ---: | ---: |
| Smooth cut contact | 0.001517145196 | 2.91038305e-11 |
| Imported cut field | 0.001246740285 | 0.001246740358 |

The smooth cut candidate passes the 1e-8 Pa equilibrium bound at every
accepted clock. Its worst independently checked inventory residual is
9.78541e-15 and cumulative absolute reset correction is 3.91492e-17 J,
versus 1.25361e-11 J with pressure reconstruction alone.

For the imported candidate, worst inventory error is 9.07818e-15 and
cumulative absolute reset correction is 2.12981e-11 J, essentially unchanged
from 2.13014e-11 J without velocity consistency. Pressure changes by only
7.28e-11 Pa. Thus the velocity correction preserves the isothermal composition
contact even with planar cut cells, but does not remove the residual in the
nonuniform-temperature imported control. The latter must not be relabeled as
a solved pressure-equilibrium case. Temperature/internal-energy reconstruction
and nonlinear coarse/fine mixing remain distinct discriminators before
production adoption; the current result does not identify their relative
contributions.

## Verification scope

Base: dirty main at `0a0e00ebe9138145fbb0e0e745665d0d768df009`.
GNU Fortran 13.3.0 Debug; existing isolated OpenMPI registration build.
No compiler or MPI configuration change. The existing executable-stack warning
remains on the reference probe.

Commands:

- `cmake --build build/eb-eos-application-20261021/isolated-openmpi/registration --target test_eb_tree_reacting_reference --parallel 2`
- `python3 tools/run_eb_velocity_controls.py --probe build/eb-eos-application-20261021/isolated-openmpi/registration/tests/test_eb_tree_reacting_reference --input build/eb-pressure-reset-20261022/moving50.dat --output build/eb-velocity-cut-20261022`
- `ctest --test-dir build/eb-eos-application-20261021/isolated-openmpi/registration -R '^unit_eb_reactive_transfer_3d$' --output-on-failure`: 1 passed.
- `python3 -m unittest discover -s tests/python -p 'test_eb_pressure*.py'`: 19 passed, zero failed/skipped.
- Mechanism provenance, Python syntax, project contract and whitespace checks pass.

Full held-MOL integration must propagate the selected policy through prepared
hydro/transport halos, predictor/corrector stages, retries and history replay,
then mirror the policy in MPI and restart validation. This diagnostic does
not establish that integration, arbitrary EB geometry, native-reference
parity or physical validation.
