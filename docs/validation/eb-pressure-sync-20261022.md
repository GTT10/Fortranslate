# Pressure changes across coarse/fine synchronization, 2026-10-22

The diagnostic probe writes `*_stages.csv` on composite fluid leaves before
and after each gas+held reflux operation. The independent checker compares
the same cells, checks finite data, duplicate/missing cells and all expected
Euler stage numbers. These are intermediate SSPRK2 states, not accepted-step
physical clocks. Accepted-step inventories remain checked independently.
The pre-sync spatial operator already uses AMR halos; it is not a standalone
single-level advection control.

## Imported regular geometry control

The regular restricted input from the preceding localization experiment was
rerun for four 50 ns steps with the existing accepted-step dual-energy reset.
The baseline and paired accepted histories are byte-identical to the preceding
run, verifying that recording stage pressures did not change those results.
Both inventory checks and all eight stage-coverage checks pass.

At the first Euler stage, paired coarse maximum pressure offset increases
from 0.00071613499 Pa before synchronization to 0.00287121785 Pa afterward.
The largest cell change is +0.00215508285 Pa at coarse (5,1,3), adjacent to
the covered region. Fine pressures are unchanged by synchronization.
Thus synchronization contributes directly to the remaining coarse pressure
error in this imported diagnostic. It is not the only source: a nonzero
pressure offset is already present before synchronization.

This is not evidence to disable conservative reflux. The baseline coarse
maximum decreases from 0.0638300002 to 0.00287121785 Pa in the same stage.
The paired policy changes the pressure representation used by the spatial
operator, while accepted gas and held quantities still undergo conservative
flux replacement.

## Smooth composition control

`--smooth-contact` is an explicit sixth diagnostic argument; it requires
regular geometry and the pressure-evolution path. Imported states are
replaced only in this opt-in diagnostic. It initializes H2/N2 mole fractions
with a periodic cosine, X(H2)=0.2+0.05*cos(2*pi*x/L), using the exact cell
average cosine (sinc factor). Pressure is 101325 Pa, temperature 900 K,
and x velocity 10 m/s. All other velocities and species fractions are zero.
At fixed T and p, conserved species and energy densities are affine in mole
fractions; coarse cell records therefore agree with volume-averaged children.
The probe verifies restriction with the unchanged 1e-12 relative threshold.
This is a hydro-only contact, not a reacting H2/N2 validation case.

The first attempted control used 1000 K, exactly the NASA polynomial
breakpoint. Recovered temperatures ranged from 999.99999999999989 to
1000.0000000000001 K and restriction failed with residual 9.9523e-9.
That failed attempt is retained, not counted as a passing run. The diagnostic
uses 900 K to avoid crossing the coefficient boundary; no tolerance or
mechanism data was changed.

The 900 K initializer passes restriction with residual 4.40626e-16; recovered
coarse temperatures are 899.99999999999977 to 900 K. Both baseline and paired
runs complete four steps and all budget, localization and stage checks pass.
The worst independently checked inventory residual is 2.44233e-14.

For the paired control, stage 1 remains at pressure roundoff (coarse
1.46e-11 Pa, fine 2.91e-11 Pa) and synchronization changes no leaf pressure.
At stage 2, coarse pressure is still at 1.46e-11 Pa before synchronization,
but rises to 9.66370e-5 Pa afterward, at coarse (8,1,1), across the periodic
coarse/fine interface. Fine pressure already reaches 1.93274e-4 Pa before
synchronization and is unchanged by it. The first accepted 50 ns endpoint
has maximum offset 9.66372e-5 Pa. At 200 ns, baseline and paired maxima are
0.0015504578 and 0.0015502958 Pa, respectively. The accepted E-H-KE gap is
zero with the explicit reset; the paired cumulative absolute numerical
energy correction is 4.77372e-11 J.

Thus an analytically consistent initial contact alone does not prevent
pressure errors at later stages. The first non-roundoff fine signal arises
in the second spatial update, and synchronization transfers an error to
coarse interface cells. This narrows the next implementation target to
stage-updated material/held-energy state and flux consistency. The paired
initial-pressure representation is not yet a sufficient evolving-contact
solution. A pure single-level reference and full held-MOL closure are still
needed to distinguish interior reconstruction from halo reconstruction.

## Artifacts and verification

Artifacts: `build/eb-pressure-sync-20261022/`; the completed imported run,
failed 1000 K initializer, separate frozen `smooth900/` retry, commands,
return codes, input/binary/source hashes, per-cell/stage CSVs and checker JSONs.
Base: dirty main at `0a0e00ebe9138145fbb0e0e745665d0d768df009`.
Toolchain: GNU Fortran 13.3.0 Debug, existing isolated OpenMPI registration.
The existing executable-stack linker warning persists.

Completed checks: focused probe build, mechanism provenance, Python syntax,
project contract and `git diff --check`; 19 pressure diagnostic unit tests
passed, zero failed/skipped. Commands:

- `cmake --build build/eb-eos-application-20261021/isolated-openmpi/registration --target test_eb_tree_reacting_reference --parallel 2`
- `python3 -m unittest discover -s tests/python -p 'test_eb_pressure*.py'`
- `python3 -m py_compile tools/check_eb_pressure_sync.py tests/python/test_eb_pressure_sync.py`
- `python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json`
- `python3 tools/check_project_contract.py --source .`
- `git diff --check`

Production numerical policy, full prepared-MOL closure, MPI, physical
convergence and fixed-reference qualification remain unchanged/open.
