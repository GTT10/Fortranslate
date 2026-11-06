# Moving pressure-error localization, 2026-10-22

The hydro-only SSPRK2 diagnostic now writes every accepted composite fluid
leaf to `*_evolution_cells.csv`: step, node/index, volume fraction, fluid
volume, signed pressure offset from 101325 Pa, and E-H-KE gap. Covered coarse
cells and solid cells are excluded. The independent localization tool checks
finite values, duplicate leaves, unchanged coverage/geometry, complete clocks,
and agreement of cell maxima with the existing history. It reports maxima
and volume-weighted mean absolute pressure offsets by level and cut/full cell.
Coverage consistency is checked across times, not against an independent mesh
enumerator; this is not a standalone topology validation.

## Matched diagnostic

Both inputs use x velocity 10 m/s, four 50 ns SSPRK2 steps (200 ns), converged
held EOS, separate momentum slopes and the existing accepted-step dual-energy
reset. Each runs both baseline and paired pressure policies. The regular
control removes the EB geometry and recomputes covered coarse conserved
records as equal-volume averages of eight fine records. It activates formerly
solid cells; it is a geometry ablation, not a new smooth physical initializer
or a spatial-convergence sequence. Initial composite pressure offsets are at
most 4.37e-11 Pa for both inputs.

| Geometry / policy | Maximum offset at 200 ns (Pa) | Peak node, i,j,k | Worst inventory relative error |
| --- | ---: | --- | ---: |
| Cut / baseline | 0.0612345946 | 2, 3,1,8 | 1.09e-14 |
| Cut / paired | 0.0012467403 | 1, 6,2,3 | 8.67e-15 |
| Regular restricted / baseline | 0.0615728743 | 2, 3,1,8 | 1.51e-14 |
| Regular restricted / paired | 0.0100544968 | 1, 5,2,3 | 1.26e-14 |

Node 1 is the coarse level; node 2 is the fine level. The original-policy
peak lies on the fine patch boundary. Paired pressure moves the peak to the
coarse level. In the cut paired case, coarse cut/full maxima are respectively
0.0012467403 / 0.0012467191 Pa; their volume-mean absolute offsets are
0.0005407361 / 0.0005405077 Pa. Thus the remaining pressure signal is not
confined to cut cells. The regular paired peak is adjacent to the covered
coarse region (i=5, with covered indices 1:4), and is larger than the cut
case. This does not establish that removing EB worsens a matched physical
solution: the ablation activates different fluid records and changes volumes.
It does establish that the residual does not require cut cells.

The next numerical discriminator is coarse interior advection versus
coarse/fine synchronization on a consistent smooth regular initializer.
Full prepared-MOL closure, reaction/transport, MPI, native-reference parity
and spatial convergence remain open. No production numerical policy changed.

## Reproducibility and completed checks

Artifacts: `build/eb-pressure-location-20261022/`, including frozen-probe,
run.py, two inputs, command/return-code/input-hash record.json, per-cell CSVs,
independent budget/localization JSONs and source-identity.json. The completion
event was consumed and the reports analyzed in the same active turn.

Base HEAD: `0a0e00ebe9138145fbb0e0e745665d0d768df009`, dirty main.
GNU Fortran 13.3.0 Debug; existing isolated OpenMPI registration build.
Binary SHA256: `645816449ad5b5eb577f810bad89610af071b54d81e4ee5df0cb8edf8a39927b`.

Completed commands:

- `cmake --build build/eb-eos-application-20261021/isolated-openmpi/registration --target test_eb_tree_reacting_reference --parallel 2`
- `python3 build/eb-pressure-location-20261022/run.py` (2 inputs, 4 policy histories; all budget and localization checks passed)
- `python3 -m unittest discover -s tests/python -p 'test_eb_pressure*.py'` (15 passed, 0 failed, 0 skipped)
- `python3 -m py_compile tools/localize_eb_pressure_evolution.py tests/python/test_eb_pressure_location.py`
- `python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json`
- `python3 tools/check_project_contract.py --source .`
- `git diff --check`

The initial project-contract invocation omitted required `--source` and exited
with an argument error; the corrected invocation passed. The build retains
the existing executable-stack linker warning. No full matrix or MPI run is
included in these counts.
