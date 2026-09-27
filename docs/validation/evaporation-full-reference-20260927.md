# Full evaporation reference diagnostic — 2026-09-27

Implementation: `a9cc5569`; candidate binary built from `be25365087344799b8798f8ccd0aadddd9b0bee7`.
The fixture now explicitly selects the existing `pelephysics` drag law.
Abramzon–Sirignano heat/mass transfer was already enabled. No numerical kernel
or acceptance tolerance was changed.

The pinned PeleC Spray-Evaporation executable completed two 1e-7 s steps
with two SDC iterations on the prepared periodic 4³, 400 K gas / 310 K liquid
adaptation. This is not the original upstream 273 K case.
Source: PeleC `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`, PelePhysics
`6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`; complete dependency identities
are in `build/pelec-pinned-hit-source/pelef-reference-provenance.json`.
Fresh reference build: GNU 13.3.0 Release, serial, particles ON, tests OFF,
`build/pelec-pinned-spray`; configure/build logs use that path prefix.
Candidate configuration is recorded in the
[SDC integration record](nonreacting-particle-sdc-20260927.md).

## Results

All four parcels have the same values. Changes are measured from the shared
initial state; reporting relative final temperatures would hide this discrepancy.

| Quantity | Reference change | Fortran change | Change-magnitude difference |
| --- | ---: | ---: | ---: |
| Temperature (K) | 0.004727062899 | 0.003210010595 | -32.09% |
| Single-drop mass loss (kg) | 1.133273583e-16 | 1.289211782e-16 | +13.76% |
| Axial speed reduction (m/s) | 0.000697978230 | 0.000615191466 | -11.86% |

**Particle trajectory parity fails to follow from these runs.** Matching the
drag model alone does not resolve the differences. Film transport properties
and integration remain candidates for diagnosis, not established causes.
No gas-field comparison, CTest matrix, multirank or physical validation was run.

## Reproduction and artifacts

Reference executable:
`build/pelec-pinned-spray/Exec/RegTests/Spray-Evaporation/PeleC-Spray-Evaporation`.
It was run with `reference.inp` from `build/spray-evaporation-sdc-run2`, exit 0;
log `reference.log`, output `spray00002.p3d` in that directory.

```sh
python3 tools/prepare_evaporation_reference.py --bundle build/heptane-3sp-reference/heptane_3sp.json --output build/spray-evaporation-matched-drag
# From the prepared directory:
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none /usr/bin/mpiexec -n 1 /home/gtt13/projects/pelec/build/heptane-spray-mpi-debug/pelef_mpi_reactive_eb_tree_3d_selected_spray tree.nml spray.nml
# From the repository root:
python3 tools/compare_evaporation_particles.py --reference build/spray-evaporation-sdc-run2/spray00002.p3d --candidate build/spray-evaporation-matched-drag/fortran.csv.particles.rank0.csv --output build/spray-evaporation-matched-drag/comparison.json
```

One reference run and one matched-drag candidate run completed with exit 0.
The diagnostic records input hashes, CGS-to-SI conversions and per-ID changes;
successful execution does not assert parity. Mechanism provenance, Python syntax,
project contract and `git diff --check` passed. Existing earlier candidate
artifacts retain the default Schiller–Naumann drag setting.
