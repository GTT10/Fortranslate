# MOL face-centroid stencil — 2026-09-27

## Reference and implementation

Implementation `2d113ac520599ffdd4b20e63e46b145bf2ebcc8a` separates the MOL stencil from the default CTU
connectivity-filtered interpolation. Pinned PeleC
`bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`, locally exported under
`build/pelec-pinned-hit-source`, uses `Source/EB.cpp` functions
`pc_fill_flux_interp_stencil` and `pc_apply_face_stencil` from
`Source/Diffusion.cpp` for MOL. This path uses absolute centroid offsets
without the CTU connectivity mask, and skips aperture-weighted coefficients
whose magnitude is at most `1e-14`. The Fortran MOL caller now explicitly
selects those rules. It retains intensive fluxes; divergence applies aperture
once. Default CTU callers retain their existing rules.

## Completed checks

GNU Fortran 13.3.0/Open MPI 4.1.6, existing Debug selected-heptane Simple build;
MPI/spray enabled, tests/SUNDIALS disabled. Commands:

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/heptane-simple-mpi-debug --target pelef_core pelef_mpi_reactive_eb_tree_3d_selected_spray --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH /usr/bin/gfortran -g -fcheck=all -ffree-line-length-none -I build/heptane-simple-mpi-debug/modules tests/unit/test_eb_reactive_wall_flux_3d.F90 build/heptane-simple-mpi-debug/libpelef_core.a -o build/mol-centroid-unit/check
env -u LD_LIBRARY_PATH -u LIBRARY_PATH build/mol-centroid-unit/check
python3 build/run_mol_centroid_acoustic.py
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
python3 -m py_compile build/run_mol_centroid_acoustic.py build/run_mol_lf_redist.py
git diff --check
```

Build and unit pass, including four added checks for geometry, analytic
bilinear weights, the small-coefficient cutoff and transactional rollback.
Logs: `build/mol-centroid-build.log`, `build/mol-centroid-unit/run.log`.
Six nested acoustic-MOL particle launches pass with exact 1→2 rank gas/particle
restart and bidirectional MOL/CTU mode mismatch rejection. Maximum absolute
ledger residual `2.710505431213761e-20` remains below the unchanged `1e-12`
gate. Result: `build/mol-centroid-acoustic/result.json`. Binary SHA256:
`7274a29ee263f5508c7f30b3ea0bbf92033662947f9452760ea2082ea4e84b4b`.
This is focused analytic/application evidence, not an upstream full-case run.

## LF negative-fuel diagnosis

Before this repair, read-only GDB extracted the failing spatial evaluation
from the binary recorded in [the LF application check](mol-lf-application-20260927.md).
Artifact scripts/logs: `build/mol-lf-rates.gdb`, `build/mol-lf-rates2.log`,
`build/mol-lf-centroid.gdb`, `build/mol-lf-centroid-final.log`.
The final centroid breakpoint script uses a Fortran logical condition; an
initial C-style condition failed and stopped on an unrelated first call.
Only the final log is used here.

At fine cell (2,4,4), initially zero fuel, all six un-interpolated fuel face
fluxes are exactly zero. The physical y/z face centroids have x offset 0.2,
and each borrows 20 percent of the adjacent x-cell flux. That neighbor's
lower/upper y fluxes are -/+`1.5933196794018477e-5`; its z fluxes are
`-1.5933196758775141e-5` and `1.5933196829261811e-5`.
The resulting hydro fuel rate is `-0.025493114870429559`. Particle and
external rates in the failed cell are zero; transport changes the combined
rate to `-0.025486868635953324`. A 50 ns fine step therefore proposes fuel
partial density `-1.2743434317976661e-9` and is rejected.

This locates the negative contribution in centroid interpolation of neighboring
LF diffusion, rather than evaporation or a roundoff-only centroid offset.
The pinned MOL weights produce the same interpolation for this plane;
separating MOL from CTU stencil semantics does not remove this mechanism.
The state-redistribution variant (`build/run_mol_lf_redist.py`,
`build/mol-lf-redist.log`) also fails its first root step. No positivity fix,
clipping, tolerance change or upstream transient-failure parity is claimed.

After the stencil repair, `python3 build/run_mol_centroid_lf.py` again exits 1
on its first nested LF root step (`build/mol-centroid-lf-failure.log`).
The negative-fuel case remains unresolved.
