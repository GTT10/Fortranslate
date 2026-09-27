# Serial AMR gas application LES SDC — 2026-09-27

Implementation `b1030ac5`, based on `fb15aa8b`, local main, version 0.246.0.

## Change

The opt-in serial fixed/selected AMR gas frontends now read `amr_tree_les` and
connect the serial instantaneous LES evaluator to common and recursive SDC.
The generic gas driver receives a callback, keeping its dependency independent
of the optional LES runtime. Parent temporal endpoints and stage fractions
reach LES evaluation. The timestep uses the same conservative diffusion limit
and current finest patch spacing as the MPI application.

LES configuration uses the same canonical checkpoint identity as MPI. Omitting
the group leaves the existing gas path in place. The callback requires SDC and
a supplied identity; checkpoint context retains its required final newline.

## Independently completed checks

Artifacts: `build/les-serial-amr-app-checks/` (`provenance.log`, `build.log`,
`build-final.log`, `ctest.log`, `LastTest.log`) and individual application
inputs/logs/CSVs in
`build/spray-amr-reactive-debug/les-serial-mpi-amr-sdc-app-{common,recursive}/`.
Existing fixed-mechanism Ninja Debug build, GNU Fortran 13.3.0 and Open MPI
4.1.6 (`/usr/bin/mpiexec`).

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_amr_tree_reactive_3d pelef_mpi_amr_tree_reactive_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^les_serial_mpi_amr_sdc_app_(common|recursive)$' --output-on-failure
python3 -m py_compile tools/check_les_mpi_amr_sdc_app.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**2/2 passed**, no failures/skips: common 35.81 s, recursive 31.48 s,
CTest total 67.34 s. Build, MPI ABI audit, mechanism provenance, Python syntax,
project contract and whitespace checks passed. Existing real-comparison and
executable-stack linker warnings remain. Mechanism source SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

Both tests extend the preceding two-level fixed-mechanism Smagorinsky fixture:
4-cubed root, moving no-slip walls, molecular transport and disabled chemistry.
Serial continuous gas and budget CSVs match MPI byte-for-byte. A checkpoint
saved on 2 MPI ranks continues identically in serial; a serial checkpoint
continues identically on 4 ranks. The existing LES activation, 2-to-4-rank
continuation and changed/removed-LES rejection checks also pass. Common and
recursive continuous cases advance 8 and 4 steps respectively.

## Limits

Selected frontend wiring is present but was not built/run in this configuration.
These application tests exercise Smagorinsky only; all four closures have
separate library evidence. Active chemistry, regridding, other boundary types,
installed binaries, full configuration matrices, external PeleC parity and
physical validation were not run. Particle-SDC LES routing remains open.
