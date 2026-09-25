# Serial periodic EB tree transport and recursive updates — 2026-09-26

Source `37295c9`, parent `4b27b85`, local main. Common molecular transport
now forwards periodic labels into the tree flux synchronization path. The
shared recursive hydro/transport operator also forwards these labels and
corrects periodic provisional endpoint states before descendants interpolate
their ancestor data in time. High-domain face owners supply the canonical
transfer; owner flux copies avoid aliasing for domain-spanning patches.
Final synchronization includes restricted fine transfers across the seam.
The existing private candidate transaction is retained.

## Completed checks

GNU Fortran 13.3.0, Ninja Debug, Python 3.12.3; existing tests-enabled build.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_eb_tree_transport_3d test_eb_tree_subcycle_transport_3d test_eb_tree_subcycle_hydro_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_eb_tree_(transport|subcycle_transport|subcycle_hydro)_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**3/3 passed**, no failures or skips, 16.06 s total:

- Common transport: 1.28 s; viscosity, conduction, species diffusion and
  barodiffusion enabled.
- Recursive PCM hydro: 3.04 s.
- Recursive transport: 11.69 s, including diffusion-scaled substeps.

Each adds periodic cases for three plane normals and four/five-patch trees.
The largest normalized inventory residual across the 18 periodic updates was
5.45e-15, below the unchanged 3e-14 bound. Transport limiter values remain
within [0,1]. Unpaired periodic transport labels reject while preserving tree
states, temperatures and diagnostic outputs. Existing invalid interval/deep
state rollback and nonperiodic cases remain in the tests.
Provenance, Python syntax, project contract and whitespace checks passed.

Artifacts: `build/eb-tree-periodic-recursive-checks/build.log` and `LastTest.log`.
Mechanism SHA:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

Recursive PLM/CTU periodic stepping and full reacting splits were not exercised
by these three tests. MPI periodic time updates and application input remain
pending. No full configuration matrix, installed audit, external solver parity
or physical validation was run.
