# Paired pressure representation for EB MOL

The opt-in `paired_pressure` policy supplies one pressure representation to
coarse covered spatial states and incoming fine halos. It requires
`consistent_velocity` and `separate_momentum`, held MOL, no redistribution,
and no rotating frame. Common and recursive clocks are supported. All three
options remain false by default. This policy changes the pinned spatial
algorithm and is not original-PeleC qualification.

For common clocks, at each spatial evaluation, evaluate EOS pressure from each cell's gas and
held internal energy before restriction. Average this scalar from finest to
coarsest using fluid-volume weights. Interpolate that pressure on the same
coarse stencil and offsets as the material state. Adjust total and held energy
by the same amount to represent the target pressure, retaining density,
momentum, species and the total-minus-held gap. Covered coarse cells and
same-level covered peer donors receive the identical representation. Input
states and accepted valid leaves are never modified by this halo operation.

MPI keeps the scalar owner-local, uses the existing volume-weighted field
restriction and donor exchange, and agrees policy bits and failures
collectively. A failed reconstruction leaves all output halos unchanged.
Hydro and molecular transport consume the same prepared representation.
Saved held intervals and source prefixes bind the pressure selection; restart
files carry `PELEF_EB_TREE_PAIRED_PRESSURE_V1`. A changed policy is refused.

No new valid-leaf energy forcing is introduced. Existing accepted-step
energy reset and numerical-energy accounting remain responsible for changes
to accepted conserved energy. The application budget tests independently
check this boundary. Nonisothermal mixing can still generate pressure error;
this policy is an interface representation, not a solution to arbitrary
nonisothermal material-contact evolution.

Recursive clocks prepare covered parent pressure at synchronized subtree
starts, before capturing history. Provisional spatial stages use current
level EOS pressure; ancestor donors interpolate the two endpoint EOS
pressures at the same time fraction as material data. They do not evaluate
pressure from the conservatively mixed endpoint state or restrict descendants
at a different time. The covered-state reconstruction is private and
transactional. See [recursive evidence](../validation/eb-pressure-recursive-20261022.md).

The completed variable-velocity reacting reference comparison fails the
original field gates. Its four-policy ablation places the dominant departure
at the separate-momentum prerequisite, while adding paired pressure slightly
reduces temperature difference. Constant-velocity pressure convergence does
not qualify this complete policy for variable-velocity reacting reference
agreement. See [the reacting evidence](../validation/eb-pressure-reacting-20261022.md).
