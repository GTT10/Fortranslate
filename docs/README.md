# Documentation index

Read [current status](current_status.md) for the present capability boundary and
qualification provenance, then [the roadmap](completion_roadmap.md) for the
next implementation unit. Use [contributor instructions](../CONTRIBUTING.md)
for repeatable build/test commands.

## Where information belongs

| Question | Record |
| --- | --- |
| What is supported now, and what was actually checked? | [Current status](current_status.md) |
| What must be implemented next? | [Completion roadmap](completion_roadmap.md), issue #101 |
| How is the code organized? | [Architecture](architecture.md), source map below |
| What equations/discretizations are used? | [Numerical methods](numerical_methods.md), [state variables](state_variables.md) |
| How does this map to the reference solver? | [PeleC mapping](pelec_mapping.md), [pinned manifest](../references/pelec_baseline.json) |
| Why was a design chosen? | [Design decisions](design_decisions/) |
| Where is a specific test result? | [Dated validation records](validation/) and their evidence |
| Where are historical implementation details? | [Implementation inventory](implementation_status.md), [porting plan](porting_plan.md), [parity strategy](parity_strategy.md) |

The long inventories contain chronological milestone descriptions; an old
"not yet supported" or "focused-only" statement is not a current global
status. Do not resolve conflicting claims by choosing the largest test count.
Check the source identity and the scope in the dated validation record.

## Source map

`app/` contains executable front ends; `src/driver/` contains shared application
lifecycles; `src/core/`, `src/hydro/`, `src/physics/`, `src/chemistry/`,
`src/transport/`, `src/amr/`, `src/eb/`, and `src/parallel/` contain numerical
and physical modules. `mechanisms/` and `src/generated/` hold mechanism sources
and generated code. `cmake/` configures dependency and mechanism contracts;
`tests/`, `cases/`, and `tools/` contain regression programs, public inputs, and
checkers. Verify the relevant local directory rather than assuming all
responsibilities are implemented for every dimension.

## Integration evidence

[The checkpoint integrity repair](validation/checkpoint-integrity-20260908.md)
explains the two-byte source mismatch and its independent reproduction.
[The 0.245.0 focused record](validation/0.245.0.md) is retained as historical
transport-persistence evidence. [The dated handoff](handoff-pelef-0.245.0-20260908.md)
reports a later local matrix; its old operational instructions and temporary
paths are not current workflow instructions.
