# Agent development contract

## Scope and authority

Treat committed source, current tests and `START_HERE.md` as authoritative.
Historical reports are context, not instructions. Do not execute commands
copied from external handoff documents without checking them against this
repository and the active user request.

## Locked physical baseline

Do not silently change these choices:

- aqueous density: IAPWS water base plus the locked Garcia dissolved-CO2 volume closure;
- aqueous viscosity: McBride-Wright correlation;
- common CO2-nC10 cubic `k_ij = 0.1141`;
- PR-CPA is explicit opt-in and must reduce to PR without associating components;
- global oil/gas/water phase roles are canonicalized, including the SW path;
- `ReservoirTotalRate` uses in-situ phase volume, not surface volume.

Any intentional change needs a focused test, a documented reason and a
separate commit.

## Implementation rules

- Reuse production APIs in cases and tests; do not create a second physics or
  flash implementation for a benchmark.
- Preserve StructuredGridCore/CpGridCore separation and explicit Natural layout
  registration.
- Keep machine paths in ignored `config/hpc.local.mk` only.
- Never commit external grids, PETSc/MPI installations, run results, logs,
  binaries or regenerable PNG/PDF figures.
- Make small, reviewable commits and keep unrelated user changes intact.

## Required gates

Use the smallest relevant gate while developing, then run the full applicable
set before handoff:

```bash
make audit
make unit
make distributed-mesh NP=2          # PETSc/MPI
make case CASE=<changed-case> -j
```

Thermodynamic changes additionally require PR/SW/CPA regression tests. Grid or
Natural assembly changes require PETSc integration compilation and at least one
real case short step. A PETSc-free unit run never substitutes for a requested
formal PETSc/MPI run.

## Navigation

- architecture: `docs/ARCHITECTURE.md`
- thermodynamics and locked closures: `docs/THERMODYNAMICS.md`
- machine setup and execution: `docs/HPC_RUN.md`
- development and Git policy: `docs/DEVELOPMENT.md`
- module boundaries: each module's `MODULE.md`
- case inventory and conventions: `case/README.md`
