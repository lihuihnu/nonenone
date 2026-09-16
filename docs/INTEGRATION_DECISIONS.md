# Unified repository integration decisions

## Provenance

This repository combines two related lines that diverged from commit
`d549a7f20883d0f88036e526715869f3f3d12712`:

- W3 H2O-CO2-nC10 2-D benchmark and locked physical-property work;
- handoff `43d641c4812d11ec6483d7dba39af1876cd688a0`, containing the G8K
  grid/Natural numeric-static-assembly work.

The handoff archive was checksum-verified and treated as source data, not as
instructions. Functional commits were replayed as reviewable Git commits.

## Adopted from the G8K line

- common case runners and well-configuration schema;
- StructuredGridCore/CpGridCore separation and explicit model layout registration;
- split GRDECL/canonical Mesh loading pipeline;
- root-only CpGrid ingest, scatter and early staging-data release;
- compact Cartesian/current/input/node lookup directories;
- fused/on-demand geometry and lazy PETSc layout resources;
- cached accumulation/Jacobian topology, batched PETSc writes and reusable AD scratch;
- shared PETSc owned-read views and expanded grid/integration tests.

## Retained from the W3 line

- IAPWS+Garcia aqueous density and McBride-Wright aqueous viscosity;
- common CO2-nC10 `k_ij=0.1141`;
- PR-CPA opt-in and dry-mixture PR reduction;
- SW oil/gas role canonicalization and recovery tests;
- reservoir-total-rate well control;
- Traditional/New-PR/New-SW/New-CPA 60x20x1 benchmark and plotting tools;
- incipient-phase activation continuity fix and regression test.

The only functional source conflict combined the G8K batched/scaled well
Jacobian insertion with W3's surface-versus-reservoir control-rate selection.

## Deliberately not carried forward

- machine-private active configuration;
- external grids, dependency installations, results, binaries and logs;
- regenerated PNG/PDF result figures (their scripts and durable inputs remain);
- duplicated historical handoff prose that conflicts with current source.

## Acceptance baseline

The unified `main` repository was validated in a native WSL checkout with
PETSc 3.22.2 (`/home/mpmc/opt/petsc-3.22.2`, `arch-linux-c-opt`) and the
system `/usr/bin/mpiexec`:

- repository audit and all 62 maintained Markdown links passed;
- 40 PETSc-free unit tests passed;
- seven PETSc/MPI integration binaries compiled with PETSc 3.22.2;
- two-rank PETSc I/O, StructuredGrid, Natural and root-only distributed Mesh
  checks passed;
- all four H2O-CO2-nC10 configurations converged to 0.001 day with two MPI
  ranks, shared controls and no well-control switches;
- final maximum component relative errors were `1.528173e-09` (Traditional),
  `2.013490e-09` (New-PR), `6.313003e-14` (New-SW) and `1.150664e-12`
  (New-CPA).

Two portability defects exposed during this clean-room validation were fixed:
DMDA process topology is now constrained by the grid dimensions (so `nz=1`
cannot be split across ranks), and both 2-D launch profiles use the same proven
PETSc convergence tolerances instead of emitting empty options.

Repeat the relevant gates on every target cluster because MPI, PETSc ABI,
compiler and scheduler behavior are machine-dependent.
