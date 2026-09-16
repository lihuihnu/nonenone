# Machine configuration

All source code is platform-neutral. Machine paths, launchers and modules are
selected through one ignored override file:

```bash
cp config/wsl.example.mk config/hpc.local.mk      # current WSL
cp config/local.example.mk config/hpc.local.mk    # another Linux workstation
# or use config/hpc.mk as the target-cluster baseline
```

Edit `config/hpc.local.mk`, then run:

```bash
make print-config
make doctor-local       # Linux / WSL
make doctor-hpc         # cluster
```

`config/hpc.local.mk` is intentionally ignored by Git. It is loaded before
the tracked `config/hpc.mk`, whose `?=` assignments only fill missing values.
Use `:=` in the private file for deterministic overrides.

Rules:

- normal `mpiexec`/`mpirun` profiles must set `MPI_LAUNCHER_ARGS :=`;
- `--mpi=pmi2` belongs to the target `yhrun` profile only;
- external PETSc, MPI, METIS, fmt, VTK and grid data are never copied into the
  repository;
- run `make doctor` after moving the repository to another machine.
