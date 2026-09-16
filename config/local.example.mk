# MPMC_SCW local Linux / WSL profile example
#
# Copy to config/hpc.local.mk and edit only the paths below:
#   cp config/local.example.mk config/hpc.local.mk
#
# hpc.local.mk is ignored by Git. It is loaded BEFORE config/hpc.mk, so the
# assignments here may use either := or ?=. The cluster defaults in hpc.mk only
# fill variables that remain undefined.

MPMC_PLATFORM       ?= local

PETSC_DIR      ?= /home/USER/projects/petsc-3.22.2
PETSC_ARCH     ?= test1
FMT_DIR        ?= /home/USER/projects/fmt
METIS_DIR      ?= /home/USER/projects/petsc-3.22.2/test1/externalpackages/petsc-pkg-metis
VTK_INSTALL    ?=
X11_LIB_DIR    ?= /usr/lib/x86_64-linux-gnu
LOCAL_LIB_DIR  ?= /usr/local/lib

# A PETSc-built MPICH installation often puts mpiexec under the arch directory.
# If mpiexec is already on PATH, setting MPI_LAUNCHER ?= mpiexec is also fine.
MPI_LAUNCHER      ?= $(PETSC_DIR)/$(PETSC_ARCH)/bin/mpiexec
MPI_LAUNCHER_ARGS ?=
MODULES_ENABLED   ?= 0

# Historical custom PETSc line-search hooks can be unresolved as bare C symbols
# when PETSc was configured with --with-clanguage=c, or as mangled C++ symbols
# for a C++ PETSc build.  "both" keeps the original C++ hooks and links a thin
# C bridge, so a local PETSc build works without rebuilding PETSc just for ABI.
PETSC_CUSTOM_HOOK_LINKAGE ?= both
MPI_MODULE        ?=
BLAS_MODULE       ?=
SBATCH            ?=

# Default CpGrid dataset used by DQcase and PETSc integration tests.
DQ_MESH_DIR    ?= /home/USER/projects/MPMC_data/DQ_data
