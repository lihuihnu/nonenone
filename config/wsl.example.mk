# Verified WSL profile used for the unified-repository acceptance gate.
# Copy to config/hpc.local.mk and adjust paths if this WSL installation moves.

MPMC_PLATFORM := local

PETSC_DIR  := /home/mpmc/opt/petsc-3.22.2
PETSC_ARCH := arch-linux-c-opt
FMT_DIR    := /usr
METIS_DIR  := /home/mpmc/opt/metis-compat

VTK_INSTALL   :=
X11_LIB_DIR   := /usr/lib/x86_64-linux-gnu
LOCAL_LIB_DIR := /usr/local/lib

MPI_LAUNCHER      := /usr/bin/mpiexec
MPI_LAUNCHER_ARGS :=
MODULES_ENABLED   := 0

PETSC_CUSTOM_HOOK_LINKAGE := both
MPI_MODULE  :=
BLAS_MODULE :=
SBATCH      :=

DQ_MESH_DIR := /home/mpmc/data/DQ_data
