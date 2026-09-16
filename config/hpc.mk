# MPMC_SCW internal HPC profile
#
# This repository is maintained for the user's target supercomputer. Keep all
# machine-specific dependency paths here instead of duplicating them in case/
# and test/ Makefiles. For temporary/private changes, copy this file to
# config/hpc.local.mk; that file is ignored by Git and loaded BEFORE this
# profile so local/private ?= assignments can override these defaults.

MPMC_PLATFORM       ?= hpc

PETSC_DIR      ?= /vol8/home/hnu_yhj/jianglei/MPMC_Test/thirdparpties/petsc-3.22.2
PETSC_ARCH     ?= arch-linux-cxx-opt
FMT_DIR        ?= /vol8/home/hnu_yhj/jianglei/MPMC_Test/thirdparpties/fmt
METIS_DIR      ?= /vol8/home/hnu_yhj/jianglei/thirdparpties/metis-5.1.0
VTK_INSTALL    ?= /vol8/home/hnu_yhj/ruili/software/vtk_install
X11_LIB_DIR    ?= /usr/lib/aarch64-linux-gnu
LOCAL_LIB_DIR  ?= /vol8/home/hnu_yhj/yanzhang/local/lib

MPI_MODULE     ?= mpich/mpi-x-pmi2
BLAS_MODULE    ?= LapackCBLAS/3.5.0
MPI_LAUNCHER      ?= yhrun
MPI_TYPE          ?= pmi2
MPI_LAUNCHER_ARGS ?= --mpi=$(MPI_TYPE)
MODULES_ENABLED   ?= 1
# Custom PETSc hook ABI: target cluster PETSc is built as C++ and references
# mangled updateState/updateSol symbols supplied by each case.
PETSC_CUSTOM_HOOK_LINKAGE ?= cxx
SBATCH            ?= sbatch

# Default CpGrid dataset used by DQcase integration tests.
DQ_MESH_DIR    ?= /vol8/home/hnu_yhj/lihui/petsc-3.22.2/MPMC_LH/data/DQ_data/
