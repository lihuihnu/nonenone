# =============================================================================
# MPMC_SCW case manager
#
#   make -j                       compile all cases
#   make DQcase -j                compile one case
#   make list                     list detected cases
#   make prepare CASE=DQcase      create dirs; refresh managed run.sh from active profile
#   make reset-run CASE=DQcase    replace custom/legacy run.sh with managed template
#   make submit CASE=DQcase       sbatch the existing/generated run.sh
#   make clean CASE=DQcase        remove build/run data, but preserve run.sh
#   make clean-all                clean every case, but preserve all run.sh files
#
# New case convention:
#   <case>/<case>.cpp
#   <case>/case_config.hpp
#   <case>/well_config.hpp
# No central case list needs editing.
# =============================================================================

# PETSc 自带 rules 也定义目标；显式指定默认目标，保证 `make -j` 就是编译全部算例。
.DEFAULT_GOAL := all

ROOT := ..
BUILD_DIR := .build

CASE_CONFIGS := $(wildcard */case_config.hpp)
CASES := $(sort $(patsubst %/case_config.hpp,%,$(CASE_CONFIGS)))

ifeq ($(strip $(CASES)),)
$(error No case found. Each case needs <name>/case_config.hpp)
endif

# -----------------------------------------------------------------------------
# External dependencies. Machine-specific paths live in one project profile.
# config/hpc.local.mk is optional, ignored by Git, and loaded before defaults.
# -----------------------------------------------------------------------------
-include ../config/hpc.local.mk
include ../config/hpc.mk

METIS_INC := -isystem $(METIS_DIR)/include
METIS_LIB := -L$(METIS_DIR)/install/lib -lmetis -Wl,-rpath,$(METIS_DIR)/install/lib

# Default MPI ranks used when run.sh is generated. Override, e.g.
#   make prepare CASE=DQcase TASKS=8
# Case-specific launcher defaults, aligned with the supplied MPMC_LH run scripts.
TASKS_3p6c_original ?= 1
PARTITION_3p6c_original ?= mt_module
SNES_MAX_IT_3p6c_original ?= 100
SNES_ATOL_3p6c_original ?= 1.e-8
SNES_RTOL_3p6c_original ?= 1.e-100
SNES_STOL_3p6c_original ?= 1.e-100
KSP_ATOL_3p6c_original ?= 1.e-12
KSP_RTOL_3p6c_original ?= 1.e-100

TASKS_DQcase ?= 8
PARTITION_DQcase ?= ft_module
SNES_MAX_IT_DQcase ?= 20
SNES_ATOL_DQcase ?= 1.e-4
SNES_RTOL_DQcase ?= 1.e-100
SNES_STOL_DQcase ?=
KSP_ATOL_DQcase ?= 1.e-10
KSP_RTOL_DQcase ?= 1.e-100

# First fully-compositional O/G/W reservoir smoke case.  Keep the default
# launcher small; the purpose is functional validation before scalability.
TASKS_3p4c_pr_reservoir ?= 1
PARTITION_3p4c_pr_reservoir ?= mt_module
SNES_MAX_IT_3p4c_pr_reservoir ?= 30
SNES_ATOL_3p4c_pr_reservoir ?= 1.e-7
SNES_RTOL_3p4c_pr_reservoir ?= 1.e-8
SNES_STOL_3p4c_pr_reservoir ?= 1.e-100
KSP_ATOL_3p4c_pr_reservoir ?= 1.e-10
KSP_RTOL_3p4c_pr_reservoir ?= 1.e-8

# Three-dimensional PR/SW/CPA comparison benchmark: 6912 cells and one
# common 5-component O/G/W problem. Eight ranks are a practical default.
TASKS_three_eos_3d_compare ?= 8
PARTITION_three_eos_3d_compare ?= ft_module
SNES_MAX_IT_three_eos_3d_compare ?= 50
SNES_ATOL_three_eos_3d_compare ?= 1.e-6
SNES_RTOL_three_eos_3d_compare ?= 1.e-8
SNES_STOL_three_eos_3d_compare ?= 1.e-100
KSP_ATOL_three_eos_3d_compare ?= 1.e-10
KSP_RTOL_three_eos_3d_compare ?= 1.e-8

# Compact five-component comparison with independently selected PR/SW/CPA
# parameter packages. Two local ranks are sufficient for the 2000-cell grid.
TASKS_five_component_eos_tuned_compare ?= 2
PARTITION_five_component_eos_tuned_compare ?= mt_module
SNES_MAX_IT_five_component_eos_tuned_compare ?= 50
SNES_ATOL_five_component_eos_tuned_compare ?= 1.e-6
SNES_RTOL_five_component_eos_tuned_compare ?= 1.e-8
SNES_STOL_five_component_eos_tuned_compare ?= 1.e-100
KSP_ATOL_five_component_eos_tuned_compare ?= 1.e-10
KSP_RTOL_five_component_eos_tuned_compare ?= 1.e-8

# Same five-physical-component PR problem using the legacy oil/gas
# compositional model plus an independent immiscible-water equation.
TASKS_five_component_eos_tuned_legacy_compare ?= 2
PARTITION_five_component_eos_tuned_legacy_compare ?= mt_module
SNES_MAX_IT_five_component_eos_tuned_legacy_compare ?= 50
SNES_ATOL_five_component_eos_tuned_legacy_compare ?= 1.e-6
SNES_RTOL_five_component_eos_tuned_legacy_compare ?= 1.e-8
SNES_STOL_five_component_eos_tuned_legacy_compare ?= 1.e-100
KSP_ATOL_five_component_eos_tuned_legacy_compare ?= 1.e-10
KSP_RTOL_five_component_eos_tuned_legacy_compare ?= 1.e-8

# Ma et al. (2021) five-component true O/G/W mutual-solubility benchmark.
# The 40-cell reservoir wrapper is a functional literature-fluid validation,
# not a scalability case, so keep the default launcher serial/conservative.
TASKS_ma2021_5c_three_phase_reservoir ?= 1
PARTITION_ma2021_5c_three_phase_reservoir ?= mt_module
SNES_MAX_IT_ma2021_5c_three_phase_reservoir ?= 40
SNES_ATOL_ma2021_5c_three_phase_reservoir ?= 1.e-7
SNES_RTOL_ma2021_5c_three_phase_reservoir ?= 1.e-8
SNES_STOL_ma2021_5c_three_phase_reservoir ?= 1.e-100
KSP_ATOL_ma2021_5c_three_phase_reservoir ?= 1.e-10
KSP_RTOL_ma2021_5c_three_phase_reservoir ?= 1.e-8


# Panfili et al. (2025) use case #3: 10-component live-water SW case.
# 47x99x15 structured surrogate gives ~70k cells and >2M primary DOFs, so
# default to a modest multi-node launcher rather than serial.
TASKS_panfili2025_case3_fullphysics ?= 16
PARTITION_panfili2025_case3_fullphysics ?= ft_module
SNES_MAX_IT_panfili2025_case3_fullphysics ?= 50
SNES_ATOL_panfili2025_case3_fullphysics ?= 1.e-6
SNES_RTOL_panfili2025_case3_fullphysics ?= 1.e-8
SNES_STOL_panfili2025_case3_fullphysics ?= 1.e-100
KSP_ATOL_panfili2025_case3_fullphysics ?= 1.e-10
KSP_RTOL_panfili2025_case3_fullphysics ?= 1.e-8

# Common 60x20x1 H2O-CO2-nC10 comparison. Traditional and New-PR/SW/CPA
# intentionally use the same launcher and nonlinear/linear tolerances.
TASKS_h2o_co2_nc10_2d_traditional ?= 2
PARTITION_h2o_co2_nc10_2d_traditional ?= mt_module
SNES_MAX_IT_h2o_co2_nc10_2d_traditional ?= 50
SNES_ATOL_h2o_co2_nc10_2d_traditional ?= 1.e-50
SNES_RTOL_h2o_co2_nc10_2d_traditional ?= 1.e-8
SNES_STOL_h2o_co2_nc10_2d_traditional ?=
KSP_ATOL_h2o_co2_nc10_2d_traditional ?= 1.e-50
KSP_RTOL_h2o_co2_nc10_2d_traditional ?= 1.e-5

TASKS_h2o_co2_nc10_2d_benchmark ?= 2
PARTITION_h2o_co2_nc10_2d_benchmark ?= mt_module
SNES_MAX_IT_h2o_co2_nc10_2d_benchmark ?= 50
SNES_ATOL_h2o_co2_nc10_2d_benchmark ?= 1.e-50
SNES_RTOL_h2o_co2_nc10_2d_benchmark ?= 1.e-8
SNES_STOL_h2o_co2_nc10_2d_benchmark ?=
KSP_ATOL_h2o_co2_nc10_2d_benchmark ?= 1.e-50
KSP_RTOL_h2o_co2_nc10_2d_benchmark ?= 1.e-5

# Compact 60-cell strict-SCW kerogen-product screening cases.  Both staged
# fluids use the same conservative serial solver defaults for direct comparison.
TASKS_scw_kerogen_squalane_1d ?= 1
PARTITION_scw_kerogen_squalane_1d ?= mt_module
SNES_MAX_IT_scw_kerogen_squalane_1d ?= 50
SNES_ATOL_scw_kerogen_squalane_1d ?= 2.e-6
SNES_RTOL_scw_kerogen_squalane_1d ?= 1.e-8
SNES_STOL_scw_kerogen_squalane_1d ?= 1.e-100
KSP_ATOL_scw_kerogen_squalane_1d ?= 1.e-12
KSP_RTOL_scw_kerogen_squalane_1d ?= 1.e-8

TASKS_scw_kerogen_lmh_1d ?= 1
PARTITION_scw_kerogen_lmh_1d ?= mt_module
SNES_MAX_IT_scw_kerogen_lmh_1d ?= 50
SNES_ATOL_scw_kerogen_lmh_1d ?= 2.e-6
SNES_RTOL_scw_kerogen_lmh_1d ?= 1.e-8
SNES_STOL_scw_kerogen_lmh_1d ?= 1.e-100
KSP_ATOL_scw_kerogen_lmh_1d ?= 1.e-12
KSP_RTOL_scw_kerogen_lmh_1d ?= 1.e-8

# Homogeneous 60x20x1 extensions of the kerogen-product screening cases.
# Two ranks are a practical local default for the 1200-cell transverse-flow
# grid; tolerances intentionally match the one-dimensional controls.
TASKS_scw_kerogen_squalane_2d ?= 2
PARTITION_scw_kerogen_squalane_2d ?= mt_module
SNES_MAX_IT_scw_kerogen_squalane_2d ?= 50
SNES_ATOL_scw_kerogen_squalane_2d ?= 2.e-6
SNES_RTOL_scw_kerogen_squalane_2d ?= 1.e-8
SNES_STOL_scw_kerogen_squalane_2d ?= 1.e-100
KSP_ATOL_scw_kerogen_squalane_2d ?= 1.e-12
KSP_RTOL_scw_kerogen_squalane_2d ?= 1.e-8

TASKS_scw_kerogen_lmh_2d ?= 2
PARTITION_scw_kerogen_lmh_2d ?= mt_module
SNES_MAX_IT_scw_kerogen_lmh_2d ?= 50
SNES_ATOL_scw_kerogen_lmh_2d ?= 2.e-6
SNES_RTOL_scw_kerogen_lmh_2d ?= 1.e-8
SNES_STOL_scw_kerogen_lmh_2d ?= 1.e-100
KSP_ATOL_scw_kerogen_lmh_2d ?= 1.e-12
KSP_RTOL_scw_kerogen_lmh_2d ?= 1.e-8

# Sun et al. (2024) Exp. 12 derived 48-cell independent-SCW/PR sandpack.
TASKS_sun2024_exp12_scw_co2_nc16_1d ?= 1
PARTITION_sun2024_exp12_scw_co2_nc16_1d ?= mt_module
SNES_MAX_IT_sun2024_exp12_scw_co2_nc16_1d ?= 50
SNES_ATOL_sun2024_exp12_scw_co2_nc16_1d ?= 1.e-50
SNES_RTOL_sun2024_exp12_scw_co2_nc16_1d ?= 1.e-8
SNES_STOL_sun2024_exp12_scw_co2_nc16_1d ?=
KSP_ATOL_sun2024_exp12_scw_co2_nc16_1d ?= 1.e-50
KSP_RTOL_sun2024_exp12_scw_co2_nc16_1d ?= 1.e-5

# Reproducible Sun-2024 factorial.  Lab-scale component inventories are much
# smaller than 1 kg, so the generic 1e-7 absolute SNES tolerance can accept a
# time step without updating the state and accumulate visible mass error.
TASKS_sun2024_scw_co2_nc16_factorial_1d ?= 1
PARTITION_sun2024_scw_co2_nc16_factorial_1d ?= mt_module
SNES_MAX_IT_sun2024_scw_co2_nc16_factorial_1d ?= 50
SNES_ATOL_sun2024_scw_co2_nc16_factorial_1d ?= 1.e-12
SNES_RTOL_sun2024_scw_co2_nc16_factorial_1d ?= 1.e-8
SNES_STOL_sun2024_scw_co2_nc16_factorial_1d ?= 1.e-12
KSP_ATOL_sun2024_scw_co2_nc16_factorial_1d ?= 1.e-12
KSP_RTOL_sun2024_scw_co2_nc16_factorial_1d ?= 1.e-8

# Five-point CO2/SCW displacement matrix on the shared 60x20x1 nC10 grid.
TASKS_h2o_co2_bsb_lumped_3d_lab ?= 1
PARTITION_h2o_co2_bsb_lumped_3d_lab ?= mt_module
SNES_MAX_IT_h2o_co2_bsb_lumped_3d_lab ?= 50
SNES_ATOL_h2o_co2_bsb_lumped_3d_lab ?= 1.e-12
SNES_RTOL_h2o_co2_bsb_lumped_3d_lab ?= 1.e-8
SNES_STOL_h2o_co2_bsb_lumped_3d_lab ?= 1.e-12
KSP_ATOL_h2o_co2_bsb_lumped_3d_lab ?= 1.e-12
KSP_RTOL_h2o_co2_bsb_lumped_3d_lab ?= 1.e-8

# Goals that only manage files/scripts do not require a local PETSc installation.
NO_PETSC_GOALS := list help prepare prepare-% reset-run submit print print-config clean clean-all
BUILD_GOALS := $(filter-out $(NO_PETSC_GOALS),$(MAKECMDGOALS))
ifeq ($(strip $(MAKECMDGOALS)),)
BUILD_GOALS := all
endif

ifneq ($(strip $(BUILD_GOALS)),)
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
endif

# -----------------------------------------------------------------------------
# Project sources. Natural/model/ad/output/adaptive are header-only here;
# common and grid contain the compiled implementation units.
# -----------------------------------------------------------------------------
COMMON_SRCS := $(wildcard $(ROOT)/common/src/*.cpp)
GRID_SRCS   := $(wildcard $(ROOT)/grid/src/*.cpp) $(wildcard $(ROOT)/grid/src/*.C)
GRID_PRIVATE_HEADERS := $(wildcard $(ROOT)/grid/src/*.hpp)
COMMON_OBJS := $(patsubst $(ROOT)/common/src/%.cpp,$(BUILD_DIR)/common/%.o,$(COMMON_SRCS))
GRID_OBJS   := $(patsubst $(ROOT)/grid/src/%.cpp,$(BUILD_DIR)/grid/%.o,$(filter %.cpp,$(GRID_SRCS))) \
               $(patsubst $(ROOT)/grid/src/%.C,$(BUILD_DIR)/grid/%.o,$(filter %.C,$(GRID_SRCS)))
PROJECT_OBJS := $(COMMON_OBJS) $(GRID_OBJS)

# The target cluster PETSc is C++-linked, while a local PETSc build may expose
# bare C undefined symbols for the historical updateState/updateSol hooks.
# `both` adds a thin C bridge without removing the existing C++ symbols.
VALID_PETSC_CUSTOM_HOOK_LINKAGES := cxx both
ifeq ($(filter $(PETSC_CUSTOM_HOOK_LINKAGE),$(VALID_PETSC_CUSTOM_HOOK_LINKAGES)),)
$(error PETSC_CUSTOM_HOOK_LINKAGE must be one of: $(VALID_PETSC_CUSTOM_HOOK_LINKAGES))
endif
PETSC_HOOK_BRIDGE_OBJS :=
ifeq ($(PETSC_CUSTOM_HOOK_LINKAGE),both)
PETSC_HOOK_BRIDGE_OBJS += $(BUILD_DIR)/support/petsc_custom_hooks_c_bridge.o
endif
CASE_OBJS := $(addprefix $(BUILD_DIR)/cases/,$(addsuffix .o,$(CASES)))

# Header dependency tracking intentionally uses Make prerequisites rather than
# compiler-generated .d files.  PETSc installations differ in how
# PETSC_CXXCOMPILE_SINGLE forwards GCC dependency flags; keeping those flags out
# of the PETSc compiler wrapper is the most portable choice.
#
# A public-header change therefore rebuilds the compiled project/case objects.
# The project is small enough that this coarse-grained dependency is preferable
# to fragile compiler-specific dependency-file plumbing.
HEADER_DIRS := \
    $(ROOT)/ad/include \
    $(ROOT)/common/include \
    $(ROOT)/grid/include \
    $(ROOT)/indices/include \
    $(ROOT)/models/include \
    $(ROOT)/output/include \
    $(ROOT)/AdaptiveTimeStepper/include \
    include
PUBLIC_HEADERS := $(shell find $(HEADER_DIRS) -type f -name '*.hpp' 2>/dev/null) \
                  $(shell find $(HEADER_DIRS) -type f -name '*.h' 2>/dev/null)

INCLUDES := \
    -I$(ROOT)/ad/include \
    -I$(ROOT)/common/include \
    -I$(ROOT)/grid/include \
    -I$(ROOT)/indices/include \
    -I$(ROOT)/models/include \
    -I$(ROOT)/output/include \
    -I$(ROOT)/AdaptiveTimeStepper/include \
    -I$(CURDIR)/include \
    -I$(FMT_DIR)/include \
    $(METIS_INC)

CFLAGS   += -fvisibility=default -fPIC
CXXFLAGS += -fvisibility=default -fPIC -DFMT_HEADER_ONLY -std=gnu++17
CXXFLAGS += -O2 -Wall -Wextra -Wpedantic -Werror $(INCLUDES)

.PHONY: all list help prepare reset-run submit print print-config clean clean-all $(CASES)

all: $(CASES)
	@echo "============================================================"
	@echo "[DONE] all cases compiled: $(CASES)"
	@echo "============================================================"

list:
	@echo "Detected cases:"
	@$(foreach c,$(CASES),echo "  - $(c)";)

help:
	@echo "make -j                         compile all cases"
	@echo "make <case> -j                  compile one case"
	@echo "make prepare CASE=<case>        create dirs; refresh managed run.sh from current profile"
	@echo "make reset-run CASE=<case>      explicitly regenerate run.sh from template"
	@echo "make submit CASE=<case>         submit with sbatch"
	@echo "make clean CASE=<case>          clean build/run data; preserve run.sh"
	@echo "make clean-all                  clean every case + .build; preserve run.sh"
	@echo "make print-config               print resolved local/HPC machine profile"
	@echo "make list                       list auto-detected cases"

# -----------------------------------------------------------------------------
# Script/directory generation
# -----------------------------------------------------------------------------

print-config:
	@echo "MPMC_PLATFORM          = $(MPMC_PLATFORM)"
	@echo "PETSC_DIR         = $(PETSC_DIR)"
	@echo "PETSC_ARCH        = $(PETSC_ARCH)"
	@echo "FMT_DIR           = $(FMT_DIR)"
	@echo "METIS_DIR         = $(METIS_DIR)"
	@echo "VTK_INSTALL       = $(VTK_INSTALL)"
	@echo "MPI_LAUNCHER      = $(MPI_LAUNCHER)"
	@echo "MPI_LAUNCHER_ARGS = $(MPI_LAUNCHER_ARGS)"
	@echo "MODULES_ENABLED   = $(MODULES_ENABLED)"
	@echo "PETSC_CUSTOM_HOOK_LINKAGE = $(PETSC_CUSTOM_HOOK_LINKAGE)"
	@echo "DQ_MESH_DIR       = $(DQ_MESH_DIR)"

prepare:
	@if [ -z "$(CASE)" ]; then echo "Usage: make prepare CASE=<case>"; exit 2; fi
	@$(MAKE) --no-print-directory prepare-$(CASE) TASKS="$(TASKS)"

prepare-%:
	@if [ ! -f "$*/case_config.hpp" ]; then echo "Unknown case: $*"; exit 2; fi
	@mkdir -p "$*/bin" "$*/logs" "$*/results"
	@managed=0; \
	  if [ -f "$*/run.sh" ] && grep -q '^# MPMC_SCW_GENERATED_RUN_SCRIPT=1$$' "$*/run.sh"; then managed=1; fi; \
	  if [ -f "$*/run.sh" ] && [ "$$managed" -eq 0 ]; then \
	    echo "[KEEP]  $*/run.sh (custom/legacy script; use 'make reset-run CASE=$*' to adopt the managed template)"; \
	  else \
	    tasks="$(TASKS)"; \
	    if [ -z "$$tasks" ]; then tasks="$($(strip TASKS_$*))"; fi; \
	    if [ -z "$$tasks" ]; then tasks=1; fi; \
	    partition="$($(strip PARTITION_$*))"; \
	    snes_max_it="$($(strip SNES_MAX_IT_$*))"; \
	    snes_atol="$($(strip SNES_ATOL_$*))"; \
	    snes_rtol="$($(strip SNES_RTOL_$*))"; \
	    snes_stol="$($(strip SNES_STOL_$*))"; \
	    ksp_atol="$($(strip KSP_ATOL_$*))"; \
	    ksp_rtol="$($(strip KSP_RTOL_$*))"; \
	    : "$${snes_max_it:=50}"; \
	    : "$${snes_atol:=1.e-7}"; \
	    : "$${snes_rtol:=1.e-8}"; \
	    : "$${ksp_atol:=1.e-10}"; \
	    : "$${ksp_rtol:=1.e-8}"; \
	    default_mesh_dir=""; \
	    if [ "$*" = "DQcase" ]; then default_mesh_dir="$(DQ_MESH_DIR)"; fi; \
	    vtk_lib=""; if [ -n "$(VTK_INSTALL)" ]; then vtk_lib="$(VTK_INSTALL)/lib"; fi; \
	    metis_lib=""; if [ -n "$(METIS_DIR)" ]; then metis_lib="$(METIS_DIR)/install/lib"; fi; \
	    sed \
	      -e 's|@CASE@|$*|g' \
	      -e 's|@TASKS@|'"$$tasks"'|g' \
	      -e 's|@PARTITION@|'"$$partition"'|g' \
	      -e 's|@SNES_MAX_IT@|'"$$snes_max_it"'|g' \
	      -e 's|@SNES_ATOL@|'"$$snes_atol"'|g' \
	      -e 's|@SNES_RTOL@|'"$$snes_rtol"'|g' \
	      -e 's|@SNES_STOL@|'"$$snes_stol"'|g' \
	      -e 's|@KSP_ATOL@|'"$$ksp_atol"'|g' \
	      -e 's|@KSP_RTOL@|'"$$ksp_rtol"'|g' \
	      -e 's|@CASE_DIR@|$(abspath $*)|g' \
	      -e 's|@CASE_ROOT@|$(abspath .)|g' \
	      -e 's|@MPMC_PLATFORM@|$(MPMC_PLATFORM)|g' \
	      -e 's|@PETSC_DIR@|$(PETSC_DIR)|g' \
	      -e 's|@PETSC_ARCH@|$(PETSC_ARCH)|g' \
	      -e 's|@VTK_LIB@|'"$$vtk_lib"'|g' \
	      -e 's|@METIS_LIB@|'"$$metis_lib"'|g' \
	      -e 's|@LOCAL_LIB@|$(LOCAL_LIB_DIR)|g' \
	      -e 's|@X11_LIB@|$(X11_LIB_DIR)|g' \
	      -e 's|@MODULES_ENABLED@|$(MODULES_ENABLED)|g' \
	      -e 's|@MPI_MODULE@|$(MPI_MODULE)|g' \
	      -e 's|@BLAS_MODULE@|$(BLAS_MODULE)|g' \
	      -e 's|@MPI_LAUNCHER@|$(MPI_LAUNCHER)|g' \
	      -e 's|@MPI_LAUNCHER_ARGS@|$(MPI_LAUNCHER_ARGS)|g' \
	      -e 's|@DEFAULT_MESH_DIR@|'"$$default_mesh_dir"'|g' \
	      run.sh.in > "$*/run.sh"; \
	    chmod +x "$*/run.sh"; \
	    echo "[CREATE] $*/run.sh"; \
	  fi
	@echo "[READY] $*: bin/ logs/ results/ run.sh"

# reset-run deliberately replaces custom/legacy run.sh with the managed template.
reset-run:
	@if [ -z "$(CASE)" ]; then echo "Usage: make reset-run CASE=<case> [TASKS=N]"; exit 2; fi
	@if [ ! -f "$(CASE)/case_config.hpp" ]; then echo "Unknown case: $(CASE)"; exit 2; fi
	@rm -f "$(CASE)/run.sh"
	@$(MAKE) --no-print-directory prepare-$(CASE) TASKS="$(TASKS)"
	@echo "[RESET] $(CASE)/run.sh regenerated from run.sh.in"

submit:
	@if [ -z "$(CASE)" ]; then echo "Usage: make submit CASE=<case>"; exit 2; fi
	@if [ "$(MPMC_PLATFORM)" = "local" ]; then echo "submit is disabled for MPMC_PLATFORM=local; use 'make run CASE=$(CASE) NP=N'"; exit 2; fi
	@$(MAKE) --no-print-directory $(CASE) TASKS="$(TASKS)"
	@$(MAKE) --no-print-directory prepare-$(CASE) TASKS="$(TASKS)"
	@if [ -z "$(SBATCH)" ]; then echo "SBATCH is not configured"; exit 2; fi
	@echo "[SBATCH] $(CASE)/run.sh"
	@cd "$(CASE)" && $(SBATCH) run.sh

print:
	@if [ -z "$(CASE)" ]; then echo "Usage: make print CASE=<case>"; exit 2; fi
	@echo "case       = $(CASE)"
	@echo "source     = $(CASE)/$(CASE).cpp"
	@echo "binary     = $(CASE)/bin/$(CASE)"
	@echo "run script = $(CASE)/run.sh"

# -----------------------------------------------------------------------------
# Compilation
# -----------------------------------------------------------------------------
$(BUILD_DIR)/common/%.o: $(ROOT)/common/src/%.cpp $(PUBLIC_HEADERS)
	@mkdir -p $(dir $@)
	@echo "[CXX][common] $<"
	@$(PETSC_CXXCOMPILE_SINGLE) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/grid/%.o: $(ROOT)/grid/src/%.cpp $(PUBLIC_HEADERS) $(GRID_PRIVATE_HEADERS)
	@mkdir -p $(dir $@)
	@echo "[CXX][grid] $<"
	@$(PETSC_CXXCOMPILE_SINGLE) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/grid/%.o: $(ROOT)/grid/src/%.C $(PUBLIC_HEADERS) $(GRID_PRIVATE_HEADERS)
	@mkdir -p $(dir $@)
	@echo "[CXX][grid] $<"
	@$(PETSC_CXXCOMPILE_SINGLE) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/support/petsc_custom_hooks_c_bridge.o: src/petsc_custom_hooks_c_bridge.cpp $(PUBLIC_HEADERS)
	@mkdir -p $(dir $@)
	@echo "[CXX][support] $<"
	@$(PETSC_CXXCOMPILE_SINGLE) $(CXXFLAGS) -c $< -o $@

# One compile/link rule is generated for every auto-detected case.
define CASE_RULE
$(BUILD_DIR)/cases/$(1).o: $(1)/$(1).cpp $$(wildcard $(1)/*.hpp) $(PUBLIC_HEADERS)
	@mkdir -p $$(dir $$@)
	@echo "[CXX][case] $$<"
	@$$(PETSC_CXXCOMPILE_SINGLE) $$(CXXFLAGS) -I$(1) -c $$< -o $$@

$(1): $(1)/bin/$(1)

$(1)/bin/$(1): $$(PROJECT_OBJS) $$(PETSC_HOOK_BRIDGE_OBJS) $(BUILD_DIR)/cases/$(1).o | prepare-$(1)
	@echo "[LD][case] $$@"
	@$$(CXXLINKER) -o $$@ $$(PROJECT_OBJS) $$(PETSC_HOOK_BRIDGE_OBJS) $(BUILD_DIR)/cases/$(1).o \
	    $$(METIS_LIB) $$(PETSC_LIB) $$(LDFLAGS)
	@echo "[DONE] $(1) -> $$@"
endef
$(foreach c,$(CASES),$(eval $(call CASE_RULE,$(c))))

# PETSc's rules may also define clean::, therefore use double-colon deliberately.
clean::
	@if [ -n "$(CASE)" ]; then \
	    echo "[CLEAN] $(CASE)"; \
	    rm -rf "$(CASE)/bin" "$(CASE)/logs" "$(CASE)/results"; \
	    rm -f "$(BUILD_DIR)/cases/$(CASE).o"; \
	  else \
	    echo "Usage: make clean CASE=<case>  (or make clean-all)"; \
	  fi

clean-all::
	@echo "[CLEAN] all generated case files"
	@$(foreach c,$(CASES),rm -rf "$(c)/bin" "$(c)/logs" "$(c)/results";)
	@rm -rf "$(BUILD_DIR)"
