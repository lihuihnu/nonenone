# ============================================================================
# MPMC_SCW unified test module
#
# Test levels:
#   unit        : pure C++ tests; no PETSc/MPI runtime is required
#   config      : compile-only checks for template/configuration combinations
#   integration : PETSc/MPI/METIS runtime tests
#   vtk         : VTK output integration test
#
# Common commands:
#   make unit
#   make run-unit
#   make full
#   make run-full NP=1 MESH_DIR=/path/to/DQ_data
# ============================================================================

.DEFAULT_GOAL := unit

# Recipes that use shell pipelines are executed by bash.
SHELL := /bin/bash

PROJECT_ROOT := ..

# ----------------------------------------------------------------------------
# External dependencies
# ----------------------------------------------------------------------------
-include $(PROJECT_ROOT)/config/hpc.local.mk
include $(PROJECT_ROOT)/config/hpc.mk

# PETSc make rules are only needed for the cluster-side targets. This keeps
# `make unit` usable on a normal development machine without PETSc installed.
PETSC_GOALS := full config integration vtk run-integration run-vtk run-full run-distributed-mesh
ifneq ($(strip $(filter $(PETSC_GOALS),$(MAKECMDGOALS))),)
include $(PETSC_DIR)/lib/petsc/conf/variables
include $(PETSC_DIR)/lib/petsc/conf/rules
endif

# ----------------------------------------------------------------------------
# MPMC include paths
# ----------------------------------------------------------------------------
MPMC_INCLUDES := \
	-I$(PROJECT_ROOT)/ad/include \
	-I$(PROJECT_ROOT)/common/include \
	-I$(PROJECT_ROOT)/grid/include \
	-I$(PROJECT_ROOT)/indices/include \
	-I$(PROJECT_ROOT)/models/include \
	-I$(PROJECT_ROOT)/output/include \
	-I$(PROJECT_ROOT)/AdaptiveTimeStepper/include \
	-I$(PROJECT_ROOT)/case/include \
	-I$(PROJECT_ROOT)/test/include \
	-I$(PROJECT_ROOT)/tools/include

COMMON_CXXFLAGS := \
	-std=c++17 \
	-O2 \
	-Wall \
	-Wextra \
	-Wpedantic \
	-Werror \
	-fvisibility=default \
	-fPIC \
	-DFMT_HEADER_ONLY \
	-I$(FMT_DIR)/include \
	$(MPMC_INCLUDES)

# ----------------------------------------------------------------------------
# Test directories
# ----------------------------------------------------------------------------
UNIT_DIR        := src/unit
CONFIG_DIR      := src/compile
INTEGRATION_DIR := src/integration
SUPPORT_DIR     := src/support
BENCHMARK_DIR   := src/benchmark
BUILD_DIR       := build
BIN_DIR         := bin
LOG_DIR         := logs
RESULT_DIR      := results

# ----------------------------------------------------------------------------
# 1. Unit tests: pure C++ module checks
# ----------------------------------------------------------------------------
UNIT_TEST_NAMES := \
	ad_test \
	ad_dimension_test \
	ad_common_integration_test \
	common_test \
	grdecl_test \
	mesh_canonical_data_test \
	cartesian_cell_directory_test \
	cpgrid_geometry_kernel_test \
	structured_cartesian_geometry_test \
	grdecl_rock_payload_test \
	rock_property_validation_test \
	indices_test \
	natural_core_test \
	natural_extended_test \
	three_phase_flash_test \
	phase_role_canonicalization_test \
	sw_phase_ordering_test \
	cpa_single_nonaqueous_phase_test \
	sw_flash_recovery_test \
	sw_heringer_temperature_transition_test \
	thermodynamic_failure_modes_test \
	lab_pseudo3d_geometry_test \
	lab_flow_control_contract_test \
	three_phase_ad_test \
	three_phase_reservoir_case_test \
	three_eos_3d_case_test \
	scw_co2_nc4_nc16_3d_case_test \
	scw_kerogen_squalane_case_test \
	scw_kerogen_lmh_case_test \
	scw_kerogen_squalane_2d_case_test \
	scw_kerogen_lmh_2d_case_test \
	sun2024_factorial_case_test \
	bsb_displacement_design_test \
	five_component_eos_tuned_case_test \
	legacy_3p6c_case_config_test \
	legacy_dq_case_config_test \
	component_mass_balance_test \
	ma2021_case_test \
	panfili2025_case3_test \
	natural_conservation_test \
	well_core_test \
	well_control_test \
	well_natural_test \
	adaptive_timestep_test \
	adaptive_well_test \
	output_test \
	output_feature_test \
	output_well_test \
	tools_phase_diagram_test \
	thermo_empirical_parameters_test \
	eos_parameter_regression_test

UNIT_BINS := $(addprefix $(BIN_DIR)/,$(UNIT_TEST_NAMES))

# Unit tests are header-heavy.  Keep the dependency deliberately coarse so a
# Natural/AD/indices/case-support header change can never reuse a stale test
# executable and report a false PASS.
UNIT_PUBLIC_HEADERS := $(shell find \
	$(PROJECT_ROOT)/ad/include \
	$(PROJECT_ROOT)/common/include \
	$(PROJECT_ROOT)/grid/include \
	$(PROJECT_ROOT)/indices/include \
	$(PROJECT_ROOT)/models/include \
	$(PROJECT_ROOT)/output/include \
	$(PROJECT_ROOT)/AdaptiveTimeStepper/include \
	$(PROJECT_ROOT)/case \
	$(PROJECT_ROOT)/tools/include \
	-type f \( -name '*.hpp' -o -name '*.h' -o -name '*.inc' \) 2>/dev/null)

CXX ?= g++

.PHONY: unit run-unit benchmark-performance

PERFORMANCE_BENCHMARKS := \
	$(BIN_DIR)/cpgrid_geometry_benchmark \
	$(BIN_DIR)/thermo_secondary_benchmark \
	$(BIN_DIR)/three_phase_flash_benchmark \
	$(BIN_DIR)/cell_property_benchmark \
	$(BIN_DIR)/cpa_phase_benchmark

benchmark-performance: $(PERFORMANCE_BENCHMARKS)
	@set -e; for benchmark in $(PERFORMANCE_BENCHMARKS); do \
		echo "[RUN][benchmark] $$benchmark"; \
		$$benchmark; \
	done

$(BIN_DIR)/cpgrid_geometry_benchmark: $(BENCHMARK_DIR)/cpgrid_geometry_benchmark.cpp $(CPGRID_GEOMETRY_UNIT_SRCS) $(UNIT_PUBLIC_HEADERS) | dirs
	@echo "[CXX][benchmark] $< + fused CpGrid geometry kernel"
	@$(CXX) -I$(CPGRID_GEOMETRY_STUB_DIR) $(COMMON_CXXFLAGS) $< $(CPGRID_GEOMETRY_UNIT_SRCS) -o $@

$(BIN_DIR)/%_benchmark: $(BENCHMARK_DIR)/%_benchmark.cpp $(UNIT_PUBLIC_HEADERS) | dirs
	@echo "[CXX][benchmark] $<"
	@$(CXX) $(COMMON_CXXFLAGS) $< -o $@

unit: $(UNIT_BINS)
	@echo "[PASS] unit tests compiled: $(words $(UNIT_TEST_NAMES))"

$(BIN_DIR)/%: $(UNIT_DIR)/%.cpp $(UNIT_PUBLIC_HEADERS) | dirs
	@echo "[CXX][unit] $<"
	@$(CXX) $(COMMON_CXXFLAGS) $< -o $@

GRDECL_UNIT_SRCS := \
	$(PROJECT_ROOT)/grid/src/grdecl.cpp \
	$(PROJECT_ROOT)/grid/src/grdecl_tokenizer.cpp \
	$(PROJECT_ROOT)/grid/src/grdecl_parser.cpp \
	$(PROJECT_ROOT)/grid/src/grdecl_geometry.cpp \
	$(PROJECT_ROOT)/grid/src/grdecl_active.cpp \
	$(PROJECT_ROOT)/grid/src/grdecl_builder.cpp
GRDECL_UNIT_DEPS := $(GRDECL_UNIT_SRCS) $(PROJECT_ROOT)/grid/src/grdecl_detail.hpp

$(BIN_DIR)/grdecl_test: $(UNIT_DIR)/grdecl_test.cpp $(GRDECL_UNIT_DEPS) $(UNIT_PUBLIC_HEADERS) | dirs
	@echo "[CXX][unit] $< + GRDECL implementation"
	@$(CXX) $(COMMON_CXXFLAGS) $< $(GRDECL_UNIT_SRCS) -o $@

MESH_CANONICAL_UNIT_SRCS := \
	$(PROJECT_ROOT)/grid/src/mesh_mrst_input.cpp \
	$(PROJECT_ROOT)/grid/src/mesh_grdecl_input.cpp
MESH_CANONICAL_UNIT_DEPS := $(MESH_CANONICAL_UNIT_SRCS) $(PROJECT_ROOT)/grid/src/mesh_canonical_data.hpp

$(BIN_DIR)/mesh_canonical_data_test: $(UNIT_DIR)/mesh_canonical_data_test.cpp $(MESH_CANONICAL_UNIT_DEPS) $(UNIT_PUBLIC_HEADERS) | dirs
	@echo "[CXX][unit] $< + canonical mesh input adapters"
	@$(CXX) $(COMMON_CXXFLAGS) -I$(PROJECT_ROOT)/grid/src $< $(MESH_CANONICAL_UNIT_SRCS) -o $@

CPGRID_GEOMETRY_UNIT_SRCS := \
	$(PROJECT_ROOT)/grid/src/point.cpp \
	$(PROJECT_ROOT)/grid/src/geometry.cpp \
	$(PROJECT_ROOT)/grid/src/face.cpp \
	$(PROJECT_ROOT)/grid/src/polyhedron.cpp
CPGRID_GEOMETRY_STUB_DIR := $(SUPPORT_DIR)/cpgrid_geometry_stub

$(BIN_DIR)/cpgrid_geometry_kernel_test: $(UNIT_DIR)/cpgrid_geometry_kernel_test.cpp $(CPGRID_GEOMETRY_UNIT_SRCS) $(UNIT_PUBLIC_HEADERS) | dirs
	@echo "[CXX][unit] $< + fused CpGrid geometry kernel"
	@$(CXX) -I$(CPGRID_GEOMETRY_STUB_DIR) $(COMMON_CXXFLAGS) $< $(CPGRID_GEOMETRY_UNIT_SRCS) -o $@

run-unit: unit
	@set -e; \
	for test_bin in $(UNIT_BINS); do \
		echo "------------------------------------------------------------"; \
		echo "[RUN][unit] $$test_bin"; \
		$$test_bin; \
	done
	@echo "============================================================"
	@echo "[PASS] all unit tests"
	@echo "============================================================"

# ----------------------------------------------------------------------------
# 2. Compile-only configuration matrix
# ----------------------------------------------------------------------------
# These tests instantiate important model combinations without running a case.
CONFIG_SOURCES := $(wildcard $(CONFIG_DIR)/*.cpp)
CONFIG_OBJECTS := $(patsubst $(CONFIG_DIR)/%.cpp,$(BUILD_DIR)/config/%.o,$(CONFIG_SOURCES))

PETSC_TEST_FLAGS := \
	$(COMMON_CXXFLAGS) \
	-I$(METIS_DIR)/include \
	-I$(METIS_DIR)/install/include

.PHONY: config

config: $(CONFIG_OBJECTS)
	@echo "[PASS] configuration matrix compiled: $(words $(CONFIG_OBJECTS))"

$(BUILD_DIR)/config/%.o: $(CONFIG_DIR)/%.cpp $(UNIT_PUBLIC_HEADERS) | dirs
	@echo "[CXX][config] $<"
	@$(PETSC_CXXCOMPILE_SINGLE) $(PETSC_TEST_FLAGS) -c $< -o $@

# ----------------------------------------------------------------------------
# 3. PETSc/MPI integration tests
# ----------------------------------------------------------------------------
INTEGRATION_TEST_NAMES := \
	common_petsc_test \
	structured_grid_test \
	distributed_mesh_test \
	cpgrid_test \
	natural_structured_test \
	natural_cpgrid_test \
	output_cpgrid_test

INTEGRATION_BINS := $(addprefix $(BIN_DIR)/,$(INTEGRATION_TEST_NAMES))
INTEGRATION_OBJECTS := $(patsubst %,$(BUILD_DIR)/integration/%.o,$(INTEGRATION_TEST_NAMES))

# The project currently has compiled sources only in common/ and grid/.
COMMON_PETSC_OBJECT := $(BUILD_DIR)/project/common_petsc_io.o
GRID_SOURCE_DIR := $(PROJECT_ROOT)/grid/src
GRID_SOURCES := $(wildcard $(GRID_SOURCE_DIR)/*.cpp)
GRID_PRIVATE_HEADERS := $(wildcard $(GRID_SOURCE_DIR)/*.hpp)
GRID_OBJECTS := $(patsubst $(GRID_SOURCE_DIR)/%.cpp,$(BUILD_DIR)/project/grid_%.o,$(GRID_SOURCES))
PROJECT_OBJECTS := $(COMMON_PETSC_OBJECT) $(GRID_OBJECTS)

# The PETSc build used on the cluster is customized and leaves two C++ symbols
# (updateState/updateSol) to the application.  A test-only weak shim satisfies
# ordinary PETSc executables; Natural/CpGrid provides strong real adapters.
PETSC_HOOKS_OBJECT := $(BUILD_DIR)/support/petsc_custom_hooks.o

METIS_LIB := \
	-L$(METIS_DIR)/install/lib \
	-lmetis \
	-Wl,-rpath,$(METIS_DIR)/install/lib

.PHONY: integration run-integration

integration: $(INTEGRATION_BINS)
	@echo "[PASS] PETSc/MPI integration tests compiled: $(words $(INTEGRATION_TEST_NAMES))"

$(BUILD_DIR)/integration/%.o: $(INTEGRATION_DIR)/%.cpp $(UNIT_PUBLIC_HEADERS) | dirs
	@echo "[CXX][integration] $<"
	@$(PETSC_CXXCOMPILE_SINGLE) $(PETSC_TEST_FLAGS) -c $< -o $@

$(PETSC_HOOKS_OBJECT): $(SUPPORT_DIR)/petsc_custom_hooks.cpp $(UNIT_PUBLIC_HEADERS) | dirs
	@echo "[CXX][support] $<"
	@$(PETSC_CXXCOMPILE_SINGLE) $(PETSC_TEST_FLAGS) -c $< -o $@

$(COMMON_PETSC_OBJECT): $(PROJECT_ROOT)/common/src/petsc_io.cpp $(UNIT_PUBLIC_HEADERS) | dirs
	@echo "[CXX][project] $<"
	@$(PETSC_CXXCOMPILE_SINGLE) $(PETSC_TEST_FLAGS) -c $< -o $@

$(BUILD_DIR)/project/grid_%.o: $(GRID_SOURCE_DIR)/%.cpp $(UNIT_PUBLIC_HEADERS) $(GRID_PRIVATE_HEADERS) | dirs
	@echo "[CXX][project] $<"
	@$(PETSC_CXXCOMPILE_SINGLE) $(PETSC_TEST_FLAGS) -c $< -o $@

$(BIN_DIR)/common_petsc_test: $(BUILD_DIR)/integration/common_petsc_test.o $(PROJECT_OBJECTS) $(PETSC_HOOKS_OBJECT)
	@echo "[LD] $@"
	@$(CXXLINKER) -o $@ $^ $(METIS_LIB) $(PETSC_LIB) $(LDFLAGS)

$(BIN_DIR)/structured_grid_test: $(BUILD_DIR)/integration/structured_grid_test.o $(PROJECT_OBJECTS) $(PETSC_HOOKS_OBJECT)
	@echo "[LD] $@"
	@$(CXXLINKER) -o $@ $^ $(METIS_LIB) $(PETSC_LIB) $(LDFLAGS)

$(BIN_DIR)/distributed_mesh_test: $(BUILD_DIR)/integration/distributed_mesh_test.o $(PROJECT_OBJECTS) $(PETSC_HOOKS_OBJECT)
	@echo "[LD] $@"
	@$(CXXLINKER) -o $@ $^ $(METIS_LIB) $(PETSC_LIB) $(LDFLAGS)

$(BIN_DIR)/cpgrid_test: $(BUILD_DIR)/integration/cpgrid_test.o $(PROJECT_OBJECTS) $(PETSC_HOOKS_OBJECT)
	@echo "[LD] $@"
	@$(CXXLINKER) -o $@ $^ $(METIS_LIB) $(PETSC_LIB) $(LDFLAGS)

$(BIN_DIR)/natural_structured_test: $(BUILD_DIR)/integration/natural_structured_test.o $(PROJECT_OBJECTS) $(PETSC_HOOKS_OBJECT)
	@echo "[LD] $@"
	@$(CXXLINKER) -o $@ $^ $(METIS_LIB) $(PETSC_LIB) $(LDFLAGS)

$(BIN_DIR)/natural_cpgrid_test: $(BUILD_DIR)/integration/natural_cpgrid_test.o $(PROJECT_OBJECTS) $(PETSC_HOOKS_OBJECT)
	@echo "[LD] $@"
	@$(CXXLINKER) -o $@ $^ $(METIS_LIB) $(PETSC_LIB) $(LDFLAGS)

$(BIN_DIR)/output_cpgrid_test: $(BUILD_DIR)/integration/output_cpgrid_test.o $(PROJECT_OBJECTS) $(PETSC_HOOKS_OBJECT)
	@echo "[LD] $@"
	@$(CXXLINKER) -o $@ $^ $(METIS_LIB) $(PETSC_LIB) $(LDFLAGS)

# Default DQ mesh used by the current MPMC_LH validation environment.
MESH_DIR ?= $(DQ_MESH_DIR)
NP ?= 1
PARTITION ?= mt_module
RUN_AND_LOG := ./scripts/run_and_log.sh

# yhrun needs an explicit partition on a cluster login node. Local mpiexec/
# mpirun must not receive Slurm/yhrun-specific arguments.
ifeq ($(MPMC_PLATFORM),hpc)
  ifeq ($(strip $(SLURM_JOB_ID)),)
    MPI_PARTITION_ARG := -p $(PARTITION)
  else
    MPI_PARTITION_ARG :=
  endif
else
  MPI_PARTITION_ARG :=
endif
TEST_LD_LIBRARY_PATH := $(LOCAL_LIB_DIR):$(PETSC_DIR)/$(PETSC_ARCH)/lib:$(VTK_INSTALL)/lib:$(METIS_DIR)/install/lib:$(X11_LIB_DIR):$(LD_LIBRARY_PATH)
MPI_RUN = env LD_LIBRARY_PATH="$(TEST_LD_LIBRARY_PATH)" $(MPI_LAUNCHER) $(MPI_PARTITION_ARG) $(MPI_LAUNCHER_ARGS) -n $(NP)

.PHONY: distributed-mesh run-distributed-mesh

distributed-mesh: $(BIN_DIR)/distributed_mesh_test
	@echo "[PASS] root-only distributed Mesh integration test compiled"

run-distributed-mesh: distributed-mesh
	@mkdir -p $(LOG_DIR)
	@echo "[RUN] root-only distributed Mesh"
	@$(RUN_AND_LOG) $(LOG_DIR)/distributed_mesh_test.log $(MPI_RUN) $(BIN_DIR)/distributed_mesh_test
	@echo "[PASS] root-only distributed Mesh integration test"

run-integration: integration
	@mkdir -p $(LOG_DIR) $(RESULT_DIR)
	@echo "[RUN] common PETSc I/O"
	@$(RUN_AND_LOG) $(LOG_DIR)/common_petsc_test.log $(MPI_RUN) $(BIN_DIR)/common_petsc_test
	@echo "[RUN] StructuredGrid"
	@$(RUN_AND_LOG) $(LOG_DIR)/structured_grid_test.log $(MPI_RUN) $(BIN_DIR)/structured_grid_test
	@echo "[RUN] root-only distributed Mesh"
	@$(RUN_AND_LOG) $(LOG_DIR)/distributed_mesh_test.log $(MPI_RUN) $(BIN_DIR)/distributed_mesh_test
	@echo "[RUN] Natural + StructuredGrid"
	@$(RUN_AND_LOG) $(LOG_DIR)/natural_structured_test.log $(MPI_RUN) $(BIN_DIR)/natural_structured_test
	@echo "[RUN] CpGrid: $(MESH_DIR)"
	@$(RUN_AND_LOG) $(LOG_DIR)/cpgrid_test.log $(MPI_RUN) $(BIN_DIR)/cpgrid_test -mesh_dir "$(MESH_DIR)"
	@echo "[RUN] Natural + CpGrid"
	@$(RUN_AND_LOG) $(LOG_DIR)/natural_cpgrid_test.log $(MPI_RUN) $(BIN_DIR)/natural_cpgrid_test -mesh_dir "$(MESH_DIR)"
	@echo "[RUN] Output + CpGrid"
	@$(RUN_AND_LOG) $(LOG_DIR)/output_cpgrid_test.log $(MPI_RUN) $(BIN_DIR)/output_cpgrid_test -mesh_dir "$(MESH_DIR)"
	@echo "============================================================"
	@echo "[PASS] all PETSc/MPI integration tests"
	@echo "============================================================"

# ----------------------------------------------------------------------------
# 4. Optional VTK integration test
# ----------------------------------------------------------------------------
VTK_TEST := output_vtk_test
VTK_OBJECT := $(BUILD_DIR)/integration/$(VTK_TEST).o
VTK_BIN := $(BIN_DIR)/$(VTK_TEST)

VTK_INCLUDE := -isystem $(VTK_INSTALL)/include/vtk-9.3
VTK_LIBS := \
	-L$(VTK_INSTALL)/lib \
	-lvtkCommonCore-9.3 \
	-lvtkCommonDataModel-9.3 \
	-lvtkCommonExecutionModel-9.3 \
	-lvtkFiltersCore-9.3 \
	-lvtkIOCore-9.3 \
	-lvtkIOXML-9.3 \
	-lvtkIOLegacy-9.3 \
	-lvtkParallelCore-9.3 \
	-lvtksys-9.3 \
	-Wl,-rpath,$(VTK_INSTALL)/lib \
	-Wl,-rpath,$(X11_LIB_DIR)

.PHONY: vtk run-vtk

vtk: $(VTK_BIN)
	@echo "[PASS] VTK integration test compiled"

$(VTK_OBJECT): $(INTEGRATION_DIR)/$(VTK_TEST).cpp | dirs
	@echo "[CXX][vtk] $<"
	@$(PETSC_CXXCOMPILE_SINGLE) $(PETSC_TEST_FLAGS) $(VTK_INCLUDE) -c $< -o $@

$(VTK_BIN): $(VTK_OBJECT) $(PROJECT_OBJECTS) $(PETSC_HOOKS_OBJECT)
	@echo "[LD] $@"
	@$(CXXLINKER) -o $@ $^ $(VTK_LIBS) $(METIS_LIB) $(PETSC_LIB) $(LDFLAGS)

run-vtk: vtk
	@mkdir -p $(LOG_DIR)
	@echo "[RUN] Output + VTK"
	@$(RUN_AND_LOG) $(LOG_DIR)/output_vtk_test.log $(MPI_RUN) $(VTK_BIN) -mesh_dir "$(MESH_DIR)"
	@echo "============================================================"
	@echo "[PASS] VTK integration test"
	@echo "============================================================"

# ----------------------------------------------------------------------------
# Aggregate targets
# ----------------------------------------------------------------------------
.PHONY: all full run-full dirs print-config clean

all: unit

full: unit config integration vtk
	@echo "============================================================"
	@echo "[PASS] all MPMC_SCW test programs compiled"
	@echo "============================================================"

run-full: full
	@$(MAKE) --no-print-directory run-unit
	@$(MAKE) --no-print-directory run-integration NP="$(NP)" MESH_DIR="$(MESH_DIR)" PARTITION="$(PARTITION)"
	@$(MAKE) --no-print-directory run-vtk NP="$(NP)" MESH_DIR="$(MESH_DIR)" PARTITION="$(PARTITION)"
	@echo "============================================================"
	@echo "[PASS] MPMC_SCW full test suite"
	@echo "============================================================"

# mkdir is intentionally centralized here; test sources never create build dirs.
dirs:
	@mkdir -p \
		$(BIN_DIR) \
		$(BUILD_DIR)/config \
		$(BUILD_DIR)/integration \
		$(BUILD_DIR)/project \
		$(BUILD_DIR)/support \
		$(LOG_DIR) \
		$(RESULT_DIR)

print-config:
	@echo "PROJECT_ROOT  = $(PROJECT_ROOT)"
	@echo "MPMC_PLATFORM = $(MPMC_PLATFORM)"
	@echo "PETSC_DIR    = $(PETSC_DIR)"
	@echo "PETSC_ARCH   = $(PETSC_ARCH)"
	@echo "FMT_DIR      = $(FMT_DIR)"
	@echo "METIS_DIR    = $(METIS_DIR)"
	@echo "VTK_INSTALL  = $(VTK_INSTALL)"
	@echo "X11_LIB_DIR  = $(X11_LIB_DIR)"
	@echo "LOCAL_LIB_DIR= $(LOCAL_LIB_DIR)"
	@echo "MPI_LAUNCHER = $(MPI_LAUNCHER)"
	@echo "MPI_ARGS     = $(MPI_LAUNCHER_ARGS)"
	@echo "PARTITION    = $(PARTITION)"
	@echo "PMI2_FILTER  = $${MPMC_FILTER_PMI2_WARNINGS:-1}"
	@echo "MESH_DIR     = $(MESH_DIR)"
	@echo "NP           = $(NP)"

clean::
	@rm -rf $(BUILD_DIR) $(BIN_DIR) $(LOG_DIR) $(RESULT_DIR)
	@rm -rf output_core_test_results output_feature_matrix_results output_well_test_results
	@rm -rf output_cpgrid_test_results
	@rm -f common_vector.m common_vector.bin common_matrix.m output_vtk_test.vtu
	@mkdir -p $(LOG_DIR) $(RESULT_DIR)
	@echo "[CLEAN] test build, logs and generated results removed"
