# ============================================================================
# MPMC_SCW tools module
#
# Pure C++ offline utilities that reuse production thermodynamic/model kernels
# without requiring PETSc/MPI. Generated data are written below
# example/<name>/results and are never required by the reservoir simulator.
# ============================================================================
.DEFAULT_GOAL := all
SHELL := /bin/bash

PROJECT_ROOT := ..
CXX ?= g++
PYTHON ?= python3

MPMC_INCLUDES := \
	-I$(PROJECT_ROOT)/ad/include \
	-I$(PROJECT_ROOT)/common/include \
	-I$(PROJECT_ROOT)/indices/include \
	-I$(PROJECT_ROOT)/models/include \
	-I$(PROJECT_ROOT)/tools/include

CXXFLAGS := -std=c++17 -O2 -Wall -Wextra -Wpedantic -Werror $(MPMC_INCLUDES)

BIN_DIR := bin
EXAMPLES := pr_c2_c3_nc5 h2o_co2_ch4_nc16_three_eos three_eos_public_benchmark five_component_common_fit cpa_dme_co2_water_vlle heringer2025_bsb_three_phase scw_binary_calibration scw_pseudocomponent_binary_regression scw_kerogen_0d_pvt_acceptance
EXAMPLE_BINS := $(addprefix $(BIN_DIR)/,$(EXAMPLES))

PUBLIC_HEADERS := $(shell find include $(PROJECT_ROOT)/models/include $(PROJECT_ROOT)/indices/include \
	$(PROJECT_ROOT)/common/include $(PROJECT_ROOT)/ad/include -type f \( -name '*.hpp' -o -name '*.h' \) 2>/dev/null)

.PHONY: all examples example run-example run-three-eos run-public-benchmark run-public-phase-diagrams redraw-public-phase-diagrams run-five-component-common-fit run-cpa-dme-vlle run-heringer2025-bsb run-scw-binary-calibration run-scw-pseudocomponent-binary-regression run-scw-kerogen-0d-pvt-acceptance require-scw-kerogen-0d-pvt-acceptance plot-cpa-dme-vlle plot-cpa-dme-vlle-by-pressure reference reference-three-eos clean

all: examples

examples: $(EXAMPLE_BINS)
	@echo "[PASS] tools examples compiled: $(EXAMPLES)"

# Backward-compatible target for the original v32 example.
example: $(BIN_DIR)/pr_c2_c3_nc5
	@echo "[PASS] tools example compiled: $<"

$(BIN_DIR)/%: example/%/main.cpp $(PUBLIC_HEADERS)
	@mkdir -p $(BIN_DIR)
	@echo "[CXX][tools] $<"
	@$(CXX) $(CXXFLAGS) $< -o $@

run-example: $(BIN_DIR)/pr_c2_c3_nc5
	@rm -rf example/pr_c2_c3_nc5/results
	@mkdir -p example/pr_c2_c3_nc5/results
	@$(BIN_DIR)/pr_c2_c3_nc5 example/pr_c2_c3_nc5/results
	@echo "[PASS] v32 hydrocarbon example data generated"

run-three-eos: $(BIN_DIR)/h2o_co2_ch4_nc16_three_eos
	@rm -rf example/h2o_co2_ch4_nc16_three_eos/results
	@mkdir -p example/h2o_co2_ch4_nc16_three_eos/results
	@$(BIN_DIR)/h2o_co2_ch4_nc16_three_eos example/h2o_co2_ch4_nc16_three_eos/results
	@echo "[PASS] H2O/CO2/CH4/nC16 three-EOS example data generated"

run-public-benchmark: $(BIN_DIR)/three_eos_public_benchmark
	@rm -rf example/three_eos_public_benchmark/results
	@mkdir -p example/three_eos_public_benchmark/results
	@$(PYTHON) example/three_eos_public_benchmark/run_benchmark.py \
		--mpmc $(CURDIR)/$(BIN_DIR)/three_eos_public_benchmark \
		--reference-dir $(CURDIR)/example/three_eos_public_benchmark/reference_data \
		--output-dir $(CURDIR)/example/three_eos_public_benchmark/results
	@echo "[PASS] public three-EOS benchmark completed"

run-public-phase-diagrams: $(BIN_DIR)/three_eos_public_benchmark $(BIN_DIR)/h2o_co2_ch4_nc16_three_eos
	@mkdir -p example/three_eos_public_benchmark/results/phase_diagrams
	@$(PYTHON) example/three_eos_public_benchmark/run_phase_diagrams.py \
		--mpmc $(CURDIR)/$(BIN_DIR)/three_eos_public_benchmark \
		--same-fluid-executable $(CURDIR)/$(BIN_DIR)/h2o_co2_ch4_nc16_three_eos \
		--reference-dir $(CURDIR)/example/three_eos_public_benchmark/reference_data \
		--output-dir $(CURDIR)/example/three_eos_public_benchmark/results/phase_diagrams
	@echo "[PASS] public three-EOS phase diagrams completed"

redraw-public-phase-diagrams:
	@$(PYTHON) example/three_eos_public_benchmark/run_phase_diagrams.py \
		--reference-dir $(CURDIR)/example/three_eos_public_benchmark/reference_data \
		--output-dir $(CURDIR)/example/three_eos_public_benchmark/results/phase_diagrams \
		--redraw-only
	@echo "[PASS] public three-EOS phase diagrams redrawn from saved production results"

run-five-component-common-fit: $(BIN_DIR)/five_component_common_fit
	@$(PYTHON) example/five_component_common_fit/prepare_dataset.py \
		--reference-dir $(CURDIR)/example/three_eos_public_benchmark/reference_data \
		--output-dir $(CURDIR)/example/five_component_common_fit/data
	@rm -rf example/five_component_common_fit/results
	@mkdir -p example/five_component_common_fit/results
	@$(BIN_DIR)/five_component_common_fit \
		example/five_component_common_fit/data/common_calibration_data.csv \
		example/five_component_common_fit/results
	@echo "[PASS] common-data PR/SW/CPA parameter fits completed"

run-cpa-dme-vlle: $(BIN_DIR)/cpa_dme_co2_water_vlle
	@rm -rf example/cpa_dme_co2_water_vlle/results
	@mkdir -p example/cpa_dme_co2_water_vlle/results
	@$(BIN_DIR)/cpa_dme_co2_water_vlle \
		example/cpa_dme_co2_water_vlle/reference_data/laursen2002_co2_dme_water_vlle.csv \
		example/cpa_dme_co2_water_vlle/results
	@echo "[PASS] public CO2/DME/H2O CPA VLLE example completed"

run-heringer2025-bsb: $(BIN_DIR)/heringer2025_bsb_three_phase
	@rm -rf example/heringer2025_bsb_three_phase/results
	@mkdir -p example/heringer2025_bsb_three_phase/results
	@$(BIN_DIR)/heringer2025_bsb_three_phase \
		example/heringer2025_bsb_three_phase/results
	@echo "[PASS] Heringer 2025 H2O/CO2/BSB-oil three-phase benchmark completed"

run-scw-binary-calibration: $(BIN_DIR)/scw_binary_calibration
	@test ! -e example/scw_binary_calibration/results || \
		(echo "Refusing to overwrite example/scw_binary_calibration/results"; exit 2)
	@$(BIN_DIR)/scw_binary_calibration \
		../case/scw_kerogen_common/calibration/squalane_target_window.csv \
		example/scw_binary_calibration/reference_data/squalane_pvt_viscosity.csv \
		example/scw_binary_calibration/results
	@echo "[PASS] zero-flow SCW binary PVT/BIP/viscosity calibration completed"

run-scw-pseudocomponent-binary-regression: $(BIN_DIR)/scw_pseudocomponent_binary_regression
	@test ! -e example/scw_pseudocomponent_binary_regression/results || \
		(echo "Refusing to overwrite example/scw_pseudocomponent_binary_regression/results"; exit 2)
	@set +e; \
	$(BIN_DIR)/scw_pseudocomponent_binary_regression \
		../case/scw_kerogen_lumped_flow_study/binary_pr_calibration/binary_systems.csv \
		../case/scw_kerogen_lumped_flow_study/binary_pr_calibration/binary_observations.csv \
		example/scw_pseudocomponent_binary_regression/results; \
	status=$$?; \
	if [ $$status -eq 3 ]; then \
		echo "[PASS] generic binary harness executed; reservoir gate remains blocked by design"; \
		exit 0; \
	fi; \
	exit $$status


run-scw-kerogen-0d-pvt-acceptance: $(BIN_DIR)/scw_kerogen_0d_pvt_acceptance
	@rm -rf example/scw_kerogen_0d_pvt_acceptance/results
	@$(BIN_DIR)/scw_kerogen_0d_pvt_acceptance \
		../case/scw_kerogen_lumped_flow_study \
		example/scw_kerogen_0d_pvt_acceptance/results
	@echo "[PASS] PR and CPA 0D PVT acceptance gate passed"

require-scw-kerogen-0d-pvt-acceptance: run-scw-kerogen-0d-pvt-acceptance
	@echo "[PASS] hard prerequisite satisfied: zero-flow PVT is self-consistent"

plot-cpa-dme-vlle: run-cpa-dme-vlle
	@$(PYTHON) example/cpa_dme_co2_water_vlle/plot_results.py \
		--input example/cpa_dme_co2_water_vlle/results/phase_comparison.csv \
		--output-dir example/cpa_dme_co2_water_vlle/results/figures
	@echo "[PASS] bilingual CPA VLLE parity figures generated"

# Preserve the existing single-panel figure and write the pressure-resolved
# comparison under a separate stem.
plot-cpa-dme-vlle-by-pressure: $(BIN_DIR)/cpa_dme_co2_water_vlle
	@mkdir -p example/cpa_dme_co2_water_vlle/results
	@$(BIN_DIR)/cpa_dme_co2_water_vlle \
		example/cpa_dme_co2_water_vlle/reference_data/laursen2002_co2_dme_water_vlle.csv \
		example/cpa_dme_co2_water_vlle/results
	@$(PYTHON) example/cpa_dme_co2_water_vlle/plot_results_by_pressure.py \
		--input example/cpa_dme_co2_water_vlle/results/phase_comparison.csv \
		--output-dir example/cpa_dme_co2_water_vlle/results/figures
	@echo "[PASS] pressure-resolved CPA VLLE parity figures generated"

reference: $(BIN_DIR)/pr_c2_c3_nc5
	@rm -rf example/pr_c2_c3_nc5/reference_output
	@mkdir -p example/pr_c2_c3_nc5/reference_output
	@$(BIN_DIR)/pr_c2_c3_nc5 example/pr_c2_c3_nc5/reference_output
	@echo "[PASS] hydrocarbon reference phase-diagram data regenerated"

reference-three-eos: $(BIN_DIR)/h2o_co2_ch4_nc16_three_eos
	@rm -rf example/h2o_co2_ch4_nc16_three_eos/reference_output
	@mkdir -p example/h2o_co2_ch4_nc16_three_eos/reference_output
	@$(BIN_DIR)/h2o_co2_ch4_nc16_three_eos example/h2o_co2_ch4_nc16_three_eos/reference_output
	@echo "[PASS] three-EOS aqueous reference data regenerated"

clean:
	@rm -rf $(BIN_DIR) example/pr_c2_c3_nc5/results example/h2o_co2_ch4_nc16_three_eos/results example/three_eos_public_benchmark/results example/five_component_common_fit/results example/cpa_dme_co2_water_vlle/results example/heringer2025_bsb_three_phase/results example/scw_binary_calibration/results example/scw_pseudocomponent_binary_regression/results example/scw_kerogen_0d_pvt_acceptance/results
