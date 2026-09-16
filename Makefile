# =============================================================================
# MPMC_SCW root command router
# =============================================================================
.DEFAULT_GOAL := help
SHELL := /bin/bash

CASE ?=
NP ?=
MESH_DIR ?=
EOS ?=
RESULT_DIR ?=
RUN_ARGS ?=
PACKAGE ?=
PACKAGE_ARGS ?=
CXX ?= g++

.PHONY: help doctor doctor-local doctor-hpc doctor-full print-config check-docs audit unit unit-gcc unit-clang cases case prepare reset-run run submit integration distributed-mesh \
        full tools tool-three-eos tool-public-benchmark clean \
        distclean package verify-package status

help:
	@echo "MPMC_SCW local/HPC commands"
	@echo ""
	@echo "  make doctor                         check the configured local/HPC profile"
	@echo "  make doctor-local                   check local Linux/WSL essentials"
	@echo "  make doctor-hpc                     check cluster scheduler/module essentials"
	@echo "  make doctor-full                    also require VTK/X11 full-test dependencies"
	@echo "  make print-config                   print resolved machine configuration"
	@echo "  make unit [CXX=g++]                 run 40 PETSc-free unit tests"
	@echo "  make unit-gcc                       clean GCC release gate"
	@echo "  make unit-clang                     run unit tests with clang++"
	@echo "  make cases -j                       compile all reservoir cases"
	@echo "  make case CASE=<name> -j            compile one case"
	@echo "  make prepare CASE=<name>            create case runtime dirs/run.sh"
	@echo "  make reset-run CASE=<name>          replace custom/legacy run.sh with managed template"
	@echo "  make run CASE=<name> NP=2 EOS=pr    compile + prepare + interactive local/HPC run"
	@echo "  make submit CASE=<name>             compile and sbatch one case"
	@echo "  make integration NP=2 MESH_DIR=<p>  run PETSc/MPI integration without VTK"
	@echo "  make distributed-mesh NP=4          run self-contained root-only Mesh MPI test"
	@echo "  make full NP=1 MESH_DIR=<path>      run full PETSc/MPI/VTK test gate"
	@echo "  make tools                           compile all offline tools"
	@echo "  make tool-three-eos                  run H2O/CO2/CH4/nC16 PR/SW/CPA scan"
	@echo "  make tool-public-benchmark           run independent public-data EOS/flash benchmark"
	@echo "  make check-docs                      verify local Markdown links"
	@echo "  make audit                           check repository cleanliness rules"
	@echo "  make clean                           remove normal generated files"
	@echo "  make package                         create Git-enabled internal zip"
	@echo "  make verify-package PACKAGE=<zip>    verify checksums, source and Git history"
	@echo ""
	@$(MAKE) --no-print-directory -C case list

doctor:
	@./scripts/hpc_doctor.sh auto

doctor-local:
	@./scripts/hpc_doctor.sh local

doctor-hpc:
	@./scripts/hpc_doctor.sh hpc

doctor-full:
	@./scripts/hpc_doctor.sh full

print-config:
	@$(MAKE) --no-print-directory -C case print-config

check-docs:
	@python3 scripts/check_markdown_links.py

audit:
	@./scripts/repo_audit.sh

unit:
	@$(MAKE) --no-print-directory -C test run-unit CXX="$(CXX)"

unit-gcc:
	@$(MAKE) --no-print-directory -C test clean
	@$(MAKE) --no-print-directory -C test run-unit CXX=g++

unit-clang:
	@$(MAKE) --no-print-directory -C test clean
	@$(MAKE) --no-print-directory -C test run-unit CXX=clang++

cases:
	@$(MAKE) --no-print-directory -C case all

case:
	@if [ -z "$(CASE)" ]; then echo "Usage: make case CASE=<name>"; exit 2; fi
	@$(MAKE) --no-print-directory -C case "$(CASE)"

prepare:
	@if [ -z "$(CASE)" ]; then echo "Usage: make prepare CASE=<name>"; exit 2; fi
	@$(MAKE) --no-print-directory -C case prepare CASE="$(CASE)"

reset-run:
	@if [ -z "$(CASE)" ]; then echo "Usage: make reset-run CASE=<name>"; exit 2; fi
	@$(MAKE) --no-print-directory -C case reset-run CASE="$(CASE)"

run:
	@if [ -z "$(CASE)" ]; then echo "Usage: make run CASE=<name> [NP=2] [EOS=pr] [RESULT_DIR=...] [MESH_DIR=...] [RUN_ARGS='...']"; exit 2; fi
	@$(MAKE) --no-print-directory -C case "$(CASE)"
	@$(MAKE) --no-print-directory -C case prepare CASE="$(CASE)"
	@cd "case/$(CASE)" && \
	  NP="$(NP)" EOS="$(EOS)" RESULT_DIR="$(RESULT_DIR)" MESH_DIR="$(MESH_DIR)" \
	  bash run.sh $(RUN_ARGS)

submit:
	@if [ -z "$(CASE)" ]; then echo "Usage: make submit CASE=<name>"; exit 2; fi
	@$(MAKE) --no-print-directory -C case submit CASE="$(CASE)"

integration:
	@if [ -z "$(MESH_DIR)" ]; then echo "Usage: make integration [NP=1] MESH_DIR=/path/to/DQ_data"; exit 2; fi
	@np="$(NP)"; if [ -z "$$np" ]; then np=1; fi; \
	  $(MAKE) --no-print-directory -C test run-integration NP="$$np" MESH_DIR="$(MESH_DIR)"

distributed-mesh:
	@np="$(NP)"; if [ -z "$$np" ]; then np=2; fi; \
	  $(MAKE) --no-print-directory -C test run-distributed-mesh NP="$$np"

full:
	@if [ -z "$(MESH_DIR)" ]; then echo "Usage: make full [NP=1] MESH_DIR=/path/to/DQ_data"; exit 2; fi
	@np="$(NP)"; if [ -z "$$np" ]; then np=1; fi; \
	  $(MAKE) --no-print-directory -C test run-full NP="$$np" MESH_DIR="$(MESH_DIR)"

tools:
	@$(MAKE) --no-print-directory -C tools all CXX="$(CXX)"

tool-three-eos:
	@$(MAKE) --no-print-directory -C tools run-three-eos CXX="$(CXX)"

tool-public-benchmark:
	@$(MAKE) --no-print-directory -C tools run-public-benchmark CXX="$(CXX)"

clean:
	@$(MAKE) --no-print-directory -C test clean
	@$(MAKE) --no-print-directory -C tools clean
	@$(MAKE) --no-print-directory -C case clean-all

distclean: clean
	@find tools/example -type d -name reference_output -prune -exec rm -rf {} +
	@find . -type d -name __pycache__ -prune -exec rm -rf {} +
	@find . -type f \( -name '*.pyc' -o -name '*.pyo' \) -delete

package:
	@if [ -n "$(PACKAGE)" ]; then \
	  ./scripts/package_internal.sh $(PACKAGE_ARGS) "$(PACKAGE)"; \
	else \
	  ./scripts/package_internal.sh $(PACKAGE_ARGS); \
	fi

verify-package:
	@if [ -z "$(PACKAGE)" ]; then echo "Usage: make verify-package PACKAGE=/path/to/package.zip"; exit 2; fi
	@./scripts/verify_handoff_package.sh "$(PACKAGE)"

status:
	@git status --short --branch
	@git describe --tags --always --dirty
