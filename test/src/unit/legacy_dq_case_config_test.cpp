/**
 * @file legacy_dq_case_config_test.cpp
 * @brief DQ 角点网格正式算例的 PETSc-free 配置与流体系统编译门禁。
 */
#include "../../../case/DQcase/case_config.hpp"
#include "../../../case/DQcase/well_config.hpp"

#include <case/well_factory.hpp>
#include <case/fluid_factory.hpp>
#include <indices/indices.hpp>

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace
{
using ModelConfig = MPMC::CompositionalModelConfig<
    CaseConfig::Model::numberOfComponents,
    CaseConfig::Model::hasWater,
    CaseConfig::Model::hasWells,
    CaseConfig::Model::enableDissolution,
    CaseConfig::Model::enableAdsorption,
    CaseConfig::Model::enableLandTrapping>;
using Indices = MPMC::ScalarIndices<ModelConfig>;

void require(bool condition, const char *message)
{
    if (!condition)
        throw std::runtime_error(message);
}
}

int main()
{
    try
    {
        {
            std::vector<MPMC::WellPerforation<int>> perforations{{0, 1.0e-12}};
            const auto well = MPMC::cases::makeWell<Indices, int>(
                WellConfig::wells.front(), 0, std::move(perforations));
            require(well.id == WellConfig::wells.front().id,
                    "case well factory must preserve the configured well id");
            require(well.control == MPMC::cases::detail::toNaturalControl(
                        WellConfig::wells.front().control),
                    "case well factory must preserve the configured well control");
        }

        static_assert(Indices::numComponents == 8);
        static_assert(Indices::hasWater);
        static_assert(Indices::hasWellUnknown);
        static_assert(Indices::hasAqueousCO2Dissolution);
        static_assert(Indices::hasLandTrapping);

        auto fluid = MPMC::cases::makeFluidSystem<Indices, CaseConfig::Config>();
        const auto phase = fluid.eos.phaseResult(
            CaseConfig::InitialState::pressure,
            CaseConfig::Fluid::temperature,
            CaseConfig::InitialState::oilComposition,
            MPMC::CompositionalPhase::Oil);
        require(std::isfinite(phase.compressibility) && phase.compressibility > 0.0,
                "DQ EOS phase result must remain finite");
        const double rho = fluid.eos.molarDensity(
            CaseConfig::InitialState::pressure,
            CaseConfig::Fluid::temperature,
            CaseConfig::InitialState::oilComposition,
            phase.compressibility);
        require(std::isfinite(rho) && rho > 0.0,
                "DQ EOS molar density must remain finite");
        std::cout << "Legacy DQ case configuration: ALL PASS\n";
        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Legacy DQ case configuration: FAIL: " << e.what() << '\n';
        return 1;
    }
}
