/**
 * @file scw_kerogen_lmh_case_test.cpp
 * @brief 验证超临界水驱轻中重拟组分算例的配置与初始闪蒸。
 */
#include "../../../case/scw_kerogen_lmh_1d/case_config.hpp"
#include "../../../case/scw_kerogen_lmh_1d/well_config.hpp"

#include <case/case_validation.hpp>
#include <case/fluid_factory.hpp>
#include <indices/indices.hpp>
#include <natural/thermo/three_phase_flash.hpp>

#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace
{
using ModelConfig = MPMC::CompositionalModelConfig<
    CaseConfig::Model::numberOfComponents, CaseConfig::Model::hasWater,
    CaseConfig::Model::hasWells, CaseConfig::Model::enableDissolution,
    CaseConfig::Model::enableAdsorption, CaseConfig::Model::enableLandTrapping,
    CaseConfig::Model::phaseBehavior>;
using Indices = MPMC::ScalarIndices<ModelConfig>;

void require(bool condition, const char *message)
{
    if (!condition) throw std::runtime_error(message);
}
}

int main()
{
    try
    {
        static_assert(Indices::fullyCompositionalThreePhase);
        static_assert(Indices::numComponents == 4);
        MPMC::cases::validateCaseConfig<Indices, CaseConfig::Config>();
        require(CaseConfig::InitialState::temperature >
                    ScwKerogen1D::waterCriticalTemperature,
                "temperature must be supercritical for water");
        require(WellConfig::wells[0].injectedComponent ==
                    CaseConfig::Fluid::waterComponent,
                "injector must carry pure H2O");

        auto fluid = MPMC::cases::makeFluidSystem<
            Indices, CaseConfig::Config>();
        MPMC::ThreePhaseFlashOptions options;
        options.waterComponent = CaseConfig::Fluid::waterComponent;
        MPMC::CubicThreePhaseFlash<Indices> flash(fluid.eos, options);
        const auto result = flash.flash(
            CaseConfig::InitialState::pressure,
            CaseConfig::InitialState::temperature,
            CaseConfig::InitialState::overallComposition);
        std::cout << "[INFO] PR initial phase mask="
                  << static_cast<int>(result.presence.bits()) << '\n';
        require(result.converged, "initial H2O/LMH flash must converge");
        require(result.presence.count() == 1 &&
                    result.presence.contains(MPMC::CompositionalPhase::Oil),
                "LMH baseline must start oil-rich and single phase");

        std::array<double, Indices::numComponents> reconstructed{};
        for (std::size_t phase = 0; phase < result.phaseMoleFraction.size(); ++phase)
        {
            if (result.phaseMoleFraction[phase] <= options.phaseFractionTolerance)
                continue;
            const auto role = static_cast<MPMC::CompositionalPhase>(phase);
            const auto properties = fluid.eosFlowProperties(
                CaseConfig::InitialState::pressure,
                result.composition[phase], result.compressibility[phase],
                CaseConfig::InitialState::temperature, role);
            require(std::isfinite(std::get<2>(properties)) &&
                        std::get<2>(properties) > 0.0,
                    "active-phase viscosity must be finite and positive");
            for (std::size_t component = 0;
                 component < reconstructed.size(); ++component)
                reconstructed[component] += result.phaseMoleFraction[phase] *
                    result.composition[phase][component];
        }
        for (std::size_t component = 0; component < reconstructed.size(); ++component)
            require(std::abs(reconstructed[component] -
                        CaseConfig::InitialState::overallComposition[component]) < 1.0e-9,
                    "initial flash material balance failed");

        std::cout << "[PASS] scw_kerogen_lmh_case_test\n";
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << "[FAIL] scw_kerogen_lmh_case_test: "
                  << error.what() << '\n';
        return 1;
    }
}
