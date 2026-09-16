/**
 * @file scw_kerogen_squalane_2d_case_test.cpp
 * @brief 验证水–角鲨烷二维算例的网格、井位、BIP 和初态闪蒸。
 */
#include "../../../case/scw_kerogen_squalane_2d/case_config.hpp"
#include "../../../case/scw_kerogen_squalane_2d/well_config.hpp"

#include <case/case_validation.hpp>
#include <case/fluid_factory.hpp>
#include <indices/indices.hpp>
#include <natural/thermo/three_phase_flash.hpp>

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace
{
using Config = ScwKerogenSqualane2D::Config;
using ModelConfig = MPMC::CompositionalModelConfig<
    Config::Model::numberOfComponents, Config::Model::hasWater,
    Config::Model::hasWells, Config::Model::enableDissolution,
    Config::Model::enableAdsorption, Config::Model::enableLandTrapping,
    Config::Model::phaseBehavior>;
using Indices = MPMC::ScalarIndices<ModelConfig>;

void require(bool value, const char *message)
{
    if (!value) throw std::runtime_error(message);
}

template <class FactoryConfig>
void checkFlash()
{
    auto fluid = MPMC::cases::makeFluidSystem<Indices, FactoryConfig>();
    MPMC::ThreePhaseFlashOptions options;
    options.waterComponent = Config::Fluid::waterComponent;
    const MPMC::CubicThreePhaseFlash<Indices> flash(fluid.eos, options);
    const auto result = flash.flash(
        Config::InitialState::pressure, Config::InitialState::temperature,
        Config::InitialState::overallComposition);
    std::cout << "[INFO] squalane initial phase mask="
              << static_cast<int>(result.presence.bits()) << '\n';
    require(result.converged, "2D squalane initial flash must converge");
    require(result.presence.count() == 1 &&
                result.presence.contains(MPMC::CompositionalPhase::Oil),
            "all EOS variants must share an oil-rich single-phase initial state");
}

} // namespace

int main()
{
    try
    {
        MPMC::cases::validateCaseConfig<Indices, Config>();
        static_assert(Config::Grid::nx == 60 && Config::Grid::ny == 20 && Config::Grid::nz == 1);
        static_assert(ScwKerogen2D::primaryCentreJ == 9);
        require(std::abs(ScwKerogen2D::centreOffset - 0.0025) < 1.0e-15,
                "even-grid centre offset must be half a y cell");
        require(WellConfig::wells[0].completion.j == 9 &&
                    WellConfig::wells[1].completion.j == 9,
                "both wells must be centred on opposite short edges");
        require(std::abs(
                    ScwKerogenCalibration::prSqualaneKij(653.2)
                    - 0.053233659543451335) < 1.0e-12,
                "squalane BIP must use the target tie-line composition fit");
        checkFlash<ScwKerogenSqualane2D::PrFactoryConfig>();
        checkFlash<ScwKerogenSqualane2D::SwFactoryConfig>();
        checkFlash<ScwKerogenSqualane2D::CpaFactoryConfig>();
        std::cout << "[PASS] scw_kerogen_squalane_2d_case_test\n";
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << "[FAIL] scw_kerogen_squalane_2d_case_test: "
                  << error.what() << '\n';
        return 1;
    }
}
