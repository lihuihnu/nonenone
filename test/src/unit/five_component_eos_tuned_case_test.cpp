/**
 * @file five_component_eos_tuned_case_test.cpp
 * @brief 新五组分优选参数 PR/SW/CPA 流动算例的配置与初始闪蒸预检。
 */
#include "../../../case/five_component_eos_tuned_compare/case_config.hpp"
#include "../../../case/five_component_eos_tuned_compare/well_config.hpp"

#include <case/case_validation.hpp>
#include <case/fluid_factory.hpp>
#include <indices/indices.hpp>
#include <natural/thermo/three_phase_flash.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

namespace
{
using ModelConfig = MPMC::CompositionalModelConfig<
    CaseConfig::Model::numberOfComponents,
    CaseConfig::Model::hasWater,
    CaseConfig::Model::hasWells,
    CaseConfig::Model::enableDissolution,
    CaseConfig::Model::enableAdsorption,
    CaseConfig::Model::enableLandTrapping,
    CaseConfig::Model::phaseBehavior>;
using Indices = MPMC::ScalarIndices<ModelConfig>;

void require(bool condition, const std::string &message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void near(double actual, double expected, double tolerance, const std::string &message)
{
    const double scale = std::max({1.0, std::abs(actual), std::abs(expected)});
    if (std::abs(actual - expected) > tolerance * scale)
        throw std::runtime_error(message);
}

template <class FactoryConfig>
auto flashInitial()
{
    auto fluid = MPMC::cases::makeFluidSystem<Indices, FactoryConfig>();
    MPMC::ThreePhaseFlashOptions options;
    options.waterComponent = CaseConfig::CommonFluid::waterComponent;
    const MPMC::CubicThreePhaseFlash<Indices> flash(fluid.eos, options);
    return flash.flash(
        CaseConfig::InitialState::pressure,
        CaseConfig::InitialState::temperature,
        CaseConfig::InitialState::overallComposition);
}

void checkConfiguration()
{
    static_assert(Indices::fullyCompositionalThreePhase);
    static_assert(CaseConfig::Model::numberOfComponents == 5);
    static_assert(CaseConfig::Grid::nx == 20 && CaseConfig::Grid::ny == 20 &&
                  CaseConfig::Grid::nz == 5);
    static_assert(CaseConfig::Time::adaptive);
    static_assert(CaseConfig::Time::numberOfSteps == 120);
    static_assert(CaseConfig::PrFluid::thermodynamicModel ==
                  MPMC::CubicThermodynamicModel::PengRobinson);
    static_assert(CaseConfig::SwFluid::thermodynamicModel ==
                  MPMC::CubicThermodynamicModel::SoreideWhitson);
    static_assert(CaseConfig::CpaFluid::thermodynamicModel ==
                  MPMC::CubicThermodynamicModel::CubicPlusAssociation);

    MPMC::cases::validateCaseConfig<Indices, CaseConfig::Config>();
    near(CaseConfig::SwFluid::soreideWhitsonSalinityMolality, 1.0, 1.0e-14,
         "SW salinity package");
    near(CaseConfig::PrFluid::componentVolumeTranslation[0],
         3.252997096253869e-6, 1.0e-14,
         "PR IAPWS-anchored water volume translation");
    near(CaseConfig::SwFluid::componentVolumeTranslation[0],
         3.267331277345136e-6, 1.0e-14,
         "SW IAPWS-anchored water volume translation");
    near(CaseConfig::SwFluid::soreideWhitsonAqueousVolumeTranslation[1],
         1.469715763209903e-6, 1.0e-14,
         "SW Garcia dissolved-CO2 aqueous density correction");
    near(CaseConfig::CpaFluid::cpaB[0], 14.52e-6, 1.0e-13,
         "CPA 4C-water co-volume");
    near(CaseConfig::CpaFluid::cpaAssociationVolume[0], 0.0692, 1.0e-13,
         "CPA 4C-water association volume");
    near(CaseConfig::CpaFluid::binaryInteraction[0][1],
         0.1584126867280001, 1.0e-13,
         "CPA fitted water-CO2 binary parameter");
    near(CaseConfig::CpaFluid::binaryInteraction[0][2],
         0.0204118168704000, 1.0e-13,
         "CPA fitted water-methane binary parameter");
    near(CaseConfig::CpaFluid::binaryInteraction[0][3], 0.04415, 1.0e-13,
         "CPA water-ethane binary parameter");
    near(CaseConfig::CpaFluid::binaryInteraction[0][4], 0.08750, 1.0e-13,
         "CPA water-n-butane binary parameter");
    near(CaseConfig::CpaFluid::cpaCrossAssociationEnergy[0][1], 14200.0,
         1.0e-13, "CPA B2 water-CO2 cross-association energy");
    near(CaseConfig::CpaFluid::cpaCrossAssociationVolume[0][1], 0.0162,
         1.0e-13, "CPA B2 water-CO2 cross-association volume");
    require(CaseConfig::CpaFluid::cpaAcceptorSites[1] == 1,
            "CPA B2 package requires one CO2 acceptor site");

    const double swCo2 = CaseConfig::SwFluid::soreideWhitsonAqueousWaterBip(
        CaseConfig::CommonFluid::co2Component,
        CaseConfig::InitialState::temperature,
        CaseConfig::SwFluid::soreideWhitsonSalinityMolality);
    near(swCo2,
         MPMC::SoreideWhitsonCorrelations::co2AqueousBipChabab2019(305.0, 1.0)
             + CaseConfig::SwFluid::co2AqueousBipOffset,
         1.0e-14,
         "SW must use the Chabab-2019 aqueous CO2 BIP plus fitted offset");
}

void checkGeologyAndWells()
{
    int channelCells = 0;
    int baffleCells = 0;
    int windowCells = 0;
    for (int k = 0; k < CaseConfig::Grid::nz; ++k)
        for (int j = 0; j < CaseConfig::Grid::ny; ++j)
            for (int i = 0; i < CaseConfig::Grid::nx; ++i)
            {
                if (CaseConfig::Rock::inBaffleWindow(i, j, k))
                    ++windowCells;
                else if (k == 2)
                    ++baffleCells;
                else if (CaseConfig::Rock::inChannel(i, j, k))
                    ++channelCells;
            }

    require(channelCells == 240, "straight channel must contain 240 cells");
    require(baffleCells == 384, "shale baffle must contain 384 cells");
    require(windowCells == 16, "single flow window must contain 16 cells");
    near(CaseConfig::Rock::kx(1, 4, 0) / CaseConfig::mD, 250.0, 1.0e-14,
         "injector-side channel permeability");
    near(CaseConfig::Rock::kz(0, 0, 2) / CaseConfig::mD, 0.05, 1.0e-14,
         "baffle vertical permeability");
    near(CaseConfig::Rock::kz(9, 9, 2) / CaseConfig::mD, 10.0, 1.0e-14,
         "window vertical permeability");

    require(WellConfig::wells.size() == 2,
            "case requires one injector and one producer");
    near(WellConfig::wells[0].target * CaseConfig::secondsPerDay,
         100000.0, 1.0e-14, "CO2 injector reference-volume rate");
    require(WellConfig::wells[0].completion.kCount == 2 &&
                WellConfig::wells[1].completion.kCount == 2,
            "each well must complete two layers");
}

void checkInitialFlash()
{
    const auto pr = flashInitial<CaseConfig::PrFactoryConfig>();
    const auto sw = flashInitial<CaseConfig::SwFactoryConfig>();
    const auto cpa = flashInitial<CaseConfig::CpaFactoryConfig>();
    require(pr.converged && sw.converged && cpa.converged,
            "all tuned EOS packages must flash the common initial P-T-z state");

    for (const auto *result : {&pr, &sw, &cpa})
    {
        require(result->phaseMoleFraction[0] > 1.0e-5 &&
                    result->phaseMoleFraction[1] > 1.0e-5 &&
                    result->phaseMoleFraction[2] > 1.0e-5,
                "initial state must contain oil-, gas- and water-rich phases");
        near(result->saturation[0] + result->saturation[1] + result->saturation[2],
             1.0, 1.0e-11, "initial saturation closure");
    }

    std::cout << "initial So/Sg/Sw\n"
              << "  PR : " << pr.saturation[0] << ' ' << pr.saturation[1] << ' '
              << pr.saturation[2] << '\n'
              << "  SW : " << sw.saturation[0] << ' ' << sw.saturation[1] << ' '
              << sw.saturation[2] << '\n'
              << "  CPA: " << cpa.saturation[0] << ' ' << cpa.saturation[1] << ' '
              << cpa.saturation[2] << '\n';
}

void checkCpaReducedAssociationPath()
{
    auto fluid = MPMC::cases::makeFluidSystem<
        Indices, CaseConfig::CpaFactoryConfig>();
    fluid.eos.resetThermodynamicProfile();
    const auto phase = fluid.eos.phaseResult(
        CaseConfig::InitialState::pressure,
        CaseConfig::InitialState::temperature,
        CaseConfig::InitialState::overallComposition,
        MPMC::CompositionalPhase::Water,
        false);
    require(std::isfinite(phase.compressibility),
            "CPA B2 phase evaluation must remain finite");

    const auto profile = fluid.eos.thermodynamicProfile();
    require(profile.cpaAssociationIterativeCalls > 0,
            "CPA B2 evaluation must exercise the reduced association path");
    const double averageIterations =
        static_cast<double>(profile.cpaAssociationIterations) /
        static_cast<double>(profile.cpaAssociationIterativeCalls);
    require(averageIterations < 32.0,
            "single-donor CPA site solve must not regress to slow fixed-point work");
}
} // namespace

int main()
{
    try
    {
        checkConfiguration();
        checkGeologyAndWells();
        checkInitialFlash();
        checkCpaReducedAssociationPath();
        std::cout << "[PASS] five_component_eos_tuned_case_test\n";
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << "[FAIL] five_component_eos_tuned_case_test: "
                  << error.what() << '\n';
        return 1;
    }
}
