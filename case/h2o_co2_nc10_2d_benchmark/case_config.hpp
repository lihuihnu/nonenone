/**
 * @file case_config.hpp
 * @brief 新热力学二维算例的 PR、SW 与 PR-CPA 配置选择。
 */
#pragma once

#include "benchmark_common.hpp"

#include <indices/model_config.hpp>
#include <natural/thermo/cubic_plus_association.hpp>
#include <natural/thermo/soreide_whitson.hpp>
#include <natural/thermo/thermodynamic_model.hpp>

#include <array>
#include <cstddef>

namespace CaseConfig
{

#ifndef MPMC_H2O_CO2_NC10_CASE_NAME
#define MPMC_H2O_CO2_NC10_CASE_NAME "h2o_co2_nc10_2d_benchmark"
#endif
inline constexpr char name[] = MPMC_H2O_CO2_NC10_CASE_NAME;
inline constexpr double gasConstant = 8.31446261815324;

struct Model
{
    static constexpr int numberOfComponents = 3;
    static constexpr bool hasWater = true;
    static constexpr bool hasWells = true;
    static constexpr bool enableDissolution = false;
    static constexpr bool enableAdsorption = false;
    static constexpr bool enableLandTrapping = false;
    static constexpr auto phaseBehavior =
        MPMC::PhaseBehaviorModel::FullyCompositionalThreePhase;
};

using Grid = BenchmarkCommon::Grid;
using Rock = BenchmarkCommon::Rock;

struct CommonFluid
{
    static constexpr std::size_t N = Model::numberOfComponents;
    static constexpr int waterComponent = 0;
    static constexpr int co2Component = 1;
    static constexpr int nc10Component = 2;
    inline static constexpr std::array<const char *, N> componentNames{
        "H2O", "CO2", "nC10"};
    inline static constexpr std::array<double, N> criticalTemperature{
        647.096, 304.1282, 617.70};
    inline static constexpr std::array<double, N> criticalPressure{
        22.064e6, 7.3773e6, 2.103e6};
    inline static constexpr std::array<double, N> criticalVolume{
        5.5948e-5, 9.4118e-5, 6.10e-4};
    inline static constexpr std::array<double, N> acentricFactor{
        0.3443, 0.22394, 0.4920};
    inline static constexpr std::array<double, N> molarMass{
        0.018015268, 0.0440098, 0.14228168};

    static constexpr double eosOmegaA = 0.45724;
    static constexpr double eosOmegaB = 0.07780;
    static constexpr double eosU = 2.4142135623730951;
    static constexpr double eosW = -0.4142135623730951;
    static constexpr double temperature = BenchmarkCommon::temperature;
    inline static constexpr std::array<double, 3> surfaceDensity{730.0, 1.8, 985.4040020947351};
    inline static constexpr std::array<double, 3> viscosity{2.4e-4, 1.86e-5, 4.69091e-4};
    static constexpr double waterViscosity = 4.69091e-4;
    static constexpr double waterFormationVolumeFactor = 1.0;

    static constexpr bool useIapwsGarciaAqueousVolume = true;
    static constexpr double aqueousVolumeWaterMolarMass = 0.018015268;
    static constexpr double aqueousVolumeMaximumUnsupportedMoleFraction = 1.0e-4;
    static constexpr bool useMcBrideWrightAqueousViscosity = true;
    static constexpr double aqueousViscosityMaximumUnsupportedMoleFraction = 1.0e-4;
};

struct PrFluid : CommonFluid
{
    static constexpr auto thermodynamicModel = MPMC::CubicThermodynamicModel::PengRobinson;
    static constexpr int eosModelFlag = 5;
    inline static constexpr std::array<std::array<double, N>, N> binaryInteraction{{
        {{0.0, 0.1896, 0.5000}},
        {{0.1896, 0.0, 0.1141}},
        {{0.5000, 0.1141, 0.0}}
    }};
};

struct SwFluid : CommonFluid
{
    static constexpr auto thermodynamicModel = MPMC::CubicThermodynamicModel::SoreideWhitson;
    static constexpr int eosModelFlag = 1;
    static constexpr double soreideWhitsonSalinityMolality = 0.0;
    inline static constexpr std::array<std::array<double, N>, N> binaryInteraction{{
        {{0.0, 0.1896, 0.5000}},
        {{0.1896, 0.0, 0.1141}},
        {{0.5000, 0.1141, 0.0}}
    }};
    static double soreideWhitsonAqueousWaterBip(int component, double, double)
    {
        if (component == co2Component) return -0.06609;
        if (component == nc10Component) return -0.14229;
        return 0.0;
    }
};

struct CpaFluid : CommonFluid
{
    static constexpr auto thermodynamicModel = MPMC::CubicThermodynamicModel::CubicPlusAssociation;
    static constexpr auto cpaCubicPhysicalTerm = MPMC::CpaCubicPhysicalTerm::PengRobinson;
    static constexpr int eosModelFlag = 5;
    inline static constexpr std::array<std::array<double, N>, N> binaryInteraction{{
        {{0.0, 0.10566, 0.0}},
        {{0.10566, 0.0, 0.1141}},
        {{0.0, 0.1141, 0.0}}
    }};
    inline static constexpr std::array<double, N> cpaA0{
        0.15782,
        0.45724 * gasConstant * gasConstant * criticalTemperature[1] * criticalTemperature[1] / criticalPressure[1],
        0.45724 * gasConstant * gasConstant * criticalTemperature[2] * criticalTemperature[2] / criticalPressure[2]};
    inline static constexpr std::array<double, N> cpaB{
        1.4788e-5,
        0.07780 * gasConstant * criticalTemperature[1] / criticalPressure[1],
        0.07780 * gasConstant * criticalTemperature[2] / criticalPressure[2]};
    inline static constexpr std::array<double, N> cpaC1{
        0.6736,
        0.37464 + 1.54226 * acentricFactor[1] - 0.26992 * acentricFactor[1] * acentricFactor[1],
        0.379642 + 1.48503 * acentricFactor[2] - 0.164423 * acentricFactor[2] * acentricFactor[2]
            + 0.016666 * acentricFactor[2] * acentricFactor[2] * acentricFactor[2]};
    // CO2 has no donor, hence cannot self-associate; the positive values below
    // are its B2 cross-association parameters and satisfy the site invariant.
    inline static constexpr std::array<double, N> cpaAssociationEnergy{16123.0, 8061.5, 0.0};
    inline static constexpr std::array<double, N> cpaAssociationVolume{0.069662, 0.15182, 0.0};
    inline static constexpr std::array<int, N> cpaDonorSites{2, 0, 0};
    inline static constexpr std::array<int, N> cpaAcceptorSites{2, 1, 0};
    inline static constexpr std::array<std::array<double, N>, N> cpaCrossAssociationEnergy{{
        {{0.0, 8061.5, 0.0}}, {{0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}
    }};
    inline static constexpr std::array<std::array<double, N>, N> cpaCrossAssociationVolume{{
        {{0.0, 0.15182, 0.0}}, {{0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}
    }};
    static constexpr auto cpaRadialDistribution = MPMC::CpaRadialDistribution::Simplified;
};

struct PrFactoryConfig { using Fluid = PrFluid; };
struct SwFactoryConfig { using Fluid = SwFluid; };
struct CpaFactoryConfig { using Fluid = CpaFluid; };

struct InitialState
{
    static constexpr double pressure = BenchmarkCommon::initialPressure;
    static constexpr double temperature = BenchmarkCommon::temperature;
    inline static constexpr std::array<double, CommonFluid::N> overallComposition{0.20, 0.0, 0.80};
};

struct Dissolution
{
    static constexpr int component = CommonFluid::co2Component;
    static constexpr double waterMolarMass = 0.018015268;
    static constexpr double salinityMolality = 0.0;
    static constexpr double initialWaterCO2MoleFraction = 0.0;
};
struct Land { static constexpr double constant = 0.0; };
struct Adsorption
{
    static constexpr double rockDensity = 2650.0;
    inline static constexpr std::array<double, CommonFluid::N> thetaMax{};
    inline static constexpr std::array<double, CommonFluid::N> coefficient{};
    static constexpr double standardPressure = 101325.0;
    static constexpr double standardTemperature = 288.15;
};
struct Numerics
{
    static constexpr double fugacityScalingFactor = 1.0;
    static constexpr bool useVariableBounds = false;
    static constexpr bool enableSnesStagnationGuard = true;
    static constexpr int snesStagnationMinimumIterations = 8;
    static constexpr int snesStagnationWindow = 5;
    static constexpr double snesStagnationRelativeImprovement = 1.0e-4;
};
struct Time
{
    static constexpr int numberOfSteps = 50;
    static constexpr double dtDays = 36.525;
    static constexpr bool adaptive = true;
    static constexpr double minimumDtDays = 1.0e-5;
    static constexpr double cutFactor = 0.5;
    static constexpr double growthFactor = 1.25;
    static constexpr double difficultShrinkFactor = 0.8;
    static constexpr int easyNonlinearIterations = 6;
    static constexpr int difficultNonlinearIterations = 14;
    static constexpr int maximumRetries = 16;
    static constexpr int maximumWellControlIterations = 8;
    static constexpr bool printAdaptiveSteps = true;
};
struct Output
{
    inline static constexpr const char *directory = "./results";
    static constexpr std::size_t every = 1;
    static constexpr bool printWells = true;
    static constexpr bool printWellPhaseDetails = true;
    static constexpr bool writeWellHistory = true;
    static constexpr bool printInventory = true;
    static constexpr bool enableComponentMassBalance = true;
    static constexpr bool printComponentMassBalance = true;
    static constexpr bool writeComponentMassBalance = true;
    static constexpr bool writeInventoryHistory = true;
    static constexpr bool writeMassTotals = true;
    static constexpr bool writeSolutionSnapshots = false;
    static constexpr bool writePhaseStateSnapshots = false;
    static constexpr bool printNewtonIterations = true;
    static constexpr bool writeSolverHistory = true;
    static constexpr bool printWellControlSwitches = true;
    static constexpr bool writeWellControlSwitchHistory = true;
    static constexpr bool printReservoirDiagnostics = true;
    static constexpr bool writeReservoirDiagnostics = true;
    static constexpr bool printFinalSummary = true;
    static constexpr bool writeFinalSummary = true;
};

struct Config final
{
    inline static constexpr const char *name = CaseConfig::name;
    using Model = CaseConfig::Model;
    using Grid = CaseConfig::Grid;
    using Fluid = CaseConfig::PrFluid;
    using InitialState = CaseConfig::InitialState;
    using Dissolution = CaseConfig::Dissolution;
    using Land = CaseConfig::Land;
    using Adsorption = CaseConfig::Adsorption;
    using Numerics = CaseConfig::Numerics;
    using Time = CaseConfig::Time;
    using Output = CaseConfig::Output;
    using Rock = CaseConfig::Rock;
};

} // namespace CaseConfig
