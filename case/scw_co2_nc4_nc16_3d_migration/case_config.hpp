/**
 * @file case_config.hpp
 * @brief 等温三维超临界水 + PR(CO2-nC4-nC16) 运移算例配置。
 */
#pragma once

#include "benchmark_common.hpp"

#include <indices/model_config.hpp>
#include <natural/thermo/thermodynamic_model.hpp>

#include <array>
#include <cstddef>

namespace CaseConfig
{

inline constexpr char name[] = "scw_co2_nc4_nc16_3d_migration";

struct Model
{
    // H2O is the separately conserved mobile SCW phase; the compositional
    // hydrocarbon subsystem contains CO2, light nC4 and heavy nC16.
    static constexpr int numberOfComponents = 3;
    static constexpr bool hasWater = true;
    static constexpr bool hasWells = true;
    static constexpr bool enableDissolution = false;
    static constexpr bool enableAdsorption = false;
    static constexpr bool enableLandTrapping = false;
    static constexpr auto phaseBehavior =
        MPMC::PhaseBehaviorModel::LegacyOilGasWithIndependentWater;
};

using Grid = ScwMigration3D::Grid;
using Rock = ScwMigration3D::Rock;

struct Fluid
{
    static constexpr std::size_t N = Model::numberOfComponents;
    static constexpr int co2Component = 0;
    static constexpr int nc4Component = 1;
    static constexpr int nc16Component = 2;
    inline static constexpr std::array<const char *, N> componentNames{
        "CO2", "nC4", "nC16"};

    inline static constexpr std::array<double, N> criticalTemperature{
        304.1282, 425.125, 723.0};
    inline static constexpr std::array<double, N> criticalPressure{
        7.3773e6, 3.7960e6, 1.410e6};
    inline static constexpr std::array<double, N> criticalVolume{
        9.4118e-5, 2.5492e-4, 9.0e-4};
    inline static constexpr std::array<double, N> acentricFactor{
        0.22394, 0.20081, 0.742};
    inline static constexpr std::array<double, N> molarMass{
        0.0440098, 0.0581222, 0.2264412};
    inline static constexpr std::array<std::array<double, N>, N>
        binaryInteraction{{
            {{0.0, 0.1277, 0.0900}},
            {{0.1277, 0.0, 0.0}},
            {{0.0900, 0.0, 0.0}}
        }};

    static constexpr auto thermodynamicModel =
        MPMC::CubicThermodynamicModel::PengRobinson;
    static constexpr int eosModelFlag = 5;
    static constexpr double eosOmegaA = 0.45724;
    static constexpr double eosOmegaB = 0.07780;
    static constexpr double eosU = 2.4142135623730951;
    static constexpr double eosW = -0.4142135623730951;
    static constexpr double temperature = ScwMigration3D::temperature;

    inline static constexpr std::array<double, 3> surfaceDensity{
        720.0, 1.8, 998.2};
    inline static constexpr std::array<double, 3> viscosity{
        2.8e-4, 3.5e-5, 7.5e-5};
    static constexpr double waterViscosity = 7.5e-5;
    static constexpr double waterFormationVolumeFactor = 1.0;
};

struct PrFactoryConfig { using Fluid = CaseConfig::Fluid; };

struct InitialState
{
    static constexpr double pressure = ScwMigration3D::initialPressure;
    static constexpr double temperature = ScwMigration3D::temperature;
    static constexpr bool preserveReferenceValues = false;
    static constexpr double oilSaturation =
        ScwMigration3D::initialOilSaturation;
    static constexpr double gasSaturation =
        ScwMigration3D::initialGasSaturation;
    static constexpr double waterSaturation =
        ScwMigration3D::initialWaterSaturation;

    // Heavy oil with a mobile/light nC4 fraction and a small dissolved CO2 seed.
    inline static constexpr std::array<double, Fluid::N> oilComposition{
        0.02, 0.23, 0.75};
    inline static constexpr std::array<double, Fluid::N> gasComposition{
        0.995, 0.005, 0.0};
};

struct Dissolution
{
    static constexpr int component = Fluid::co2Component;
    static constexpr double waterMolarMass = 0.018015268;
    static constexpr double salinityMolality = 0.0;
    static constexpr double initialWaterCO2MoleFraction = 0.0;
};
struct Land { static constexpr double constant = 0.0; };
struct Adsorption
{
    static constexpr double rockDensity = 2650.0;
    inline static constexpr std::array<double, Fluid::N> thetaMax{};
    inline static constexpr std::array<double, Fluid::N> coefficient{};
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
    // 100 common targets span ten years = 0.25 nominal injected PV.
    static constexpr int numberOfSteps = 100;
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
    static constexpr std::size_t every = 5;
    static constexpr bool printWells = true;
    static constexpr bool printWellPhaseDetails = true;
    static constexpr bool writeWellHistory = true;
    static constexpr bool printInventory = true;
    static constexpr bool enableComponentMassBalance = true;
    static constexpr bool printComponentMassBalance = true;
    static constexpr bool writeComponentMassBalance = true;
    static constexpr bool writeInventoryHistory = true;
    static constexpr bool writeMassTotals = true;
    static constexpr bool writeSolutionSnapshots = true;
    static constexpr bool writePhaseStateSnapshots = true;
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
    using Fluid = CaseConfig::Fluid;
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
