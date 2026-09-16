/**
 * @file case_config.hpp
 * @brief 传统基线 H2O-CO2-nC10 二维算例的模型与物性配置。
 */
#pragma once

#include "../h2o_co2_nc10_2d_benchmark/benchmark_common.hpp"

#include <indices/model_config.hpp>
#include <natural/thermo/thermodynamic_model.hpp>

#include <array>
#include <cstddef>

namespace TraditionalConfig
{
#ifndef MPMC_H2O_CO2_NC10_TRADITIONAL_CASE_NAME
#define MPMC_H2O_CO2_NC10_TRADITIONAL_CASE_NAME "h2o_co2_nc10_2d_traditional"
#endif
inline constexpr char name[] = MPMC_H2O_CO2_NC10_TRADITIONAL_CASE_NAME;

struct Model
{
    static constexpr int numberOfComponents = 2;
    static constexpr bool hasWater = true;
    static constexpr bool hasWells = true;
    static constexpr bool enableDissolution = false;
    static constexpr bool enableAdsorption = false;
    static constexpr bool enableLandTrapping = false;
    static constexpr auto phaseBehavior = MPMC::PhaseBehaviorModel::LegacyOilGasWithIndependentWater;
};
using Grid = BenchmarkCommon::Grid;
using Rock = BenchmarkCommon::Rock;

struct Fluid
{
    static constexpr std::size_t N = Model::numberOfComponents;
    static constexpr int co2Component = 0;
    static constexpr int nc10Component = 1;
    inline static constexpr std::array<const char *, N> componentNames{"CO2", "nC10"};
    inline static constexpr std::array<double, N> criticalTemperature{304.1282, 617.70};
    inline static constexpr std::array<double, N> criticalPressure{7.3773e6, 2.103e6};
    inline static constexpr std::array<double, N> criticalVolume{9.4118e-5, 6.10e-4};
    inline static constexpr std::array<double, N> acentricFactor{0.22394, 0.4920};
    inline static constexpr std::array<double, N> molarMass{0.0440098, 0.14228168};
    inline static constexpr std::array<std::array<double, N>, N> binaryInteraction{{
        {{0.0, 0.1141}}, {{0.1141, 0.0}}
    }};
    static constexpr auto thermodynamicModel = MPMC::CubicThermodynamicModel::PengRobinson;
    static constexpr int eosModelFlag = 5;
    static constexpr double eosOmegaA = 0.45724;
    static constexpr double eosOmegaB = 0.07780;
    static constexpr double eosU = 2.4142135623730951;
    static constexpr double eosW = -0.4142135623730951;
    static constexpr double temperature = BenchmarkCommon::temperature;
    inline static constexpr std::array<double, 3> surfaceDensity{730.0, 1.8, 985.4040020947351};
    inline static constexpr std::array<double, 3> viscosity{2.4e-4, 1.86e-5, 4.69091e-4};
    static constexpr double waterViscosity = 4.69091e-4;
    static constexpr double waterFormationVolumeFactor = 1.0;
};
struct PrFactoryConfig { using Fluid = TraditionalConfig::Fluid; };

struct InitialState
{
    static constexpr double pressure = BenchmarkCommon::initialPressure;
    static constexpr double temperature = BenchmarkCommon::temperature;
    static constexpr bool preserveReferenceValues = false;
    static constexpr double oilSaturation = BenchmarkCommon::targetOilSaturation;
    static constexpr double gasSaturation = 0.0;
    static constexpr double waterSaturation = BenchmarkCommon::targetWaterSaturation;
    inline static constexpr std::array<double, Fluid::N> oilComposition{0.0, 1.0};
    inline static constexpr std::array<double, Fluid::N> gasComposition{0.0, 1.0};
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
    inline static constexpr const char *directory = "./results/traditional";
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
    inline static constexpr const char *name = TraditionalConfig::name;
    using Model = TraditionalConfig::Model;
    using Grid = TraditionalConfig::Grid;
    using Fluid = TraditionalConfig::Fluid;
    using InitialState = TraditionalConfig::InitialState;
    using Dissolution = TraditionalConfig::Dissolution;
    using Land = TraditionalConfig::Land;
    using Adsorption = TraditionalConfig::Adsorption;
    using Numerics = TraditionalConfig::Numerics;
    using Time = TraditionalConfig::Time;
    using Output = TraditionalConfig::Output;
    using Rock = TraditionalConfig::Rock;
};
} // namespace TraditionalConfig
