#pragma once

#include "../scw_kerogen_common/benchmark_common.hpp"
#include "../scw_kerogen_common/bsb_reference_properties.hpp"

#include <indices/model_config.hpp>
#include <natural/thermo/thermodynamic_model.hpp>

#include <array>
#include <cstddef>

/**
 * @file case_config.hpp
 * @brief 第二阶段水与轻中重干酪根裂解产物拟组分配置。
 */
namespace CaseConfig
{

inline constexpr char name[] = "scw_kerogen_lmh_1d";

struct Model
{
    static constexpr int numberOfComponents = 4;
    static constexpr bool hasWater = true;
    static constexpr bool hasWells = true;
    static constexpr bool enableDissolution = false;
    static constexpr bool enableAdsorption = false;
    static constexpr bool enableLandTrapping = false;
    static constexpr auto phaseBehavior =
        MPMC::PhaseBehaviorModel::FullyCompositionalThreePhase;
};

using Grid = ScwKerogen1D::Grid;
using Rock = ScwKerogen1D::Rock;
using Numerics = ScwKerogen1D::Numerics;
using Time = ScwKerogen1D::Time;
using Output = ScwKerogen1D::Output;

struct Fluid
{
    static constexpr std::size_t N = Model::numberOfComponents;
    static constexpr int waterComponent = 0;
    static constexpr int lightComponent = 1;
    static constexpr int middleComponent = 2;
    static constexpr int heavyComponent = 3;
    inline static constexpr std::array<const char *, N> componentNames{
        "H2O", "Light_BSB_C1_C6", "Middle_BSB_C7_C15", "Heavy_BSB_C16_C27"};

    inline static constexpr std::array<double, N> criticalTemperature{
        BsbReference::water.criticalTemperatureK,
        BsbReference::light.criticalTemperatureK,
        BsbReference::middle.criticalTemperatureK,
        BsbReference::heavy.criticalTemperatureK};
    inline static constexpr std::array<double, N> criticalPressure{
        BsbReference::water.criticalPressurePa,
        BsbReference::light.criticalPressurePa,
        BsbReference::middle.criticalPressurePa,
        BsbReference::heavy.criticalPressurePa};
    inline static constexpr std::array<double, N> criticalVolume{
        BsbReference::water.criticalVolumeM3PerMol,
        BsbReference::light.criticalVolumeM3PerMol,
        BsbReference::middle.criticalVolumeM3PerMol,
        BsbReference::heavy.criticalVolumeM3PerMol};
    inline static constexpr std::array<double, N> acentricFactor{
        BsbReference::water.acentricFactor,
        BsbReference::light.acentricFactor,
        BsbReference::middle.acentricFactor,
        BsbReference::heavy.acentricFactor};
    inline static constexpr std::array<double, N> molarMass{
        BsbReference::water.molarMassKgPerMol,
        BsbReference::light.molarMassKgPerMol,
        BsbReference::middle.molarMassKgPerMol,
        BsbReference::heavy.molarMassKgPerMol};

    // Match the H2O-CO2-BSB PR baseline: all water-hydrocarbon BIPs are 0.5
    // and hydrocarbon-hydrocarbon BIPs are zero.
    inline static constexpr std::array<std::array<double, N>, N>
        binaryInteraction{{
            {{0.0, BsbReference::waterHydrocarbonKij,
                   BsbReference::waterHydrocarbonKij,
                   BsbReference::waterHydrocarbonKij}},
            {{BsbReference::waterHydrocarbonKij, 0.0, 0.0, 0.0}},
            {{BsbReference::waterHydrocarbonKij, 0.0, 0.0, 0.0}},
            {{BsbReference::waterHydrocarbonKij, 0.0, 0.0, 0.0}}
        }};

    static double cubicBinaryInteractionCoefficient(
        int i, int j, double temperatureK)
    {
        if (i == j)
            return 0.0;
        const int other = i == waterComponent ? j
            : (j == waterComponent ? i : -1);
        if (other < 0)
            return 0.0;
        (void)temperatureK;
        return BsbReference::waterHydrocarbonKij;
    }

    static constexpr auto thermodynamicModel =
        MPMC::CubicThermodynamicModel::PengRobinson;
    static constexpr int eosModelFlag = 1;
    static constexpr double eosOmegaA = 0.45724;
    static constexpr double eosOmegaB = 0.07780;
    static constexpr double eosU = 2.4142135623730951;
    static constexpr double eosW = -0.4142135623730951;
    static constexpr double temperature = ScwKerogen1D::temperature;

    inline static constexpr std::array<double, 3> surfaceDensity{
        750.0, 1.8, 985.4040020947351};
    inline static constexpr std::array<double, 3> viscosity{
        8.0e-4, 4.0e-5, 7.0e-5};
    static constexpr double waterViscosity = 7.0e-5;
    static constexpr double waterFormationVolumeFactor = 1.0;
    static constexpr bool useIapws2008AqueousViscosity = true;
    static constexpr double aqueousViscosityMaximumSoluteMoleFraction = 0.02;
};

struct InitialState
{
    static constexpr double pressure = ScwKerogen1D::initialPressure;
    static constexpr double temperature = ScwKerogen1D::temperature;
    // Preserve the 20/80 water/hydrocarbon split.  Within the hydrocarbon cut,
    // use the BSB L/M/H mole fractions renormalized after excluding XH.
    static constexpr double initialH2OMoleFraction =
        BsbReference::initialWaterMoleFraction;
    inline static constexpr std::array<double, Fluid::N> overallComposition{
        initialH2OMoleFraction,
        (1.0 - initialH2OMoleFraction)
            * BsbReference::lmhMoleFraction[0],
        (1.0 - initialH2OMoleFraction)
            * BsbReference::lmhMoleFraction[1],
        (1.0 - initialH2OMoleFraction)
            * BsbReference::lmhMoleFraction[2]};
};

struct Dissolution
{
    static constexpr int component = Fluid::lightComponent;
    static constexpr double waterMolarMass = BsbReference::water.molarMassKgPerMol;
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
