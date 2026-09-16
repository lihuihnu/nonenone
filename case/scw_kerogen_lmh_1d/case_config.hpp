#pragma once

#include "../scw_kerogen_common/benchmark_common.hpp"
#include "../scw_kerogen_common/calibrated_binary_properties.hpp"

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
        "H2O", "Light_nC4", "Middle_nC10", "Heavy_squalane"};

    inline static constexpr std::array<double, N> criticalTemperature{
        647.096, 425.125, 617.7, 855.9};
    inline static constexpr std::array<double, N> criticalPressure{
        22.064e6, 3.7960e6, 2.110e6, 0.7164e6};
    inline static constexpr std::array<double, N> criticalVolume{
        7.49587808839913e-5, 2.5492e-4, 6.0300e-4,
        3.05355324646748e-3};
    inline static constexpr std::array<double, N> acentricFactor{
        0.3443, 0.20081, 0.4884, 1.255};
    inline static constexpr std::array<double, N> molarMass{
        0.01801528, 0.0581222, 0.142286, 0.4228};

    // Only H2O-squalane is tied to the existing 653.2 K experimental
    // parameterization.  H2O-nC4/H2O-nC10 are screening values and all
    // hydrocarbon-hydrocarbon BIPs are zero; this case is a mechanism study,
    // not a characterized kerogen-oil PVT model.
    inline static constexpr std::array<std::array<double, N>, N>
        binaryInteraction{{
            {{0.0,    0.5091, 0.2618373654302397, 0.053233659543451335}},
            {{0.5091, 0.0,    0.0,                0.0}},
            {{0.2618373654302397, 0.0, 0.0,        0.0}},
            {{0.053233659543451335, 0.0, 0.0,      0.0}}
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
        if (other == lightComponent)
            return ScwKerogenCalibration::nonAqueousNc4Kij(temperatureK);
        if (other == middleComponent)
            return ScwKerogenCalibration::extrapolatedPrNc10Kij(temperatureK);
        return ScwKerogenCalibration::prSqualaneKij(temperatureK);
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
        800.0, 120.0, 360.0};
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
    // Preserve the original light:middle:heavy hydrocarbon ratio.  z_H2O=0.25
    // is the common PR/SW/CPA one-phase hydrocarbon-rich baseline; it is not an
    // artificial finite water saturation because PTz flash returns S_w=0.
    static constexpr double initialH2OMoleFraction = 0.25;
    inline static constexpr std::array<double, Fluid::N> overallComposition{
        initialH2OMoleFraction,
        (1.0 - initialH2OMoleFraction)
            * 0.1041666666666667,
        (1.0 - initialH2OMoleFraction)
            * 0.2604166666666667,
        (1.0 - initialH2OMoleFraction)
            * 0.6354166666666666};
};

struct Dissolution
{
    static constexpr int component = Fluid::lightComponent;
    static constexpr double waterMolarMass = 0.01801528;
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
