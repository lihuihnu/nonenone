#pragma once

#include "../five_component_eos_tuned_compare/case_config.hpp"

#include <natural/thermo/thermodynamic_model.hpp>

#include <array>
#include <cstddef>

/**
 * @file case_config.hpp
 * @brief 五种物理组分的传统油气组分模型加独立水相对照算例。
 *
 * 传统模型的 EOS 只包含 CO2/CH4/C2H6/nC4H10 四个可挥发组分；H2O
 * 由独立水守恒方程表示，不进入油气 Flash，也不能在油、气、水相之间分配。
 * 该传统对照固定使用 PR；网格、岩石、井、温压、干流体比例和时间表复用
 * fully-compositional 算例。
 */
namespace LegacyCaseConfig
{

inline constexpr char name[] = "five_component_eos_tuned_legacy_compare";
inline constexpr double bar = CaseConfig::bar;
inline constexpr double mD = CaseConfig::mD;
inline constexpr double secondsPerDay = CaseConfig::secondsPerDay;
inline constexpr double gasConstant = CaseConfig::gasConstant;

struct Model
{
    // Four EOS components plus the independent H2O conservation equation =
    // five physical components in the reservoir model.
    static constexpr int numberOfComponents = 4;
    static constexpr bool hasWater = true;
    static constexpr bool hasWells = true;
    static constexpr bool enableDissolution = false;
    static constexpr bool enableAdsorption = false;
    static constexpr bool enableLandTrapping = false;
    static constexpr auto phaseBehavior =
        MPMC::PhaseBehaviorModel::LegacyOilGasWithIndependentWater;
};

using Grid = CaseConfig::Grid;
using Rock = CaseConfig::Rock;

struct CommonFluid
{
    static constexpr std::size_t N = Model::numberOfComponents;
    static constexpr int co2Component = 0;
    static constexpr int methaneComponent = 1;
    static constexpr int ethaneComponent = 2;
    static constexpr int nButaneComponent = 3;

    inline static constexpr std::array<const char *, N> componentNames{
        "CO2", "CH4", "C2H6", "nC4H10"};
    inline static constexpr std::array<double, N> criticalTemperature{
        304.128200003, 190.564002651, 305.322, 425.125};
    inline static constexpr std::array<double, N> criticalPressure{
        7.377298373e6, 4.599200474e6, 4.872200000e6, 3.796000017e6};
    inline static constexpr std::array<double, N> criticalVolume{
        9.411848339e-5, 9.862771707e-5, 1.458387816e-4, 2.549219298e-4};
    inline static constexpr std::array<double, N> acentricFactor{
        0.22394, 0.01142, 0.099, 0.200810094644};
    inline static constexpr std::array<double, N> molarMass{
        0.0440098, 0.0160428, 0.03006904, 0.0581222};

    static constexpr double eosOmegaA = 0.45724;
    static constexpr double eosOmegaB = 0.07780;
    static constexpr double eosU = 2.414213562373095;
    static constexpr double eosW = -0.414213562373095;
    static constexpr double temperature = 305.0;

    inline static constexpr std::array<double, 3> surfaceDensity{
        620.0, 1.80, 998.0};
    inline static constexpr std::array<double, 3> viscosity{
        3.0e-4, 1.5e-5, 8.0e-4};
    static constexpr double waterViscosity = 8.0e-4;
    static constexpr double waterFormationVolumeFactor = 1.0;
};

struct PrFluid : CommonFluid
{
    static constexpr auto thermodynamicModel =
        MPMC::CubicThermodynamicModel::PengRobinson;
    static constexpr int eosModelFlag = 5;
    inline static constexpr std::array<std::array<double, N>, N> binaryInteraction{{
        {{0.0, 0.1100076406250000, 0.1743251250000000, 0.1277}},
        {{0.1100076406250000, 0.0, 0.0022815210880000, 0.0108823211788800}},
        {{0.1743251250000000, 0.0022815210880000, 0.0, 0.0000}},
        {{0.1277, 0.0108823211788800, 0.0000, 0.0}}
    }};
};

struct PrFactoryConfig { using Fluid = PrFluid; };

struct InitialState
{
    static constexpr double pressure = 60.0 * bar;
    static constexpr double temperature = CommonFluid::temperature;
    static constexpr double overallWaterMoleFraction = 0.30;
    static constexpr double waterMolarMass = 0.018015268;

    // Original wet z=[0.30,0.10,0.15,0.15,0.30], with H2O removed and the
    // remaining 0.70 dry fraction normalized to one.
    inline static constexpr std::array<double, CommonFluid::N> dryOverallComposition{
        1.0 / 7.0, 3.0 / 14.0, 3.0 / 14.0, 3.0 / 7.0};

    // The generic legacy-case validator requires a finite reference state.
    // Runtime initialization replaces these reference values with the dry PR
    // flash and an Sw that exactly reconstructs overall z_H2O=0.30 from
    // dry-phase and independent-water molar densities.
    static constexpr bool preserveReferenceValues = false;
    static constexpr double waterSaturation = 0.073297;
    static constexpr double oilSaturation = 0.593791;
    static constexpr double gasSaturation =
        1.0 - waterSaturation - oilSaturation;
    inline static constexpr auto oilComposition = dryOverallComposition;
    inline static constexpr auto gasComposition = dryOverallComposition;

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

using Numerics = CaseConfig::Numerics;
using Time = CaseConfig::Time;
using Output = CaseConfig::Output;

struct Config final
{
    inline static constexpr const char *name = LegacyCaseConfig::name;
    using Model = LegacyCaseConfig::Model;
    using Grid = LegacyCaseConfig::Grid;
    using Fluid = LegacyCaseConfig::PrFluid;
    using InitialState = LegacyCaseConfig::InitialState;
    using Dissolution = LegacyCaseConfig::Dissolution;
    using Land = LegacyCaseConfig::Land;
    using Adsorption = LegacyCaseConfig::Adsorption;
    using Numerics = LegacyCaseConfig::Numerics;
    using Time = LegacyCaseConfig::Time;
    using Output = LegacyCaseConfig::Output;
    using Rock = LegacyCaseConfig::Rock;
};

} // namespace LegacyCaseConfig
